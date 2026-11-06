"""Offline comparisons against the gReLU metric and label-transform classes."""

import pytest

np = pytest.importorskip("numpy")
torch = pytest.importorskip("torch")
pytest.importorskip("grelu")
from dnakit.predictions import grelu_p36_functions as m


@pytest.mark.parametrize("values", [[[[-1.0, 0.0, 2.0]]], [[[0.0]]], [[[1.0, 3.0], [2.0, 4.0]]]])
def test_label_transform_matches_grelu(values):
    from grelu.transforms.label_transforms import LabelTransform

    expected = LabelTransform(min_clip=0, max_clip=3, transform_func="log1p").forward(
        np.asarray(values, dtype=float)
    )
    got = m.transform_labels(values, min_clip=0, max_clip=3, transform_func="log1p", verbose=False)
    np.testing.assert_allclose(got["labels"], expected)


@pytest.mark.parametrize("values", [[[[0.0, 1.0, 2.0]]], [[[0.0]]], [[[1.0, 3.0], [2.0, 4.0]]]])
def test_mse_matches_grelu(values):
    from grelu.lightning.metrics import MSE

    array = torch.tensor(values)
    metric = MSE(num_outputs=array.shape[1], average=False)
    metric.update(array, torch.zeros_like(array))
    result = m.compute_prediction_metric(
        "mse", values, torch.zeros_like(array).tolist(), verbose=False
    )
    np.testing.assert_allclose(result["per_task"], metric.compute().flatten().numpy())
