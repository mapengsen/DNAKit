"""Actual promoter rendering and geometry checks for pyGenomeViz 0.4 and 1.x."""
from pathlib import Path
from importlib.metadata import version
import pytest

pytest.importorskip("pygenomeviz")
from dnakit.visualization.pygenomeviz_p37_promoter_functions import draw_genome_promoters
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


@pytest.mark.parametrize("position,strand", [(0, 1), (500, -1), (1000, -1)])
def test_promoter_image_and_data(tmp_path, position, strand):
    promoters = [{"position": position, "strand": strand}]
    result = draw_genome_promoters(promoters, str(tmp_path / "p.svg"), segment_size=1000, verbose=False)
    assert Path(result["image_path"]).stat().st_size > 100
    assert result["promoters"][0]["position"] == position
    assert result["promoters"][0]["strand"] == strand
    assert result["promoters"][0]["arrow_end"] == [position + strand * 20, 1.0]
    if int(version("pygenomeviz").split(".")[0]) >= 1:
        from pygenomeviz import GenomeViz
        segment = GenomeViz().add_feature_track("direct", 1000).get_segment()
        segment.add_promoter(position, strand=strand)
        expected = segment.transform_ax_method2kws_list["annotate"][0]
        assert result["promoters"][0]["arrow_end"] == pytest.approx(expected["xy"])


def test_promoter_outside_segment_is_rejected(tmp_path):
    with pytest.raises(ConfigurationError, match="outside"):
        draw_genome_promoters([{"position": 1000}], str(tmp_path / "x.svg"), segment_size=1000, verbose=False)
