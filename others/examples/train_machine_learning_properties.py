"""Train and reload an ML property model with editable Python settings.

This version adds scikit-learn/optional boosting models and three feature sources.
Run normally after editing INPUT_FILE, or use --demo for a synthetic smoke run.
"""

from __future__ import annotations

import argparse
import csv
import random
import tempfile
from pathlib import Path

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    DataMapping,
    predict_custom_properties,
    train_custom_property_model,
)

INPUT_FILE = Path("my_labels.csv")
OUTPUT_DIR = Path("runs/my-ml-property")
TASK_TYPE = "regression"  # classification / regression
MODEL = "random_forest"  # See available_machine_learning_models().
MODEL_PARAMS = {"n_estimators": 100, "max_depth": 8}
FEATURE_SOURCE = "sequence"  # sequence / embedding / numeric
EMBEDDING_MODEL = "dnabert2"  # Used only with feature_source="embedding".
FEATURE_FIELDS = ("length", "gc_fraction")  # Used only for numeric input.


def train(source: Path, output_dir: Path) -> None:
    config = CustomPropertyTrainingConfig(
        task_type=TASK_TYPE,
        model=MODEL,
        finetune_method="machine_learning",
        model_params=MODEL_PARAMS,
        feature_source=FEATURE_SOURCE,
        embedding_model=EMBEDDING_MODEL,
        mapping=DataMapping(
            sequence_field=None if FEATURE_SOURCE == "numeric" else "sequence",
            feature_fields=FEATURE_FIELDS if FEATURE_SOURCE == "numeric" else (),
            label_fields=("target",),
            id_field="id",
        ),
        kmer_k=3,
        standardize_features=True,
        workers=6,
        seed=42,
        output_dir=output_dir,
        # For an embedding checkpoint that needs its own code, explicitly set
        # allow_remote_code=True after checking the checkpoint's source.
    )
    result = train_custom_property_model(source, config=config)
    print("测试集指标:", result.metrics["test"])
    print("模型保存至:", result.output_dir)
    prediction = predict_custom_properties(
        source,
        model_dir=result.output_dir,
        trust_model=True,  # This artifact was just trained by this script.
    )
    print("重新加载后的前三条预测:", prediction.predictions[:3])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", action="store_true", help="Use synthetic demonstration data.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    if not args.demo:
        train(INPUT_FILE, args.output_dir)
        return
    generator = random.Random(42)
    with tempfile.TemporaryDirectory(prefix="dnakit-ml-example-") as directory:
        source = Path(directory) / "demo.csv"
        with source.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("id", "sequence", "length", "gc_fraction", "target"))
            for index in range(120):
                sequence = "".join(generator.choices("ACGT", k=generator.randint(40, 160)))
                gc = (sequence.count("G") + sequence.count("C")) / len(sequence)
                target = gc * 3 + len(sequence) / 100
                if TASK_TYPE == "classification":
                    target = int(gc >= 0.5)
                writer.writerow((f"demo-{index}", sequence, len(sequence), gc, target))
        print("演示数据为人工生成, 仅用于验证训练与加载流程。")
        train(source, args.output_dir)


if __name__ == "__main__":
    main()
