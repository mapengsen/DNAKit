"""通过 Python 配置自定义性质训练。无需命令行参数。

本版示例: 支持训练参数与六个大模型服务商的 Python 配置。大模型辅助默认关闭。
先安装 dnakit[training]。修改下方文件与配置后在 IDE 中运行本脚本。
"""

from __future__ import annotations

from pathlib import Path

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    LLMSchemaConfig,
    train_custom_property_model,
)

# 路径相对于运行时工作目录。也可以填写绝对路径。
INPUT_FILE = Path("my_labels.csv")  # 也支持 xlsx、jsonl 等内置读取器支持的文件。
OUTPUT_DIR = Path("runs/my-property")  # 使用新的输出目录以避免覆盖已有结果。

# 如需大模型分析文件结构: 将 None 改为下面任意一个 provider 并填写模型名。
# deepseek / gemini / qwen / kimi / openrouter / openai / custom
# 密钥可在 PyCharm 的运行配置中设置环境变量。不要把真实密钥提交到代码仓库。
# 默认环境变量: DEEPSEEK_API_KEY / GEMINI_API_KEY / DASHSCOPE_API_KEY /
# MOONSHOT_API_KEY / OPENROUTER_API_KEY / OPENAI_API_KEY / LLM_API_KEY
LLM_PROVIDER: str | None = None
LLM_MODEL: str | None = None  # 填写服务商控制台中可用的模型 ID。
LLM_BASE_URL: str | None = None  # 使用预设地址。custom 必填。可覆盖地域或网关地址。
LLM_EXTRA_BODY: dict[str, object] = {}  # 例如 Qwen 可配置 {"enable_thinking": False}。


def main() -> None:
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(f"请先将 INPUT_FILE 改为已有的数据文件: {INPUT_FILE.resolve()}")

    llm = None
    if LLM_PROVIDER is not None:
        if not LLM_MODEL:
            raise ValueError("启用大模型辅助前请填写 LLM_MODEL。")
        llm = LLMSchemaConfig(
            provider=LLM_PROVIDER,
            model=LLM_MODEL,
            base_url=LLM_BASE_URL,
            include_sample_values=False,  # 此示例只发送字段/类型结构。不发送样本值。
            extra_body=LLM_EXTRA_BODY,
        )

    config = CustomPropertyTrainingConfig(
        task_type="regression",  # 分类任务改为 "classification"。
        model="dnabert2",  # 可调用 available_custom_training_models() 查看支持的模型。
        finetune_method="full",  # 线性探测改为 "linear_probe"。
        output_dir=OUTPUT_DIR,
        workers=6,  # Notebook 或不支持多进程的环境中设为 0。
        batch_size=32,  # 单个选定设备的批大小。当前训练器不支持多卡 DDP。
        learning_rate=1e-4,
        beta1=0.9,
        beta2=0.95,
        weight_decay=0.1,
        seed=42,
        train_ratio=0.8,
        validation_ratio=0.1,
        test_ratio=0.1,
        epochs=10,
        device="auto",  # 也可指定 "cpu"、"cuda:0" 等。
        allow_remote_code=True,  # DNABERT-2 需要此项。只加载你信任的 checkpoint。
        llm=llm,
        # 已知字段时可传 mapping=DataMapping(...)。默认不使用外部大模型。
    )
    result = train_custom_property_model(INPUT_FILE, config=config)
    print("测试集指标:", result.metrics["test"])
    print("模型、指标及训练配置已保存至:", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
