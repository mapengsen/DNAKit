"""在普通 Python 流程中调用自定义大模型提示词。

本版更新: 用 PROVIDER 切换六家服务商。保留自定义提示词与 JSON 后处理。
在 IDE 的环境变量中设置 DNAKIT_LLM_API_KEY 和 DNAKIT_LLM_MODEL 后运行。
本脚本会向所配置的 API 发送一次人工示例数据请求。可能产生费用。
"""

from __future__ import annotations

import os

from dnakit.llm import LLMConfig, call_llm

PROVIDER = "deepseek"  # 可选 gemini、qwen、kimi、openrouter、openai、custom。
BASE_URL: str | None = None  # 使用预设。地域/业务空间专属地址可在这里覆盖。


def main() -> None:
    config = LLMConfig(
        provider=PROVIDER,
        base_url=BASE_URL,
        model=os.environ["DNAKIT_LLM_MODEL"],  # 服务商提供的实际模型名称。
        api_key_env="DNAKIT_LLM_API_KEY",  # 也支持 api_key=...。不要将密钥提交到仓库。
        max_output_tokens=1_024,
        timeout_seconds=60,
        # extra_body={"enable_thinking": False},  # 仅为支持此参数的 Qwen 模型配置。
    )

    # 这里可以替换为任意 Python 处理步骤产生的文本、列表或字典。
    # 只传入你同意发给服务商的数据。此接口不会替你自动抽样。
    records = [
        {"id": "sample-a", "note": "The sample passed the initial quality check."},
        {"id": "sample-b", "note": "The measurement needs to be repeated."},
    ]
    result = call_llm(
        "将每条记录的 note 翻译为中文。返回 JSON 对象: "
        '{"items": [{"id": "保留原 id", "summary": "中文翻译"}]}。不要新增记录。',
        input_data=records,
        output_format="json",
        config=config,
    )

    # JSON 模式只保证语法与顶层对象类型。业务字段仍由你的代码验证。
    data = result.data
    items = data.get("items") if data is not None else None
    if not isinstance(items, list):
        raise ValueError("模型返回的 JSON 缺少 items 列表")
    summaries: dict[str, str] = {}
    for item in items:
        if (
            not isinstance(item, dict)
            or not isinstance(item.get("id"), str)
            or not isinstance(item.get("summary"), str)
            or item["id"] in summaries
        ):
            raise ValueError("模型返回了无效或重复的记录")
        summaries[item["id"]] = item["summary"]
    if set(summaries) != {record["id"] for record in records}:
        raise ValueError("模型返回的记录 ID 与输入不一致")

    # 继续你的 Python 流程。不执行模型返回的任何代码。
    for record in records:
        print(record["id"], summaries[record["id"]])
    print("Token usage:", dict(result.usage))


if __name__ == "__main__":
    main()
