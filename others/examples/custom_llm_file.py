"""本版新增: 指定本地文本文件路径并交给大模型处理。

在 IDE 中设置 DNAKIT_LLM_API_KEY 和 DNAKIT_LLM_MODEL 后运行。
默认读取仓库 LICENSE。一次请求会发送完整文件内容和文件名。可能产生费用。
"""

from __future__ import annotations

import os
from pathlib import Path

from dnakit.llm import LLMConfig, call_llm, read_llm_file

PROVIDER = "deepseek"  # 也可选择 gemini、qwen、kimi、openrouter、openai。
SOURCE_FILE = Path(__file__).resolve().parents[2] / "LICENSE"
FILE_ENCODING = "utf-8-sig"
MAX_FILE_BYTES = 64_000


def main() -> None:
    # 本地预览不调用 API。可先确认这是你希望共享的文件内容。
    preview = read_llm_file(SOURCE_FILE, encoding=FILE_ENCODING, max_bytes=MAX_FILE_BYTES)
    print("Local preview:", preview[:200])

    result = call_llm(
        "用中文概括文件内容。只根据文件原文回答。",
        file_path=SOURCE_FILE,
        file_encoding=FILE_ENCODING,
        max_file_bytes=MAX_FILE_BYTES,
        config=LLMConfig(
            provider=PROVIDER,
            model=os.environ["DNAKIT_LLM_MODEL"],
            api_key_env="DNAKIT_LLM_API_KEY",
        ),
    )
    print(result.text)


if __name__ == "__main__":
    main()
