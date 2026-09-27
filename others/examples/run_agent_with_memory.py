# ruff: noqa: RUF001, RUF002
"""用本地文件缓存恢复 Agent，并检索用户指定的资料目录。

本版新增：--optimization-mode repair，允许兼容性源码修复并生成可投稿材料。
模型名称和密钥从调用参数/环境变量读取，脚本不内置模型或凭据。
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from dnakit.llm import LLMConfig
from dnakit.tools import EmbeddingConfig, PythonExecutionConfig, RAGConfig, run_agent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True, help="任务或知识库问题；恢复时使用同一文本")
    parser.add_argument("--cache-dir", type=Path, default=Path("agent_cache"), help="本地结果缓存")
    parser.add_argument("--reuse-results", action="store_true", help="恢复同一任务已验证的历史结果")
    parser.add_argument("--rag-dir", type=Path, help="本地 UTF-8 文本资料目录")
    parser.add_argument("--rag-index-dir", type=Path, help="可选的 RAG 索引目录")
    parser.add_argument("--rag-mode", choices=("keyword", "semantic", "hybrid"), default="keyword")
    parser.add_argument("--model", default=os.environ.get("DNAKIT_LLM_MODEL"))
    parser.add_argument("--provider", default=os.environ.get("DNAKIT_LLM_PROVIDER", "deepseek"))
    parser.add_argument("--base-url", default=os.environ.get("DNAKIT_LLM_BASE_URL"))
    parser.add_argument("--max-steps", type=int, default=50, help="本次模型调用总预算，含开发审查")
    parser.add_argument("--allow-side-effects", action="store_true", help="允许任务所需写入等操作")
    parser.add_argument("--workspace", type=Path, help="显式启用 Python 执行并指定工作目录")
    parser.add_argument("--generate-tools", action="store_true", help="验证并增加主包 Python 工具")
    parser.add_argument(
        "--optimize-tools",
        action="store_true",
        help="启用工具优化；使用 --optimization-mode 选择范围",
    )
    parser.add_argument(
        "--optimization-mode",
        choices=("description_only", "repair"),
        default="description_only",
        help="description_only 优化说明；repair 允许补全或修复本地 DNAKit 实现",
    )
    parser.add_argument("--quiet", action="store_true", help="关闭进度输出")
    args = parser.parse_args()
    if not args.model:
        parser.error("请设置 DNAKIT_LLM_MODEL 或传入 --model")
    if args.workspace and not args.allow_side_effects:
        parser.error("--workspace 启用 Python 执行，需要同时指定 --allow-side-effects")
    if (args.generate_tools or args.optimize_tools) and not args.allow_side_effects:
        parser.error("工具生成/优化需要 --allow-side-effects 和可写的 DNAKit 源码检出")
    if args.optimization_mode != "description_only" and not args.optimize_tools:
        parser.error("--optimization-mode repair 需要 --optimize-tools")
    if not args.rag_dir and (args.rag_index_dir or args.rag_mode != "keyword"):
        parser.error("RAG 配置需要 --rag-dir")
    embedding = None
    if args.rag_dir and args.rag_mode != "keyword":
        model = os.environ.get("DNAKIT_EMBEDDING_MODEL")
        endpoint = os.environ.get("DNAKIT_EMBEDDING_BASE_URL")
        if not model or not endpoint:
            parser.error("语义/混合检索需要 DNAKIT_EMBEDDING_MODEL 和 DNAKIT_EMBEDDING_BASE_URL")
        embedding = EmbeddingConfig(model=model, base_url=endpoint)
    rag = (
        RAGConfig(
            directory=args.rag_dir,
            index_dir=args.rag_index_dir,
            mode=args.rag_mode,
            embedding_config=embedding,
        )
        if args.rag_dir
        else None
    )
    result = run_agent(
        args.request,
        config=LLMConfig(
            model=args.model,
            provider=args.provider,
            base_url=args.base_url,
            api_key_env="DNAKIT_LLM_API_KEY",
        ),
        cache_dir=args.cache_dir,
        reuse_results=args.reuse_results,
        rag=rag,
        max_steps=args.max_steps,
        allow_side_effects=args.allow_side_effects,
        generate_tools=args.generate_tools,
        optimize_tools=args.optimize_tools,
        optimization_mode=args.optimization_mode,
        code_execution=PythonExecutionConfig(args.workspace) if args.workspace else None,
        verbose=not args.quiet,
    )
    print(result.text)
    print(f"状态：{result.status}；复用成功步骤：{result.reused_steps}")
    print(f"缓存目录：{result.cache_dir}")
    for item in result.knowledge_context:
        print(f"资料来源：{item['source']}，片段：{item['chunk_id']}")


if __name__ == "__main__":
    main()
