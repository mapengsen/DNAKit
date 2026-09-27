# DNAKit examples

本目录保存可复现示例的输入和预期工件。

当前包含人工可核验的固定演示夹具和一个严格工作流：

- `fixed_demo.fasta`：三条短 DNA 序列；
- `fixed_demo_expected.json`：固定预期结果。
- `advanced_workflow.yml`：`dnakit-workflow-v1` 固定输入端到端配置。
- `train_custom_properties.py`：只需编辑 Python 配置的自定义分类/回归训练示例，包含 workers、每设备批大小、AdamW 参数和可复现划分；修改输入路径后可直接在 IDE 中运行。
- `train_machine_learning_properties.py`：机器学习分类/回归，支持序列特征、DNA embedding、数值列及可选提升树库。修改顶部配置后运行；`python others/examples/train_machine_learning_properties.py --demo` 可用人工数据验证训练和重新加载。
- `custom_llm.py`：在任意 Python 流程中使用自定义提示词处理文本/JSON，取得结果后继续本地处理。修改 `PROVIDER` 可选择 DeepSeek、Gemini、Qwen、Kimi、OpenRouter 或 OpenAI，`BASE_URL` 可覆盖地域地址；示例统一读取 `DNAKIT_LLM_API_KEY` 与 `DNAKIT_LLM_MODEL`。运行会发送人工数据并产生一次 API 请求，不自动读取文件。
- `custom_llm_file.py`：通过 `file_path` 指定本地文本文件，先本地预览，再发送一次模型请求。默认输入为仓库 `LICENSE`；可修改 `SOURCE_FILE`、编码及文件大小限制，密钥和模型名通过环境变量配置。
- `run_agent_with_memory.py`：使用 `--cache-dir` 保存本地 JSON 检查点，`--reuse-results` 恢复已完成步骤，`--rag-dir` 指定资料目录并先检索再交给模型。读取 `DNAKIT_LLM_MODEL`、`DNAKIT_LLM_API_KEY`，默认显示进度；运行 `python others/examples/run_agent_with_memory.py --help` 查看参数。

文档中的 CLI 和 Python 工作流直接使用这些固定输入；预期工件用于人工核验，不会被示例命令覆盖。可执行命令见 `docs/examples/index.md`，端到端 Notebook 见 `others/notebooks/00_skeleton_check.ipynb` 和 `others/notebooks/01_advanced_workflow.ipynb`。

- `bioinformatics_api.py`：原生 DNA 编码、MMD、性质命中率、同义密码子约束优化、建树、黏端连接和 Python DAG 的离线演示。先安装 `.[bio-analysis]`，再执行 `python others/examples/bioinformatics_api.py`；不调用外部程序或网络服务。
