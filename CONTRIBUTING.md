# 贡献指南

DNAKit 当前处于本地开发阶段，尚未建立公开 GitHub 协作流程。贡献应先保持小范围、可测试并与需求追踪矩阵对应。

## 开发环境

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
python -m pip install -e ".[dev,docs]"
```

环境必须保持 `PYTHONNOUSERSITE=1`，避免用户级包影响验证。

## 修改要求

- 支持 Python 3.10 及以上版本，使用类型注解。
- 公共行为需要 pytest 测试；边界行为需要明确错误和参数验证。
- 算法必须标记为 `native`、`adapter`、`reimplementation` 或 `novel`。
- adapter 必须记录后端名称、版本、参数、许可证提示和缺失后端行为。
- 随机流程必须记录 seed；外部数据库必须记录版本和校验信息。
- 不提交 NUPACK、Dashing、数据库快照、构建产物、缓存或其他不可分发工件。
- 不虚构尚未实现的 API、benchmark 或正确性结论。

## 本地检查

```bash
ruff check src tests
ruff format --check src tests
mypy
pytest
pytest --nbmake notebooks/00_skeleton_check.ipynb
mkdocs build --strict
python -m build
twine check dist/*
```

代码、测试、文档、追踪矩阵和变更记录应在同一次修改中保持一致。

## 提交内容

每项实现应说明：

1. 对应的需求 ID；
2. 公共 API 和 CLI 行为；
3. 输入边界、异常和数值容差；
4. 测试与对照工具；
5. 算法来源、引用、许可证和 provenance；
6. 用户可见的变更记录和示例。

## 发布限制

项目采用 MIT 许可证，但当前仍不创建正式发布、不上传 TestPyPI/PyPI、不推送 GitHub，也不部署网站。发布操作必须在项目所有者明确批准后单独执行。
