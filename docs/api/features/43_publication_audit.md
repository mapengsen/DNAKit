# 43) 320 项 DNA/基因工具的论文引用审计

这一页按照论文投稿要求审计 `public_dna_tools.md` 中的 320 个去重项目。审计将**同行评审论文、预印本、软件归档、项目文档和缺少可靠来源**分开；许可证能否再分发、功能是否经过科学验证也分别记录。DOI 不会自动被视为期刊论文。

## 1) 生成完整审计

- **作用：** 对 320 个项目逐条输出可复核的引用状态、许可证状态和待补材料。
- **处理方法：** 读取冻结的工具清单、已核验 DOI 记录和项目来源；不在线猜测 DOI。
- **API：** `dnakit.integrations.audit_publication_evidence()`。
- **输入：** 无输入时审计全部项目；可用 `query`、`manuscript_category` 或 `journal_citation_status` 筛选。
- **输出：** 每行包含 `journal_citation_status`、`manuscript_category`、`citation_evidence`、`license_status`、`required_metadata` 和 `scientific_validation`。
- **示例代码：**

```python
from dnakit.integrations import audit_publication_evidence

rows = audit_publication_evidence(verbose=False)
print(len(rows))
print(rows[0]["journal_citation_status"])
```

- **示例结果：**

```text
320
journal_article
```

## 2) 查看五类投稿处理结果

- **作用：** 把项目分成“可直接写入方法”“需要补版本/访问日期”“仅项目引用”“缺少可靠论文/DOI”“许可证受限”五类。
- **处理方法：** 许可证受限优先标记；引用状态和许可证状态不会互相覆盖。
- **API：** `dnakit.integrations.publication_audit_summary()`。
- **输入：** 无输入。
- **输出：** 返回五类计数、严格期刊状态计数、许可证计数和 `ready_for_journal`。
- **示例代码：**

```python
from dnakit.integrations import publication_audit_summary

summary = publication_audit_summary(verbose=False)
print(summary["manuscript_categories"])
print(summary["journal_citation_status"])
```

- **示例结果：**

```text
{'manuscript_ready': 84, 'needs_version_or_access_date': 4,
 'project_only': 47, 'missing_reliable_citation': 150,
 'license_restricted': 35}
{'journal_article': 86, 'preprint': 2, 'software_archive': 2,
 'project_only': 71, 'source_check_required': 5, 'missing': 154}
```

## 3) 验证单个项目是否有期刊论文

- **作用：** 检查一个工具能否使用同行评审论文作为方法引用。
- **处理方法：** `journal_article` 只表示已找到同行评审期刊论文；`preprint`、`software_archive` 和 `project_only` 不会被提升为期刊论文。
- **API：** `dnakit.integrations.describe_publication_audit(name[必须])`。
- **输入：** 清单中的项目名称，例如 `"BWA"`、`"pyfaidx"` 或 `"dnaio"`。
- **输出：** 返回 DOI、论文类别、建议的版本/访问日期字段、许可证状态和科学验证边界。
- **示例代码：**

```python
from dnakit.integrations import describe_publication_audit

for name in ("BWA", "pyfaidx", "dnaio"):
    row = describe_publication_audit(name, verbose=False)
    print(name, row["journal_citation_status"], row["manuscript_category"])
```

- **示例结果：**

```text
BWA journal_article manuscript_ready
pyfaidx preprint needs_version_or_access_date
dnaio software_archive needs_version_or_access_date
```

## 4) 论文中仍需记录的字段

- **作用：** 避免只写软件名称而无法复现实验。
- **处理方法：** 期刊论文项目至少补充精确软件版本；预印本和软件归档还要注明来源类型及访问日期；项目-only 项目必须保留官方项目 URL。
- **API：** 每条结果的 `required_metadata` 字段。
- **输入：** 审计结果中的一行。
- **输出：** 例如 `exact_software_version`、`access_date`、`label_non_journal_source`。
- **示例代码：**

```python
row = describe_publication_audit("Centrifuge", verbose=False)
print(row["required_metadata"])
```

- **示例结果：**

```text
['exact_software_version', 'access_date', 'label_non_journal_source']
```

## 5) 许可证和科学验证边界

- DOI 或论文只证明引用来源，不授予源码、模型、数据库或远程服务的再分发权。
- `license_status` 只反映当前清单中的许可证证据；`restricted_or_unverified` 必须人工复核许可条款。
- `scientific_validation` 固定为 `not_assessed_by_publication_audit`；论文引用审计不会替代基准数据、阳性/阴性对照、统计检验或独立复现。
