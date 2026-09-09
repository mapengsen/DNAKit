DNAKit Galaxy tools
===================

该目录包含四个 Galaxy Tool Shed 包装器：

* ``DNAKit convert``：FASTA/FASTQ 格式转换；
* ``DNAKit deduplicate``：FASTA 序列去重；
* ``DNAKit split``：可复现地划分训练、验证和测试集；
* ``DNAKit report``：生成自包含 HTML 报告。

包装器固定依赖 ``dnakit=0.1.1``，该版本已在 Bioconda 正式发布。

本地静态检查：

.. code-block:: bash

   planemo lint galaxy/dnakit
   planemo shed_lint galaxy/dnakit

本地运行 Galaxy 测试：

.. code-block:: bash

   planemo test galaxy/dnakit

上传需要目标 Tool Shed 账号和对应站点的 API Key。Test Tool Shed 与
Main Tool Shed 的密钥不可混用。将正式站密钥保存到
``GALAXY_TOOL_SHED_API_KEY`` 环境变量，测试通过后上传 Main Tool Shed：

.. code-block:: bash

   planemo shed_update galaxy/dnakit \
     --shed_target toolshed \
     --shed_key_from_env GALAXY_TOOL_SHED_API_KEY \
     --force_repository_creation

该命令会在仓库不存在时创建仓库，并上传工具和元数据。测试站上传时，
使用测试站密钥并将 ``--shed_target`` 改为 ``testtoolshed``。
