DNAKit Galaxy tools
===================

该目录包含四个 Galaxy Tool Shed 包装器：

* ``DNAKit convert``：FASTA/FASTQ 格式转换；
* ``DNAKit deduplicate``：FASTA 序列去重；
* ``DNAKit split``：可复现地划分训练、验证和测试集；
* ``DNAKit report``：生成自包含 HTML 报告。

包装器固定依赖 ``dnakit=0.1.1``。正式上传 Galaxy Tool Shed 前，必须先让
这个版本能够从 Bioconda 解析。

本地静态检查：

.. code-block:: bash

   planemo lint galaxy/dnakit
   planemo shed_lint galaxy/dnakit

本地运行 Galaxy 测试：

.. code-block:: bash

   planemo test galaxy/dnakit

上传需要 Galaxy Tool Shed 账号和 API Key。先将密钥保存到
``GALAXY_TOOL_SHED_API_KEY`` 环境变量，再创建 Test Tool Shed 仓库：

.. code-block:: bash

   planemo shed_create galaxy/dnakit \
     --shed_target testtoolshed \
     --shed_key_from_env GALAXY_TOOL_SHED_API_KEY

测试通过后，再将 ``--shed_target`` 改为 ``toolshed``，创建 Main Tool Shed 仓库。
