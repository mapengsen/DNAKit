# Biotite 外部程序操作（P30）

这些接口执行真实上游程序。ViennaRNA 2.7.2 / Biotite 2.0.0 的 11 个 RNA 接口已完成真实对比：40 例一致，4 例单碱基输入由上游自身失败。MMseqs2、SRA 和 tantan 尚未完成真实程序结果验证。RNA 坐标必须为有限数值，异常上游输出会明确报错。

Sources: [Biotite 2 application API](https://www.biotite-python.org/latest/apidoc/biotite.application.html), [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump).

## mmseqs_align

**作用：** Re-align the hits of a prefilter or search result with Smith-Waterman (`mmseqs align`).

**处理方法：** Re-align the hits of a prefilter or search result with Smith-Waterman (`mmseqs align`).

Databases are given as path prefixes. ``result_db`` holds the candidate hits (for example
the output of ``mmseqs_prefilter``); the new ``alignment_db`` holds scored alignments that
``mmseqs_convertalis`` can export. An option left at None is not passed.

Parameters:
- query_db, target_db, result_db: existing databases; alignment_db: new alignment database.
- e_value, min_seq_id, coverage, cov_mode, min_aln_len: hit filters as in ``mmseqs_search``.
- alignment_mode: 0 automatic, 1 score and end position, 2 also start and coverage, 3 also
  sequence identity.
- alignment_output_mode: 0 automatic, 1-3 as above, 4 ungapped only, 5 score only.
- add_backtrace, realign, gap_open, gap_extend, max_accept, max_rejected, zdrop: as in
  ``mmseqs_search``.
- wrapped_scoring: double the nucleotide query to allow diagonals that wrap around its ends
  (circular sequences).
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_align(query_db: 'str', target_db: 'str', result_db: 'str', alignment_db: 'str', *, e_value: 'float | None' = None, min_seq_id: 'float | None' = None, coverage: 'float | None' = None, cov_mode: 'int | None' = None, min_aln_len: 'int | None' = None, alignment_mode: 'int | None' = None, alignment_output_mode: 'int | None' = None, add_backtrace: 'bool | None' = None, realign: 'bool | None' = None, gap_open: 'int | None' = None, gap_extend: 'int | None' = None, max_accept: 'int | None' = None, max_rejected: 'int | None' = None, zdrop: 'int | None' = None, wrapped_scoring: 'bool | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_convertmsa

**作用：** Convert Stockholm multiple alignments into an MMseqs2 MSA database (`mmseqs convertmsa`).

**处理方法：** Convert Stockholm multiple alignments into an MMseqs2 MSA database (`mmseqs convertmsa`).

Parameters:
- msa_file: Stockholm file (``.sto``, optionally gzip-compressed) with one or many
  alignments.
- msa_db: prefix of the new MSA database.
- identifier_field: Stockholm comment used as the MSA identifier: 0 ID, 1 AC (default 0;
  when the comment is missing the name of the first sequence is used).
- mmseqs_path: path of the mmseqs executable; default is the one on PATH.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_convertmsa(msa_file: 'str', msa_db: 'str', *, identifier_field: 'int | None' = None, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_createindex

**作用：** Precompute the k-mer index of a sequence database (`mmseqs createindex`).

**处理方法：** Precompute the k-mer index of a sequence database (`mmseqs createindex`).

The index is written next to the database as ``<sequence_db>.idx*`` and speeds up later
searches against it. The call refuses a database that already has an index.

Parameters:
- sequence_db: prefix of an existing sequence database.
- kmer_length: -k, k-mer length (0 chooses automatically).
- sensitivity: -s, sensitivity the index is built for (default 7.5).
- mask: mask low-complexity regions with tantan (default on).
- search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
  nucleotide alignment.
- split: number of index chunks (0 chooses automatically).
- spaced_kmer_mode, spaced_kmer_pattern: spaced k-mers (default on) and a custom 0/1 pattern.
- max_seqs: maximum prefilter results per query that the index is sized for.
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_createindex(sequence_db: 'str', *, kmer_length: 'int | None' = None, sensitivity: 'float | None' = None, mask: 'bool | None' = None, search_type: 'int | None' = None, split: 'int | None' = None, spaced_kmer_mode: 'bool | None' = None, spaced_kmer_pattern: 'str | None' = None, max_seqs: 'int | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_createsubdb

**作用：** Extract a subset of the entries of a database into a new database (`mmseqs createsubdb`).

**处理方法：** Extract a subset of the entries of a database into a new database (`mmseqs createsubdb`).

Parameters:
- subset: text file with one database key (or FASTA identifier, see ``id_mode``) per line,
  or the prefix of another database whose keys are used.
- input_db: prefix of the existing database.
- output_db: prefix of the new database.
- subdb_mode: 0 copy the data, 1 soft-link the data and write a new index (default 0).
- id_mode: 0 select by database key, 1 select by FASTA identifier through the .lookup file
  (default 0).
- mmseqs_path: path of the mmseqs executable; default is the one on PATH.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_createsubdb(subset: 'str', input_db: 'str', output_db: 'str', *, subdb_mode: 'int | None' = None, id_mode: 'int | None' = None, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_databases

**作用：** Download and prepare a named reference database (`mmseqs databases`).

**处理方法：** Download and prepare a named reference database (`mmseqs databases`).

This contacts the database provider over the network and can download many gigabytes; it
is never run by the ordinary tests.

Parameters:
- name: one of UniRef100, UniRef90, UniRef50, UniProtKB, UniProtKB/TrEMBL,
  UniProtKB/Swiss-Prot, NR, NT, GTDB, PDB, PDB70, Pfam-A.full, Pfam-A.seed, Pfam-B, CDD,
  eggNOG, VOGDB, dbCAN2, SILVA, Resfinder or Kalamari (NT, SILVA, Resfinder and Kalamari
  hold nucleotide sequences).
- sequence_db: prefix of the new database.
- remove_tmp_files: delete the downloaded temporary files afterwards.
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_databases(name: 'str', sequence_db: 'str', *, remove_tmp_files: 'bool | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_msa2profile

**作用：** Turn an MSA database into a profile database (`mmseqs msa2profile`).

**处理方法：** Turn an MSA database into a profile database (`mmseqs msa2profile`).

Parameters:
- msa_db: prefix of an MSA database (from ``mmseqs_convertmsa`` or ``mmseqs_result2msa``).
- profile_db: prefix of the new profile database.
- msa_type: 0 cA3M, 1 A3M, 2 FASTA (default 2).
- match_mode: 0 keep columns with a residue in the first sequence, 1 keep columns with a
  residue in ``match_ratio`` of all sequences (default 0).
- match_ratio: column occupancy for match_mode 1 (default 0.5).
- skip_query: leave the first (query) sequence out of the profile.
- filter_msa: filter the MSA before building the profile (default on).
- max_seq_id: maximum pairwise identity kept by the filter (default 0.9).
- coverage: minimum fraction of the query covered by member sequences (default 0.0).
- gap_open, gap_extend: gap costs used for the profile.
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_msa2profile(msa_db: 'str', profile_db: 'str', *, msa_type: 'int | None' = None, match_mode: 'int | None' = None, match_ratio: 'float | None' = None, skip_query: 'bool | None' = None, filter_msa: 'bool | None' = None, max_seq_id: 'float | None' = None, coverage: 'float | None' = None, gap_open: 'int | None' = None, gap_extend: 'int | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_profile2consensus

**作用：** Write the consensus sequence of every profile into a sequence database

**处理方法：** Write the consensus sequence of every profile into a sequence database
(`mmseqs profile2consensus`).

Parameters:
- profile_db: prefix of a profile database.
- sequence_db: prefix of the new sequence database (one consensus per profile).
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_profile2consensus(profile_db: 'str', sequence_db: 'str', *, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_profile2pssm

**作用：** Export profiles as position-specific scoring matrices (`mmseqs profile2pssm`).

**处理方法：** Export profiles as position-specific scoring matrices (`mmseqs profile2pssm`).

Parameters:
- profile_db: prefix of a profile database.
- pssm_file: new text file with one PSSM per profile (a database prefix when
  ``db_output`` is true).
- db_output: return a result database instead of a text file (default off).
- threads, mmseqs_path: as in ``mmseqs_search``.

For the text output the result adds ``pssm_preview`` (the first 20 lines).

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_profile2pssm(profile_db: 'str', pssm_file: 'str', *, db_output: 'bool | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_profile2repseq

**作用：** Write the representative member sequence of every profile (`mmseqs profile2repseq`).

**处理方法：** Write the representative member sequence of every profile (`mmseqs profile2repseq`).

The representative is the profile's own query sequence as stored in the profile database.

Parameters:
- profile_db: prefix of a profile database.
- sequence_db: prefix of the new sequence database.
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_profile2repseq(profile_db: 'str', sequence_db: 'str', *, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_result2profile

**作用：** Build a profile database from search results (`mmseqs result2profile`).

**处理方法：** Build a profile database from search results (`mmseqs result2profile`).

For each query the hits of ``result_db`` (an alignment database from ``mmseqs_search``
run with backtraces, ``add_backtrace=True``) are stacked into an MSA and converted to a
profile.

Parameters:
- query_db, target_db, result_db: existing databases; profile_db: new profile database.
- e_profile: include hits with an E-value below this threshold (default 0.001).
- mask_profile: mask the query with tantan before building the profile (default on).
- filter_msa, max_seq_id, coverage: MSA filtering as in ``mmseqs_msa2profile``.
- gap_open, gap_extend: gap costs used for the profile.
- threads, mmseqs_path: as in ``mmseqs_search``.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_result2profile(query_db: 'str', target_db: 'str', result_db: 'str', profile_db: 'str', *, e_profile: 'float | None' = None, mask_profile: 'bool | None' = None, filter_msa: 'bool | None' = None, max_seq_id: 'float | None' = None, coverage: 'float | None' = None, gap_open: 'int | None' = None, gap_extend: 'int | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_search

**作用：** Search a query database against a target database (`mmseqs search`).

**处理方法：** Search a query database against a target database (`mmseqs search`).

Databases are given as path prefixes (create them with ``mmseqs_createdb``). The result is an
alignment database that ``mmseqs_convertalis`` turns into a table. An option left at None is
not passed, so the MMseqs2 default applies.

Parameters:
- query_db, target_db: prefixes of the query and target sequence (or profile) databases.
- alignment_db: prefix of the new alignment database.
- sensitivity: -s, 1.0 faster, 4.0 fast, 7.5 sensitive (default 5.7).
- kmer_length: -k, k-mer length (0 chooses automatically).
- e_value: -e, list matches below this E-value (default 0.001).
- min_seq_id: --min-seq-id, minimum sequence identity 0.0-1.0.
- coverage, cov_mode: -c and --cov-mode (0 query and target, 1 target, 2 query, 3-5 length
  criteria).
- min_aln_len: minimum alignment length.
- max_seqs: --max-seqs, results per query passing the prefilter (default 300).
- search_type: 0 automatic, 1 amino acid, 2 translated, 3 nucleotide, 4 translated
  nucleotide alignment.
- num_iterations: iterative profile search rounds (default 1); start_sens and sens_steps set
  the sensitivity ladder.
- mask: low-complexity masking with tantan in the prefilter (default on).
- strand: 0 reverse, 1 forward, 2 both (nucleotide searches only; default 1).
- add_backtrace: -a, store the alignment backtrace so that convertalis can write CIGARs.
- realign: shorter, more conservative alignments.
- gap_open, gap_extend, max_accept, max_rejected, zdrop: alignment scoring and stopping
  parameters (zdrop applies to nucleotide alignment).
- threads: CPU cores used (default 1; MMseqs2 itself defaults to all cores).
- mmseqs_path: path of the mmseqs executable; default is the one on PATH.

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_search(query_db: 'str', target_db: 'str', alignment_db: 'str', *, sensitivity: 'float | None' = None, kmer_length: 'int | None' = None, e_value: 'float | None' = None, min_seq_id: 'float | None' = None, coverage: 'float | None' = None, cov_mode: 'int | None' = None, min_aln_len: 'int | None' = None, max_seqs: 'int | None' = None, search_type: 'int | None' = None, num_iterations: 'int | None' = None, start_sens: 'float | None' = None, sens_steps: 'int | None' = None, mask: 'bool | None' = None, strand: 'int | None' = None, add_backtrace: 'bool | None' = None, realign: 'bool | None' = None, gap_open: 'int | None' = None, gap_extend: 'int | None' = None, max_accept: 'int | None' = None, max_rejected: 'int | None' = None, zdrop: 'int | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mmseqs_unpackdb

**作用：** Unpack every entry of a database into its own file (`mmseqs unpackdb`).

**处理方法：** Unpack every entry of a database into its own file (`mmseqs unpackdb`).

Typical use is an MSA database, which is unpacked into one alignment file per entry.

Parameters:
- database: prefix of an existing database.
- output_dir: new folder that receives the files.
- unpack_name_mode: 0 name files by database key, 1 by accession through the .lookup file
  (default 1).
- unpack_suffix: file suffix for the unpacked files, for example ``.a3m``; add ``.gz`` to
  write compressed files.
- threads, mmseqs_path: as in ``mmseqs_search``.

The result adds ``files`` (the unpacked file names).

**API：** `dnakit.alignment.mmseqs2_p30_cli_functions.mmseqs_unpackdb(database: 'str', output_dir: 'str', *, unpack_name_mode: 'int | None' = None, unpack_suffix: 'str | None' = None, threads: 'int' = 1, mmseqs_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## prefetch_sra_run

**作用：** Download an SRA run archive with NCBI prefetch; network access is required.

**处理方法：** Download an SRA run archive with NCBI prefetch; network access is required.

**API：** `dnakit.download.sra_p30_cli_functions.prefetch_sra_run(accession: 'str', output_dir: 'str | Path', *, executable_path: 'str | None' = None, timeout_seconds: 'float' = 3600, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## extract_sra_reads

**作用：** Extract FASTQ or FASTA with fasterq-dump; optionally use an existing local SRA archive.

**处理方法：** Extract FASTQ or FASTA with fasterq-dump; optionally use an existing local SRA archive.

``output_format`` is ``fastq`` (default) or ``fasta``. Without ``local_sra`` the
tool retrieves the run from NCBI. SRA Toolkit's default read splitting is used;
all single and paired output files are returned. No existing output is replaced.

**API：** `dnakit.download.sra_p30_cli_functions.extract_sra_reads(accession: 'str', output_dir: 'str | Path', *, output_format: 'str' = 'fastq', local_sra: 'str | Path | None' = None, threads: 'int | None' = None, executable_path: 'str | None' = None, timeout_seconds: 'float' = 3600, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

示例执行需要用户准备真实 MMseqs 数据库或测序 run；此处不提供未经执行的成功输出。

## mask_tandem_repeats

**作用：** Run tantan on DNA strings and return repeat masks and lowercase-masked sequences.

**处理方法：** Run tantan on DNA strings and return repeat masks and lowercase-masked sequences.

Bases are A/C/G/T/N (case insensitive); the original length and order are preserved.
True mask entries denote tandem repeat positions, following TantanApp.run.
The executable is external and must be installed by the caller.

**API：** `dnakit.sequence.tantan_p30_functions.mask_tandem_repeats(sequences: 'Sequence[str]', *, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.sequence.tantan_p30_functions import mask_tandem_repeats
from dnakit.exceptions import BackendUnavailableError
try:
    result = mask_tandem_repeats(["ACGTACGTACGT"], verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
external executable not installed
```

## fold_rna

**作用：** RNAfold minimum-free-energy structure of one nucleotide strand.

**处理方法：** RNAfold minimum-free-energy structure of one nucleotide strand.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna(sequence: 'str', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna("GGGAAACCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## fold_rna_pseudoknots

**作用：** RNAPKplex secondary structure, including pseudoknots.

**处理方法：** RNAPKplex secondary structure, including pseudoknots.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna_pseudoknots(sequence: 'str', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna_pseudoknots
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna_pseudoknots("GGGAAACCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## fold_rna_local

**作用：** RNALfold locally stable structures with positions in the original sequence.

**处理方法：** RNALfold locally stable structures with positions in the original sequence.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna_local(sequence: 'str', *, max_span: 'int' = 150, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna_local
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna_local("GGGAAACCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## cofold_rna

**作用：** RNAcofold minimum-free-energy structure of two strands, including intramolecular pairs.

**处理方法：** RNAcofold minimum-free-energy structure of two strands, including intramolecular pairs.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.cofold_rna(sequence1: 'str', sequence2: 'str', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import cofold_rna
from dnakit.exceptions import BackendUnavailableError
try:
    result = cofold_rna("GGGAAA", "UUUCCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## duplex_rna

**作用：** RNAduplex hybridization restricted to intermolecular base pairs.

**处理方法：** RNAduplex hybridization restricted to intermolecular base pairs.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.duplex_rna(sequence1: 'str', sequence2: 'str', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import duplex_rna
from dnakit.exceptions import BackendUnavailableError
try:
    result = duplex_rna("GGGAAA", "UUUCCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## find_rna_interactions

**作用：** RNAplex interaction-site search, returning all reported interactions and spans.

**处理方法：** RNAplex interaction-site search, returning all reported interactions and spans.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.find_rna_interactions(target: 'str', query: 'str', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import find_rna_interactions
from dnakit.exceptions import BackendUnavailableError
try:
    result = find_rna_interactions("GGGAAA", "UUUCCC", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## fold_rna_complex

**作用：** RNAmultifold minimum-free-energy structure of a multi-strand complex.

**处理方法：** RNAmultifold minimum-free-energy structure of a multi-strand complex.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna_complex(sequences: 'Sequence[str]', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna_complex
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna_complex(["GGGAAA", "UUUCCC"], verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## fold_rna_alignment

**作用：** RNAalifold consensus structure with physical and covariance energies.

**处理方法：** RNAalifold consensus structure with physical and covariance energies.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna_alignment(alignment: 'Sequence[str]', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna_alignment
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna_alignment(["GGGAAACCC", "GGGAAA-CC"], verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## fold_rna_alignment_local

**作用：** RNALalifold locally stable consensus structures and alignment-column spans.

**处理方法：** RNALalifold locally stable consensus structures and alignment-column spans.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.fold_rna_alignment_local(alignment: 'Sequence[str]', *, max_span: 'int' = 70, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import fold_rna_alignment_local
from dnakit.exceptions import BackendUnavailableError
try:
    result = fold_rna_alignment_local(["GGGAAACCC", "GGGAAA-CC"], verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## duplex_rna_alignments

**作用：** RNAaliduplex hybridization structure of two nucleotide alignments.

**处理方法：** RNAaliduplex hybridization structure of two nucleotide alignments.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.duplex_rna_alignments(alignment1: 'Sequence[str]', alignment2: 'Sequence[str]', *, temperature: 'float' = 37, executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import duplex_rna_alignments
from dnakit.exceptions import BackendUnavailableError
try:
    result = duplex_rna_alignments(["GGGAAA", "GGGAAA"], ["UUUCCC", "UUUCCC"], verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

## rna_structure_coordinates

**作用：** RNAplot 2D coordinates for each base of a dot-bracket secondary structure.

**处理方法：** RNAplot 2D coordinates for each base of a dot-bracket secondary structure.

**API：** `dnakit.secondary_structure.viennarna_p30_functions.rna_structure_coordinates(dot_bracket: 'str', *, layout: 'str' = 'naview', executable_path: 'str | None' = None, timeout_seconds: 'float' = 300, verbose: bool = True) -> 'dict[str, Any]'`

**输入：** 根据 API 提供序列、真实本地输入/数据库路径或 SRA run 编号；可选参数指定已安装的程序路径。

**输出：** 返回实际结构、能量、掩码、坐标或已核实存在的文件路径与来源信息；程序缺失时抛出 BackendUnavailableError。

**示例代码：**

```python
from dnakit.secondary_structure.viennarna_p30_functions import rna_structure_coordinates
from dnakit.exceptions import BackendUnavailableError
try:
    result = rna_structure_coordinates("(((...)))", verbose=False)
    print(result["status"])
except BackendUnavailableError:
    print("external executable not installed")
```

**示例输出：**

```text
completed
```

