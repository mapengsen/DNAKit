# DNA Cauldron 计划与骨架结果

本页记录 P35 相关的真实 DNAKit 接口。函数按所需可选依赖延迟导入；坐标、默认值与参数选择以各接口说明为准。

验证范围说明：Backbone tables report supplied insertion results; they do not perform backbone selection.

## `dnakit.molecular_lab.dnacauldron_p35_functions.select_dnacauldron_connectors`

**作用：** Select the bridging connector parts needed to close a circular assembly (``autoselect_connectors``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def select_dnacauldron_connectors(parts: dict[str, str], connectors: dict[str, str], *, enzymes: list[str] | None=None, topology: str='auto', filter_site_fragments: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Select the bridging connector parts needed to close a circular assembly (``autoselect_connectors``).

``parts`` are the main parts of the assembly (one per slot) and ``connectors`` is a pool of candidate
bridging parts; both are circular DNA strings with the enzyme sites by default. When a cycle that
contains every main part exists, the connectors it uses are returned (with the circular constructs of
the completed mix); otherwise ``ConfigurationError`` says that no construct was found.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.simulate_dnacauldron_lcr_assembly`

**作用：** Simulate a ligase cycling reaction assembly of linear parts joined by bridging oligos.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def simulate_dnacauldron_lcr_assembly(parts: dict[str, str], bridging_oligos: dict[str, str], *, circular: bool=True, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, min_parts: int=2, max_results: int=20, annotate_homologies: bool=True, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate a ligase cycling reaction assembly of linear parts joined by bridging oligos.

``parts`` are linear DNA strings; each bridging oligo must span a junction (its two halves are
checked with DNA Cauldron's ``HomologyChecker`` defaults: 15 to 80 bases) and is written as it appears
in the concatenation of the two neighbouring parts. With ``circular=True`` the circular constructs are
returned, otherwise the linear ones (at least ``min_parts`` parts). Each construct lists its parts and
the bridging oligos it used.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.find_dnacauldron_bridging_oligo`

**作用：** Find the bridging oligo that joins two fragments in a ligase cycling reaction.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def find_dnacauldron_bridging_oligo(left_sequence: str, right_sequence: str, bridging_oligos: dict[str, str], *, min_size: int=15, max_size: int | None=80, min_tm: float=0, max_tm: float | None=None, max_distance: int=0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Find the bridging oligo that joins two fragments in a ligase cycling reaction.

The oligo must occur in ``left_sequence + right_sequence`` across the junction, with both halves
accepted by the homology checker. Returns ``found`` and, when found, the oligo name and the span
``[start, end)`` of the oligo in the concatenation (``find_oligo_homology``). The first matching oligo in
name order is returned, as in DNA Cauldron.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.plot_dnacauldron_cuts`

**作用：** Plot a construct map that marks where restriction enzymes cut, to a new PNG file.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def plot_dnacauldron_cuts(record: dict[str, Any], enzymes: list[str], output_path: str, *, linear: str='auto', figure_width: float=5.0, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Plot a construct map that marks where restriction enzymes cut, to a new PNG file.

``record`` is a record dictionary; ``linear`` is ``auto`` (from the record topology), ``true`` or
``false``. Each cut position of each enzyme (Biopython ``search``, one-based) is added as a feature and
drawn on a linear or circular map (``plot_cuts``). The result gives the PNG path and the cut positions.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.list_dnacauldron_record_overhangs`

**作用：** Read the overhang sequences recorded as annotations in a construct record.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def list_dnacauldron_record_overhangs(record: dict[str, Any] | None=None, *, genbank_path: str | None=None, with_locations: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Read the overhang sequences recorded as annotations in a construct record.

Give a ``record`` dictionary or a GenBank file (``genbank_path``). Four-base ``homology`` features
(labels written when a construct is assembled with homologies annotated) are returned in position order;
when there are none, the 4-base gaps between ``From <part>`` features are used
(``list_overhangs_from_record_annotations``). ``with_locations`` also returns each start position.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.get_dnacauldron_insert_overhangs`

**作用：** Get the overhang pair of the insert of a Type-2S part by digesting it with an enzyme.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def get_dnacauldron_insert_overhangs(sequence: str, *, enzyme: str='BsmBI', standardize: bool=True, topology: str='circular', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Get the overhang pair of the insert of a Type-2S part by digesting it with an enzyme.

The part is digested into sticky-end fragments; the one fragment without a site is the insert and its
``(left, right)`` overhangs are returned. With ``standardize`` the pair is the smaller of ``(o1, o2)`` and
``(revcomp(o2), revcomp(o1))`` so that a part and its reverse-complement give the same pair
(``get_overhangs_from_record``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.anneal_dnacauldron_oligo_pair`

**作用：** Simulate the annealing of an oligo pair into a sticky-end fragment (``OligoPairAnnealing``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def anneal_dnacauldron_oligo_pair(oligo: str, other_oligo: str, *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate the annealing of an oligo pair into a sticky-end fragment (``OligoPairAnnealing``).

Both oligos are written 5' to 3'. One oligo must contain the reverse complement of the other; the
annealed part is the double-stranded ``sequence`` and the unpaired bases on each side are the sticky
ends (strand 1). Both orders are tried, so up to two products are returned; none raises
``ConfigurationError``. This is the adapter-making step of BASIC assembly.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.simulate_dnacauldron_basic_assembly`

**作用：** Simulate a BASIC assembly: BsaI adapters on each part, then sticky-end circular assembly.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def simulate_dnacauldron_basic_assembly(parts: dict[str, str], adapter_oligos: dict[str, list[str]], assembly_parts: list[str], *, max_constructs: int=40, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate a BASIC assembly: BsaI adapters on each part, then sticky-end circular assembly.

``parts`` are the linear DNA parts (with BsaI sites); ``adapter_oligos`` maps each adapter name to its
two oligos, which are annealed first (``OligoPairAnnealing``); ``assembly_parts`` is the ordered list of
triplets ``adapter, part, adapter`` (a multiple of 3 names). Each triplet is ligated through BsaI
overhangs and the adapted parts are assembled into circular constructs through their long linker ends.
Up to ``max_constructs`` constructs are returned; assembly flaws are listed in ``errors``.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.read_dnacauldron_assembly_plan`

**作用：** Import an assembly plan from a CSV or Excel spreadsheet (``AssemblyPlan.from_spreadsheet``).

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def read_dnacauldron_assembly_plan(plan_path: str, *, sheet_name: str='all', header: bool=False, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Import an assembly plan from a CSV or Excel spreadsheet (``AssemblyPlan.from_spreadsheet``).

Each row is ``assembly name, class: <type>, part, part, ...`` with the class one of
``type2s_assembly``, ``gibson_assembly``, ``BASIC_assembly``, ``biobrick_assembly``, ``oligo_annealing``
or ``lcr_assembly``; other cells may be ``parameter: value``. The rows of a CSV must have the same
number of columns (pad with empty cells). The result lists each assembly with its class, parts and
the dependency level that DNA Cauldron computes from the shared names.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.write_dnacauldron_assembly_plan`

**作用：** Export an assembly plan as a table, and optionally as a CSV file.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def write_dnacauldron_assembly_plan(assemblies: list[dict[str, Any]], *, output_path: str | None=None, assembly_class: str='type2s_assembly', verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Export an assembly plan as a table, and optionally as a CSV file.

``assemblies`` lists ``{"name", "parts"}`` (parts are part names or the names of other assemblies).
The table has one row per assembly, sorted by dependency level then name, with the columns
``assembly`` and ``parts`` (``AssemblyPlan.to_dataframe``). With ``output_path`` (a new ``.csv``
file) the plan is also written in DNA Cauldron's ``construct,parts`` format (``to_spreadsheet``).
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_functions.tabulate_dnacauldron_plan_constructs`

**作用：** Simulate an assembly plan and tabulate the constructs of every assembly.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def tabulate_dnacauldron_plan_constructs(plan_path: str, sequences_path: str, *, connectors_path: str | None=None, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Simulate an assembly plan and tabulate the constructs of every assembly.

The plan is a CSV or spreadsheet (see ``read_dnacauldron_assembly_plan``); ``sequences_path`` is a folder
of part files (FASTA, GenBank or SnapGene, the file name is the part name); ``connectors_path`` is an
optional folder whose sub-folders are connector collections. One row is returned per construct with the
columns of DNA Cauldron's construct table: construct id, parts, number of parts, size, assembly name,
dependencies and level. ``outcomes`` counts valid, errored and cancelled assemblies.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## `dnakit.molecular_lab.dnacauldron_p35_backbone_functions.tabulate_dnacauldron_backbone_choices`

**作用：** Return BackboneChoice.list_to_infos_spreadsheet as JSON records.

**处理方法：** 由模块中的显式接口校验输入并调用对应后端。

**API：**

```text
def tabulate_dnacauldron_backbone_choices(choices: Sequence[Mapping[str, Any]], *, verbose: bool=True) -> dict[str, Any]:
```

**输入：** 签名列出必选参数、可选参数与原生默认值；参数含义和输出如下。

```text
Return BackboneChoice.list_to_infos_spreadsheet as JSON records.

Each imported choice needs ``record`` (id, sequence, optional topology/features),
and may have ``already_on_backbone`` (bool), ``backbone_record``, ``final_record``
in the same record format, and ``error`` (text). This reports caller-supplied
insertion results; it does not run backbone selection or infer successful cloning.
```

**输出：** 返回可 JSON 序列化的字典；文件输出返回路径。可选依赖缺失时抛出带安装提示的 BackendUnavailableError。

## 可执行示例

**示例代码：**

```python
from dnakit.molecular_lab.dnacauldron_p35_backbone_functions import tabulate_dnacauldron_backbone_choices
r = tabulate_dnacauldron_backbone_choices([{'record': {'id': 'insert', 'sequence': 'ACGT'}, 'error': 'No compatible backbone'}], verbose=False)
print(r['choices'][0])
```

**示例输出：**

```text
{'original_record': 'insert', 'already_on_backbone': 'no', 'detected_backbone': '', 'final_record_length': 4, 'error': 'No compatible backbone'}
```
