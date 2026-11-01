# gtars refget 存储：序列、集合与别名

本页对应工作清单 P24 包中 gtars 的 refget 条目（`gtars.refget`）。refget 存储是一个本地目录，按 GA4GH sha512t24u 摘要保存序列与序列集合，并可带别名。仅在实体指定方式（摘要或别名）、单个或多个文件、元数据或碱基序列上有差别的条目，合并为带参数的同一个函数。测试版本为 gtars 0.10.0。下面的示例在同一个工作目录中按顺序运行。

## 示例文件

- **示例文件：** `a.fa`

```text
>chr1
ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
>chr2
TTGACCGGTAACCGGTTAGCATGCATGC
```

- **示例文件：** `b.fa`

```text
>chr1
ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
>chr3
GGGTTTAAACCCGGGTTTAAACCCGGG
```

## 1) `dnakit.annotation.refget_store_add_fasta()`

- **作用：** 对 FASTA 文件计算 GA4GH 摘要，把每个文件作为一个序列集合（连同序列）加入本地 refget 存储（不存在时自动创建）。
- **处理方法：** 单个文件调用 `RefgetStore.add_sequence_collection_from_fasta`，多个文件调用 `add_sequence_collections_from_fastas`；`use_collection_object=True` 时先用 `digest_fasta` 得到集合对象再用 `add_sequence_collection` 存入（结果相同）。已存在的集合默认跳过，`force=True` 时覆盖。
- **API：** `dnakit.annotation.refget_store_add_fasta(store, fasta_files, *, force=False, use_collection_object=False, verbose=True)`。
- **输入：** 存储目录、FASTA 路径列表、`force`、`use_collection_object`。
- **输出：** `collections`（文件、摘要、是否新增）、`collection_count`、`sequence_count`、`provenance`。
- **示例代码：**

```python
import shutil
from dnakit.annotation.gtars_p24_refget_functions import refget_store_add_fasta

shutil.rmtree('store', ignore_errors=True)
r = refget_store_add_fasta('store', ['a.fa', 'b.fa'], verbose=False)
for c in r['collections']:
    print(c['digest'], c['added'])
print(r['collection_count'], r['sequence_count'])
```

- **示例输出：**

```text
Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i True
KRN51Y8frKwPAQ7KRcxt2cnE3YAhTlHo True
2 3
```

## 2) `dnakit.annotation.refget_store_add_sequence()`

- **作用：** 对单条 DNA 序列计算摘要并存入存储（不属于任何集合）。
- **处理方法：** `gtars.refget.digest_sequence` 把序列转为大写并计算 sha512t24u、MD5 摘要和字母表；`RefgetStore.add_sequence` 存入记录并写出索引。
- **API：** `dnakit.annotation.refget_store_add_sequence(store, sequence, *, name=None, description=None, force=False, verbose=True)`。
- **输入：** 存储目录、序列、可选的 `name`、`description`、`force`。
- **输出：** `name`、`length`、`alphabet`、`sha512t24u`、`md5`、`provenance`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_store_add_sequence

r = refget_store_add_sequence('store', 'acgtnnacgt', name='x', verbose=False)
print(r['sha512t24u'], r['md5'], r['length'], r['alphabet'])
```

- **示例输出：**

```text
Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN 1c4ab0b3062cfd287a4c8b1d77ac5bd8 10 dna3bit
```

## 3) `dnakit.annotation.refget_collections_list()`

- **作用：** 列出 refget 存储中的序列集合，可按某个属性摘要过滤并分页。
- **处理方法：** 调用 `RefgetStore.list_collections(page, page_size, filters)`；过滤条件为一个 seqcol 属性（names、lengths、sequences 等）及其摘要。
- **API：** `dnakit.annotation.refget_collections_list(store, *, attribute=None, attribute_digest=None, page=0, page_size=100, verbose=True)`。
- **输入：** 存储目录、`attribute` 与 `attribute_digest`（需同时给出）、`page`、`page_size`。
- **输出：** `collections`（摘要、序列数及各属性摘要）、`pagination`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list

r = refget_collections_list('store', verbose=False)
for c in r['collections']:
    print(c['digest'], c['n_sequences'])
print(r['pagination'])
```

- **示例输出：**

```text
KRN51Y8frKwPAQ7KRcxt2cnE3YAhTlHo 2
Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i 2
{'page': 0, 'page_size': 100, 'total': 2}
```

## 4) `dnakit.annotation.refget_collection_digests()`

- **作用：** 以 seqcol 1 级（属性摘要）或 2 级（完整的名称、长度、序列摘要数组）返回一个集合。
- **处理方法：** `level=1` 调用 `RefgetStore.get_collection_level1`，`level=2` 调用 `get_collection_level2`。
- **API：** `dnakit.annotation.refget_collection_digests(store, digest, *, level=1, verbose=True)`。
- **输入：** 存储目录、集合摘要、`level`（1 或 2）。
- **输出：** `representation`（字典）、`level`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_digests, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
print(refget_collection_digests('store', d, verbose=False)['representation']['names'])
print(refget_collection_digests('store', d, level=2, verbose=False)['representation'])
```

- **示例输出：**

```text
oK9djEQruFwX_TMmEbS2qhoPpixthiCj
{'names': ['chr1', 'chr3'], 'lengths': [33, 27], 'sequences': ['SQ.z-pMq0F5z9iTxrIjy78TXVEhhix_H5rc', 'SQ.ou6fkS6aEnM-xOCptFiIOqDyl8yAhTSy']}
```

## 5) `dnakit.annotation.refget_collection_retrieve()`

- **作用：** 按摘要或（命名空间，别名）取回一个集合：集合元数据、每条序列的元数据，并可选返回碱基序列。
- **处理方法：** 用 `get_collection_metadata` 与 `get_collection_level2` 取元数据；`include_sequences=True` 时逐条加载（`load_sequence`、`get_sequence`）并解码。`max_bases` 限制输出规模。
- **API：** `dnakit.annotation.refget_collection_retrieve(store, *, digest=None, namespace=None, alias=None, include_sequences=False, max_bases=10000000, verbose=True)`。
- **输入：** 存储目录、`digest` 或 `namespace` + `alias`、`include_sequences`、`max_bases`。
- **输出：** 集合各摘要、`total_bases`、`sequences`（名称、长度、字母表、摘要，可选 `sequence`）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_retrieve, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_collection_retrieve('store', digest=d, include_sequences=True, verbose=False)
print(r['total_bases'])
for s in r['sequences']:
    print(s['name'], s['length'], s['sequence'])
```

- **示例输出：**

```text
60
chr1 33 ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
chr3 27 GGGTTTAAACCCGGGTTTAAACCCGGG
```

## 6) `dnakit.annotation.refget_sequence_retrieve()`

- **作用：** 按摘要（sha512t24u 或 MD5）或（命名空间，别名）取回一条序列：元数据，并可选返回碱基序列。
- **处理方法：** 调用 `get_sequence_metadata`（MD5 时用 `get_sequence`）；`include_sequence=True` 时加载并解码。
- **API：** `dnakit.annotation.refget_sequence_retrieve(store, *, digest=None, namespace=None, alias=None, include_sequence=False, max_bases=10000000, verbose=True)`。
- **输入：** 存储目录、`digest` 或 `namespace` + `alias`、`include_sequence`、`max_bases`。
- **输出：** `name`、`length`、`alphabet`、`sha512t24u`、`md5`，可选 `sequence`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_sequence_retrieve

r = refget_sequence_retrieve('store', digest='Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', include_sequence=True, verbose=False)
print(r['name'], r['length'], r['sequence'])
r2 = refget_sequence_retrieve('store', digest=r['md5'], verbose=False)
print(r2['sha512t24u'] == r['sha512t24u'])
```

- **示例输出：**

```text
x 10 ACGTNNACGT
True
```

## 7) `dnakit.annotation.refget_sequence_substring()`

- **作用：** 按 0 起点、半开区间提取一条存储序列的子串。
- **处理方法：** `RefgetStore.get_substrings` 只读取所需字节。序列可由摘要、（命名空间，别名）或“集合摘要 + 序列名”指定。
- **API：** `dnakit.annotation.refget_sequence_substring(store, ranges, *, digest=None, namespace=None, alias=None, collection=None, name=None, verbose=True)`。
- **输入：** 存储目录、`ranges`（[start, end] 对）及一种序列指定方式。
- **输出：** `substrings`（start、end、sequence）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list, refget_sequence_substring

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_sequence_substring('store', [[0, 10], [5, 8], [30, 33]], collection=d, name='chr1', verbose=False)
for s in r['substrings']:
    print(s['start'], s['end'], s['sequence'])
```

- **示例输出：**

```text
0 10 ACGTACGTAG
5 8 CGT
30 33 AGC
```

## 8) `dnakit.annotation.refget_collection_compare()`

- **作用：** 用 seqcol 比较规范比较存储中的两个集合：共有/独有的属性、数组元素重合数及顺序是否一致。
- **处理方法：** 调用 `RefgetStore.compare(digest_a, digest_b)`。
- **API：** `dnakit.annotation.refget_collection_compare(store, digest_a, digest_b, *, verbose=True)`。
- **输入：** 存储目录与两个集合摘要。
- **输出：** `comparison`，含 `digests`、`attributes`、`array_elements`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_compare, refget_collections_list

cs = refget_collections_list('store', verbose=False)['collections']
r = refget_collection_compare('store', cs[0]['digest'], cs[1]['digest'], verbose=False)['comparison']
print(r['attributes'])
print(sorted(r['array_elements']['a_and_b_count'].items()))
```

- **示例输出：**

```text
{'a_only': [], 'b_only': [], 'a_and_b': ['lengths', 'name_length_pairs', 'names', 'sequences', 'sorted_name_length_pairs', 'sorted_sequences']}
[('lengths', 1), ('name_length_pairs', 1), ('names', 1), ('sequences', 1), ('sorted_name_length_pairs', 1), ('sorted_sequences', 1)]
```

## 9) `dnakit.annotation.refget_collections_by_attribute()`

- **作用：** 查找带有给定属性摘要（例如相同的名称列表）的集合。
- **处理方法：** 调用 `RefgetStore.find_collections_by_attribute`；属性可为 names、lengths、sequences、name_length_pairs、sorted_name_length_pairs、sorted_sequences。
- **API：** `dnakit.annotation.refget_collections_by_attribute(store, attribute, attribute_digest, *, verbose=True)`。
- **输入：** 存储目录、属性名、属性摘要。
- **输出：** `collections`（排序后的摘要）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_digests, refget_collections_by_attribute, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
names = refget_collection_digests('store', d, verbose=False)['representation']['names']
print(refget_collections_by_attribute('store', 'names', names, verbose=False)['collections'] == [d])
```

- **示例输出：**

```text
True
```

## 10) `dnakit.annotation.refget_collection_attribute_array()`

- **作用：** 返回属性摘要背后的数组（名称、长度或序列摘要）。
- **处理方法：** 调用 `RefgetStore.get_attribute(attr_name, attr_digest)`；摘要未知时 `values` 为 null。
- **API：** `dnakit.annotation.refget_collection_attribute_array(store, attribute, attribute_digest, *, verbose=True)`。
- **输入：** 存储目录、属性名、属性摘要。
- **输出：** `values`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_attribute_array, refget_collection_digests, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
lengths = refget_collection_digests('store', d, verbose=False)['representation']['lengths']
print(refget_collection_attribute_array('store', 'lengths', lengths, verbose=False)['values'])
```

- **示例输出：**

```text
[33, 27]
```

## 11) `dnakit.annotation.refget_alias_register()`

- **作用：** 为序列或集合登记别名：单个别名，或两列 TSV 文件批量登记。
- **处理方法：** 单个别名用 `add_sequence_alias` / `add_collection_alias`，文件用 `load_sequence_aliases` / `load_collection_aliases`（别名、不带 `SQ.` 前缀的摘要）。别名保存在存储目录中。
- **API：** `dnakit.annotation.refget_alias_register(store, kind, namespace, *, alias=None, digest=None, alias_file=None, verbose=True)`。
- **输入：** 存储目录、`kind`（sequence 或 collection）、`namespace`，以及 `alias` + `digest` 或 `alias_file`。
- **输出：** `aliases_registered`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_register, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
print(refget_alias_register('store', 'collection', 'demo', alias='v1', digest=d, verbose=False)['aliases_registered'])
print(refget_alias_register('store', 'sequence', 'demo', alias='x', digest='Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', verbose=False)['aliases_registered'])
open('aliases.tsv', 'w').write('second\tYb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN\n')
print(refget_alias_register('store', 'sequence', 'fromfile', alias_file='aliases.tsv', verbose=False)['aliases_registered'])
```

- **示例输出：**

```text
1
1
1
```

## 12) `dnakit.annotation.refget_alias_list()`

- **作用：** 列出存储的别名命名空间，或某个命名空间下的别名。
- **处理方法：** 未给 `namespace` 时调用 `list_sequence_alias_namespaces` / `list_collection_alias_namespaces`，给出时调用 `list_sequence_aliases` / `list_collection_aliases`。
- **API：** `dnakit.annotation.refget_alias_list(store, kind, *, namespace=None, verbose=True)`。
- **输入：** 存储目录、`kind`、可选 `namespace`。
- **输出：** `namespaces` 或 `aliases`（命名空间不存在时为 null）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_list

print(refget_alias_list('store', 'sequence', verbose=False)['namespaces'])
print(refget_alias_list('store', 'collection', namespace='demo', verbose=False)['aliases'])
print(refget_alias_list('store', 'collection', namespace='none', verbose=False)['aliases'])
```

- **示例输出：**

```text
['demo', 'fromfile']
['v1']
None
```

## 13) `dnakit.annotation.refget_alias_reverse_lookup()`

- **作用：** 返回指向某个序列或集合摘要的所有（命名空间，别名）对。
- **处理方法：** 调用 `get_aliases_for_sequence` / `get_aliases_for_collection`。
- **API：** `dnakit.annotation.refget_alias_reverse_lookup(store, kind, digest, *, verbose=True)`。
- **输入：** 存储目录、`kind`、摘要。
- **输出：** `aliases`（namespace、alias）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_reverse_lookup

r = refget_alias_reverse_lookup('store', 'sequence', 'Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', verbose=False)
print(r['aliases'])
```

- **示例输出：**

```text
[{'namespace': 'demo', 'alias': 'x'}, {'namespace': 'fromfile', 'alias': 'second'}]
```

## 14) `dnakit.annotation.refget_alias_remove()`

- **作用：** 删除一个序列或集合别名。
- **处理方法：** 调用 `remove_sequence_alias` / `remove_collection_alias`；别名不存在时 `removed` 为 false。
- **API：** `dnakit.annotation.refget_alias_remove(store, kind, namespace, alias, *, verbose=True)`。
- **输入：** 存储目录、`kind`、`namespace`、`alias`。
- **输出：** `removed`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_remove

print(refget_alias_remove('store', 'sequence', 'demo', 'x', verbose=False)['removed'])
print(refget_alias_remove('store', 'sequence', 'demo', 'x', verbose=False)['removed'])
```

- **示例输出：**

```text
True
False
```

## 15) `dnakit.annotation.refget_alias_pull()`

- **作用：** 把远程 refget 存储的别名表拉取到本地缓存存储（需联网）。
- **处理方法：** 先 `RefgetStore.open_remote(cache, url)`，再 `pull_aliases(namespace, strategy)`；`strategy` 为 keep-ours、keep-theirs 或 notify。常规测试不会调用此函数。
- **API：** `dnakit.annotation.refget_alias_pull(remote_url, cache_dir, *, namespace=None, strategy='keep-ours', verbose=True)`。
- **输入：** `remote_url`、`cache_dir`、可选 `namespace`、`strategy`。
- **输出：** `report`（pulled、skipped、not_found、conflicts）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_pull

try:
    refget_alias_pull('ftp://example.org/store', 'cache', verbose=False)
except Exception as exc:
    print(type(exc).__name__)
```

- **示例输出：**

```text
ConfigurationError
```

## 16) `dnakit.annotation.refget_collection_import()`

- **作用：** 把一个集合连同序列、别名和 FHR 元数据从一个本地存储复制到另一个。
- **处理方法：** 调用 `RefgetStore.import_collection(source, digest)`；两个存储都必须是磁盘存储。
- **API：** `dnakit.annotation.refget_collection_import(store, source_store, digest, *, verbose=True)`。
- **输入：** 目标存储、源存储、集合摘要。
- **输出：** `imported`（集合元数据）。
- **示例代码：**

```python
import shutil
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_import, refget_collections_list

shutil.rmtree('copy', ignore_errors=True)
d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_collection_import('copy', 'store', d, verbose=False)
print(r['imported']['digest'] == d, r['imported']['n_sequences'])
```

- **示例输出：**

```text
True 2
```

## 17) `dnakit.annotation.refget_orphan_removal_plan()`

- **作用：** 预览移除某个集合并清理孤立序列时会删除哪些序列；不做任何修改。
- **处理方法：** 调用 `RefgetStore.plan_orphan_removal(digest)`（仅供参考，不加锁）。
- **API：** `dnakit.annotation.refget_orphan_removal_plan(store, digest, *, verbose=True)`。
- **输入：** 存储目录、集合摘要。
- **输出：** `orphan_sequences`（排序后的摘要）。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list, refget_orphan_removal_plan

d = refget_collections_list('copy', verbose=False)['collections'][0]['digest']
print(refget_orphan_removal_plan('copy', d, verbose=False)['orphan_sequences'])
```

- **示例输出：**

```text
['ou6fkS6aEnM-xOCptFiIOqDyl8yAhTSy', 'z-pMq0F5z9iTxrIjy78TXVEhhix_H5rc']
```

## 18) `dnakit.annotation.refget_collection_remove()`

- **作用：** 从存储中移除一个集合及其别名和元数据；可选同时移除不再被引用的序列。
- **处理方法：** 调用 `RefgetStore.remove_collection(digest, remove_orphan_sequences)`。
- **API：** `dnakit.annotation.refget_collection_remove(store, digest, *, remove_orphan_sequences=False, verbose=True)`。
- **输入：** 存储目录、集合摘要、`remove_orphan_sequences`。
- **输出：** `removed`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_remove, refget_collections_list

d = refget_collections_list('copy', verbose=False)['collections'][0]['digest']
print(refget_collection_remove('copy', d, remove_orphan_sequences=True, verbose=False)['removed'])
print(refget_collection_remove('copy', d, verbose=False)['removed'])
print(refget_collections_list('copy', verbose=False)['collections'])
```

- **示例输出：**

```text
True
False
[]
```

## 19) `dnakit.annotation.refget_write_rgsi()`

- **作用：** 对 FASTA 文件计算摘要并写出 RGSI 索引：集合级摘要，随后每条序列一行。
- **处理方法：** 先 `digest_fasta`，再 `SequenceCollection.write_rgsi`。输出文件必须不存在。
- **API：** `dnakit.annotation.refget_write_rgsi(fasta, output, *, verbose=True)`。
- **输入：** FASTA 路径、输出路径。
- **输出：** `collection_digest`、`sequence_count`、`output`。
- **示例代码：**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_write_rgsi

import os
if os.path.exists('a.rgsi'):
    os.remove('a.rgsi')
r = refget_write_rgsi('a.fa', 'a.rgsi', verbose=False)
print(r['collection_digest'], r['sequence_count'])
print(open('a.rgsi').read().replace('\t', '|'))
```

- **示例输出：**

```text
Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i 2
##seqcol_digest=Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i
##names_digest=XEsH8IMZ09CBX17iXEWRagH50VGfARLo
##sequences_digest=JpKNkyzDHXHEV2j_6pLU7wn-z3SI_wn_
##lengths_digest=-TNjlmMFA1RgeHcONGSeNtHVdw4dGQRY
#name|length|alphabet|sha512t24u|md5|description
chr1|33|dna2bit|z-pMq0F5z9iTxrIjy78TXVEhhix_H5rc|9d0318e68b5b6254e514baf860d1ffba|
chr2|28|dna2bit|X0k13ENSst_-wZuv0dPtSr2C1_pWuUoP|561a5b2c52fe60ff79fb3cd65cef48d2|

```
