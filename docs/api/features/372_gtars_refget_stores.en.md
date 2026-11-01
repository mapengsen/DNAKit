# gtars refget stores: sequences, collections and aliases

This page covers the refget rows of gtars in work-list pack P24 (`gtars.refget`). A refget store is a local directory that holds sequences and sequence collections under their GA4GH sha512t24u digests, with optional aliases. Rows that differ only by the way an entity is named (digest or alias), by one or many files, or by metadata versus residues are merged into one function with parameters. gtars 0.10.0 was tested. Every example below runs in one working directory, in order.

## Example files

- **Example file:** `a.fa`

```text
>chr1
ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
>chr2
TTGACCGGTAACCGGTTAGCATGCATGC
```

- **Example file:** `b.fa`

```text
>chr1
ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
>chr3
GGGTTTAAACCCGGGTTTAAACCCGGG
```

## 1) `dnakit.annotation.refget_store_add_fasta()`

- **Purpose:** Digest FASTA files and add every file as a sequence collection, with its sequences, to a local refget store (created when missing).
- **Method:** One file calls `RefgetStore.add_sequence_collection_from_fasta`, several files call `add_sequence_collections_from_fastas`; `use_collection_object=True` digests each file with `digest_fasta` first and stores it with `add_sequence_collection` (same result). Existing collections are skipped unless `force=True`.
- **API:** `dnakit.annotation.refget_store_add_fasta(store, fasta_files, *, force=False, use_collection_object=False, verbose=True)`.
- **Input:** Store directory, a list of FASTA paths, `force`, `use_collection_object`.
- **Output:** `collections` (file, digest, added), `collection_count`, `sequence_count`, `provenance`.
- **Example code:**

```python
import shutil
from dnakit.annotation.gtars_p24_refget_functions import refget_store_add_fasta

shutil.rmtree('store', ignore_errors=True)
r = refget_store_add_fasta('store', ['a.fa', 'b.fa'], verbose=False)
for c in r['collections']:
    print(c['digest'], c['added'])
print(r['collection_count'], r['sequence_count'])
```

- **Example output:**

```text
Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i True
KRN51Y8frKwPAQ7KRcxt2cnE3YAhTlHo True
2 3
```

## 2) `dnakit.annotation.refget_store_add_sequence()`

- **Purpose:** Digest one DNA sequence and store it without any collection.
- **Method:** `gtars.refget.digest_sequence` uppercases the sequence and computes the sha512t24u and MD5 digests and the alphabet; `RefgetStore.add_sequence` stores the record and the index is written.
- **API:** `dnakit.annotation.refget_store_add_sequence(store, sequence, *, name=None, description=None, force=False, verbose=True)`.
- **Input:** Store directory, sequence, optional `name`, `description`, `force`.
- **Output:** `name`, `length`, `alphabet`, `sha512t24u`, `md5`, `provenance`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_store_add_sequence

r = refget_store_add_sequence('store', 'acgtnnacgt', name='x', verbose=False)
print(r['sha512t24u'], r['md5'], r['length'], r['alphabet'])
```

- **Example output:**

```text
Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN 1c4ab0b3062cfd287a4c8b1d77ac5bd8 10 dna3bit
```

## 3) `dnakit.annotation.refget_collections_list()`

- **Purpose:** List the sequence collections of a refget store, optionally only those with one attribute digest, with paging.
- **Method:** `RefgetStore.list_collections(page, page_size, filters)`; the filter is one seqcol attribute (names, lengths, sequences, ...) and its digest.
- **API:** `dnakit.annotation.refget_collections_list(store, *, attribute=None, attribute_digest=None, page=0, page_size=100, verbose=True)`.
- **Input:** Store directory, `attribute` and `attribute_digest` (together), `page`, `page_size`.
- **Output:** `collections` (digest, n_sequences and the attribute digests), `pagination`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list

r = refget_collections_list('store', verbose=False)
for c in r['collections']:
    print(c['digest'], c['n_sequences'])
print(r['pagination'])
```

- **Example output:**

```text
KRN51Y8frKwPAQ7KRcxt2cnE3YAhTlHo 2
Rzfz05pYxPktUFWK7lrdLByo_yCj9W_i 2
{'page': 0, 'page_size': 100, 'total': 2}
```

## 4) `dnakit.annotation.refget_collection_digests()`

- **Purpose:** Return a collection as seqcol level 1 (attribute digests) or level 2 (the full names, lengths and sequence digests).
- **Method:** `RefgetStore.get_collection_level1` (`level=1`) or `get_collection_level2` (`level=2`).
- **API:** `dnakit.annotation.refget_collection_digests(store, digest, *, level=1, verbose=True)`.
- **Input:** Store directory, collection digest, `level` (1 or 2).
- **Output:** `representation` (dict), `level`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_digests, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
print(refget_collection_digests('store', d, verbose=False)['representation']['names'])
print(refget_collection_digests('store', d, level=2, verbose=False)['representation'])
```

- **Example output:**

```text
oK9djEQruFwX_TMmEbS2qhoPpixthiCj
{'names': ['chr1', 'chr3'], 'lengths': [33, 27], 'sequences': ['SQ.z-pMq0F5z9iTxrIjy78TXVEhhix_H5rc', 'SQ.ou6fkS6aEnM-xOCptFiIOqDyl8yAhTSy']}
```

## 5) `dnakit.annotation.refget_collection_retrieve()`

- **Purpose:** Retrieve one collection by digest or by (namespace, alias): its metadata and the metadata of every sequence, and optionally the residues.
- **Method:** `get_collection_metadata` and `get_collection_level2` give the metadata; with `include_sequences=True` each sequence is loaded (`load_sequence`, `get_sequence`) and decoded. `max_bases` bounds the output.
- **API:** `dnakit.annotation.refget_collection_retrieve(store, *, digest=None, namespace=None, alias=None, include_sequences=False, max_bases=10000000, verbose=True)`.
- **Input:** Store directory, `digest` or `namespace` + `alias`, `include_sequences`, `max_bases`.
- **Output:** Collection digests, `total_bases`, `sequences` (name, length, alphabet, digests, optional `sequence`).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_retrieve, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_collection_retrieve('store', digest=d, include_sequences=True, verbose=False)
print(r['total_bases'])
for s in r['sequences']:
    print(s['name'], s['length'], s['sequence'])
```

- **Example output:**

```text
60
chr1 33 ACGTACGTAGCTAGCTAGGATCGATCGATTAGC
chr3 27 GGGTTTAAACCCGGGTTTAAACCCGGG
```

## 6) `dnakit.annotation.refget_sequence_retrieve()`

- **Purpose:** Retrieve one sequence by digest (sha512t24u or MD5) or by (namespace, alias): its metadata and optionally the residues.
- **Method:** `get_sequence_metadata` (or `get_sequence` for an MD5 key); with `include_sequence=True` the sequence is loaded and decoded.
- **API:** `dnakit.annotation.refget_sequence_retrieve(store, *, digest=None, namespace=None, alias=None, include_sequence=False, max_bases=10000000, verbose=True)`.
- **Input:** Store directory, `digest` or `namespace` + `alias`, `include_sequence`, `max_bases`.
- **Output:** `name`, `length`, `alphabet`, `sha512t24u`, `md5`, optional `sequence`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_sequence_retrieve

r = refget_sequence_retrieve('store', digest='Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', include_sequence=True, verbose=False)
print(r['name'], r['length'], r['sequence'])
r2 = refget_sequence_retrieve('store', digest=r['md5'], verbose=False)
print(r2['sha512t24u'] == r['sha512t24u'])
```

- **Example output:**

```text
x 10 ACGTNNACGT
True
```

## 7) `dnakit.annotation.refget_sequence_substring()`

- **Purpose:** Extract substrings of one stored sequence for 0-based half-open ranges.
- **Method:** `RefgetStore.get_substrings` reads only the needed bytes. The sequence is given by digest, by (namespace, alias) or by collection digest plus sequence name.
- **API:** `dnakit.annotation.refget_sequence_substring(store, ranges, *, digest=None, namespace=None, alias=None, collection=None, name=None, verbose=True)`.
- **Input:** Store directory, `ranges` ([start, end] pairs) and one way to name the sequence.
- **Output:** `substrings` (start, end, sequence).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list, refget_sequence_substring

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_sequence_substring('store', [[0, 10], [5, 8], [30, 33]], collection=d, name='chr1', verbose=False)
for s in r['substrings']:
    print(s['start'], s['end'], s['sequence'])
```

- **Example output:**

```text
0 10 ACGTACGTAG
5 8 CGT
30 33 AGC
```

## 8) `dnakit.annotation.refget_collection_compare()`

- **Purpose:** Compare two collections of a store with the seqcol comparison: which attributes they share, and which array elements and in which order.
- **Method:** `RefgetStore.compare(digest_a, digest_b)`.
- **API:** `dnakit.annotation.refget_collection_compare(store, digest_a, digest_b, *, verbose=True)`.
- **Input:** Store directory and two collection digests.
- **Output:** `comparison` with `digests`, `attributes`, `array_elements`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_compare, refget_collections_list

cs = refget_collections_list('store', verbose=False)['collections']
r = refget_collection_compare('store', cs[0]['digest'], cs[1]['digest'], verbose=False)['comparison']
print(r['attributes'])
print(sorted(r['array_elements']['a_and_b_count'].items()))
```

- **Example output:**

```text
{'a_only': [], 'b_only': [], 'a_and_b': ['lengths', 'name_length_pairs', 'names', 'sequences', 'sorted_name_length_pairs', 'sorted_sequences']}
[('lengths', 1), ('name_length_pairs', 1), ('names', 1), ('sequences', 1), ('sorted_name_length_pairs', 1), ('sorted_sequences', 1)]
```

## 9) `dnakit.annotation.refget_collections_by_attribute()`

- **Purpose:** Find the collections of a store that carry a given attribute digest, for example the same list of names.
- **Method:** `RefgetStore.find_collections_by_attribute`; attributes are names, lengths, sequences, name_length_pairs, sorted_name_length_pairs, sorted_sequences.
- **API:** `dnakit.annotation.refget_collections_by_attribute(store, attribute, attribute_digest, *, verbose=True)`.
- **Input:** Store directory, attribute name, attribute digest.
- **Output:** `collections` (sorted digests).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_digests, refget_collections_by_attribute, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
names = refget_collection_digests('store', d, verbose=False)['representation']['names']
print(refget_collections_by_attribute('store', 'names', names, verbose=False)['collections'] == [d])
```

- **Example output:**

```text
True
```

## 10) `dnakit.annotation.refget_collection_attribute_array()`

- **Purpose:** Return the array (names, lengths or sequence digests) that stands behind an attribute digest.
- **Method:** `RefgetStore.get_attribute(attr_name, attr_digest)`; `values` is null when the digest is unknown.
- **API:** `dnakit.annotation.refget_collection_attribute_array(store, attribute, attribute_digest, *, verbose=True)`.
- **Input:** Store directory, attribute name, attribute digest.
- **Output:** `values`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_attribute_array, refget_collection_digests, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
lengths = refget_collection_digests('store', d, verbose=False)['representation']['lengths']
print(refget_collection_attribute_array('store', 'lengths', lengths, verbose=False)['values'])
```

- **Example output:**

```text
[33, 27]
```

## 11) `dnakit.annotation.refget_alias_register()`

- **Purpose:** Register aliases for sequences or collections: one alias, or a two-column TSV file of aliases.
- **Method:** `add_sequence_alias` / `add_collection_alias` for one alias, `load_sequence_aliases` / `load_collection_aliases` for a file (alias, digest without the `SQ.` prefix). Aliases are stored in the store directory.
- **API:** `dnakit.annotation.refget_alias_register(store, kind, namespace, *, alias=None, digest=None, alias_file=None, verbose=True)`.
- **Input:** Store directory, `kind` (sequence or collection), `namespace`, and `alias` + `digest` or `alias_file`.
- **Output:** `aliases_registered`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_register, refget_collections_list

d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
print(refget_alias_register('store', 'collection', 'demo', alias='v1', digest=d, verbose=False)['aliases_registered'])
print(refget_alias_register('store', 'sequence', 'demo', alias='x', digest='Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', verbose=False)['aliases_registered'])
open('aliases.tsv', 'w').write('second\tYb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN\n')
print(refget_alias_register('store', 'sequence', 'fromfile', alias_file='aliases.tsv', verbose=False)['aliases_registered'])
```

- **Example output:**

```text
1
1
1
```

## 12) `dnakit.annotation.refget_alias_list()`

- **Purpose:** List the alias namespaces of a store, or the aliases of one namespace.
- **Method:** `list_sequence_alias_namespaces` / `list_collection_alias_namespaces`, and `list_sequence_aliases` / `list_collection_aliases` when `namespace` is given.
- **API:** `dnakit.annotation.refget_alias_list(store, kind, *, namespace=None, verbose=True)`.
- **Input:** Store directory, `kind`, optional `namespace`.
- **Output:** `namespaces` or `aliases` (null for an unknown namespace).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_list

print(refget_alias_list('store', 'sequence', verbose=False)['namespaces'])
print(refget_alias_list('store', 'collection', namespace='demo', verbose=False)['aliases'])
print(refget_alias_list('store', 'collection', namespace='none', verbose=False)['aliases'])
```

- **Example output:**

```text
['demo', 'fromfile']
['v1']
None
```

## 13) `dnakit.annotation.refget_alias_reverse_lookup()`

- **Purpose:** Return every (namespace, alias) pair that points at a sequence or collection digest.
- **Method:** `get_aliases_for_sequence` / `get_aliases_for_collection`.
- **API:** `dnakit.annotation.refget_alias_reverse_lookup(store, kind, digest, *, verbose=True)`.
- **Input:** Store directory, `kind`, digest.
- **Output:** `aliases` (namespace, alias).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_reverse_lookup

r = refget_alias_reverse_lookup('store', 'sequence', 'Yb1RIcSYZuLT4uJ6GrN-HikkLFxQxDyN', verbose=False)
print(r['aliases'])
```

- **Example output:**

```text
[{'namespace': 'demo', 'alias': 'x'}, {'namespace': 'fromfile', 'alias': 'second'}]
```

## 14) `dnakit.annotation.refget_alias_remove()`

- **Purpose:** Remove one sequence or collection alias.
- **Method:** `remove_sequence_alias` / `remove_collection_alias`; `removed` is false when the alias did not exist.
- **API:** `dnakit.annotation.refget_alias_remove(store, kind, namespace, alias, *, verbose=True)`.
- **Input:** Store directory, `kind`, `namespace`, `alias`.
- **Output:** `removed`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_remove

print(refget_alias_remove('store', 'sequence', 'demo', 'x', verbose=False)['removed'])
print(refget_alias_remove('store', 'sequence', 'demo', 'x', verbose=False)['removed'])
```

- **Example output:**

```text
True
False
```

## 15) `dnakit.annotation.refget_alias_pull()`

- **Purpose:** Pull alias tables from a remote refget store into a local cache store (network).
- **Method:** `RefgetStore.open_remote(cache, url)` followed by `pull_aliases(namespace, strategy)`; `strategy` is keep-ours, keep-theirs or notify. Ordinary tests never call it.
- **API:** `dnakit.annotation.refget_alias_pull(remote_url, cache_dir, *, namespace=None, strategy='keep-ours', verbose=True)`.
- **Input:** `remote_url`, `cache_dir`, optional `namespace`, `strategy`.
- **Output:** `report` (pulled, skipped, not_found, conflicts).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_alias_pull

try:
    refget_alias_pull('ftp://example.org/store', 'cache', verbose=False)
except Exception as exc:
    print(type(exc).__name__)
```

- **Example output:**

```text
ConfigurationError
```

## 16) `dnakit.annotation.refget_collection_import()`

- **Purpose:** Copy one collection, with its sequences, aliases and FHR metadata, from one local store into another.
- **Method:** `RefgetStore.import_collection(source, digest)`; both stores must be disk-backed.
- **API:** `dnakit.annotation.refget_collection_import(store, source_store, digest, *, verbose=True)`.
- **Input:** Destination store, source store, collection digest.
- **Output:** `imported` (collection metadata).
- **Example code:**

```python
import shutil
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_import, refget_collections_list

shutil.rmtree('copy', ignore_errors=True)
d = sorted(c['digest'] for c in refget_collections_list('store', verbose=False)['collections'])[0]
r = refget_collection_import('copy', 'store', d, verbose=False)
print(r['imported']['digest'] == d, r['imported']['n_sequences'])
```

- **Example output:**

```text
True 2
```

## 17) `dnakit.annotation.refget_orphan_removal_plan()`

- **Purpose:** Preview which sequences would be deleted if a collection were removed with orphan cleanup; nothing is changed.
- **Method:** `RefgetStore.plan_orphan_removal(digest)` (advisory, takes no lock).
- **API:** `dnakit.annotation.refget_orphan_removal_plan(store, digest, *, verbose=True)`.
- **Input:** Store directory, collection digest.
- **Output:** `orphan_sequences` (sorted digests).
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collections_list, refget_orphan_removal_plan

d = refget_collections_list('copy', verbose=False)['collections'][0]['digest']
print(refget_orphan_removal_plan('copy', d, verbose=False)['orphan_sequences'])
```

- **Example output:**

```text
['ou6fkS6aEnM-xOCptFiIOqDyl8yAhTSy', 'z-pMq0F5z9iTxrIjy78TXVEhhix_H5rc']
```

## 18) `dnakit.annotation.refget_collection_remove()`

- **Purpose:** Remove a collection, with its aliases and metadata, from a store; optionally also the sequences no other collection uses.
- **Method:** `RefgetStore.remove_collection(digest, remove_orphan_sequences)`.
- **API:** `dnakit.annotation.refget_collection_remove(store, digest, *, remove_orphan_sequences=False, verbose=True)`.
- **Input:** Store directory, collection digest, `remove_orphan_sequences`.
- **Output:** `removed`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_collection_remove, refget_collections_list

d = refget_collections_list('copy', verbose=False)['collections'][0]['digest']
print(refget_collection_remove('copy', d, remove_orphan_sequences=True, verbose=False)['removed'])
print(refget_collection_remove('copy', d, verbose=False)['removed'])
print(refget_collections_list('copy', verbose=False)['collections'])
```

- **Example output:**

```text
True
False
[]
```

## 19) `dnakit.annotation.refget_write_rgsi()`

- **Purpose:** Digest a FASTA file and write its RGSI index: the collection-level digests followed by one line per sequence.
- **Method:** `digest_fasta` then `SequenceCollection.write_rgsi`. The output file must not exist.
- **API:** `dnakit.annotation.refget_write_rgsi(fasta, output, *, verbose=True)`.
- **Input:** FASTA path, output path.
- **Output:** `collection_digest`, `sequence_count`, `output`.
- **Example code:**

```python
from dnakit.annotation.gtars_p24_refget_functions import refget_write_rgsi

import os
if os.path.exists('a.rgsi'):
    os.remove('a.rgsi')
r = refget_write_rgsi('a.fa', 'a.rgsi', verbose=False)
print(r['collection_digest'], r['sequence_count'])
print(open('a.rgsi').read().replace('\t', '|'))
```

- **Example output:**

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
