# Biopython restriction analysis, CAPS, BioSQL and genome figures

This document covers wrappers over [Biopython](https://biopython.org) 1.88 in three areas: `Bio.Restriction` and `Bio.CAPS` in `dnakit.sequence` (REBASE enzyme data shipped with Biopython; positions are 1-based, the first base after the top-strand cut), `BioSQL` in `dnakit.databases` (records in a SQLite, MySQL or PostgreSQL database built from the BioSQL schema) and `Bio.Graphics` in `dnakit.visualization` (genome diagram cross-links and chromosome ideograms drawn with ReportLab). Install with `python -m pip install "dnakit[bio-formats,reports]"`; packages are imported only when a function runs. The BioSQL examples need the BioSQL SQLite schema file, which Biopython does not ship and DNAKit does not redistribute: download `biosqldb-sqlite.sql` from https://github.com/biosql/biosql and set the environment variable `DNAKIT_BIOSQL_SQLITE_SCHEMA` to its path. The examples import from the modules `dnakit.sequence.biopython_p19_restriction_functions`, `dnakit.databases.biopython_p19_biosql_functions` and `dnakit.visualization.biopython_p19_graphics_functions`.

## Example data

`records.gb` is a two-record GenBank file (`REC1`, 30 bases, and `REC2`, 20 bases, each with one `CDS` feature) used by the BioSQL examples. The restriction examples use the 38-base sequence `GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC`, which contains two EcoRI sites, one BamHI, SmaI and KpnI site each and no HindIII site.

**Example file:** `records.gb`

```text
LOCUS       REC1                      30 bp    DNA              UNK 01-JAN-1980
DEFINITION  test REC1.
ACCESSION   ACC0
VERSION     ACC0.1
KEYWORDS    .
SOURCE      .
  ORGANISM  Testus x
            .
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="g0"
ORIGIN
        1 atggctaaag gtttcccaga actggtttaa
//
LOCUS       REC2                      20 bp    DNA              UNK 01-JAN-1980
DEFINITION  test REC2.
ACCESSION   ACC1
VERSION     ACC1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  Testus x
            .
FEATURES             Location/Qualifiers
     CDS             1..9
                     /gene="g1"
ORIGIN
        1 ggggccccaa aattttacgt
//
```

## 1) `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzyme_properties()`

- **Purpose:** Describe restriction enzymes: recognition site, cut pattern, end type, overhang, frequency and suppliers.
- **Method:** Reads the REBASE data shipped with Biopython (`Bio.Restriction.<Enzyme>`): `elucidate()` gives the cut pattern (`^` top-strand cut, `_` bottom-strand cut), `is_blunt`, `is_5overhang` and `is_3overhang` give the end type, `ovhg` and `ovhgseq` the overhang, `frequency()` the number of random bases per expected site, and `suppl`, `supplier_list` and `is_comm` the suppliers. This covers the end-type, cut-pattern, site-frequency and supplier look-ups of the Biopython enzyme classes.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzyme_properties(enzymes)`.
- **Input:** A list of enzyme names (REBASE names such as `EcoRI`).
- **Output:** One record per enzyme: `site`, `size`, `elucidate`, `end_type`, `overhang`, `overhang_sequence`, `palindromic`, `frequency`, `commercial`, `supplier_codes` and `suppliers`.
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_enzyme_properties

for item in restriction_enzyme_properties(["EcoRI", "SmaI", "KpnI", "BsaI"])["enzymes"]:
    print(item["name"], item["site"], item["elucidate"], item["end_type"], repr(item["overhang_sequence"]), item["palindromic"], item["frequency"], item["commercial"])
print(restriction_enzyme_properties(["EcoRI"])["enzymes"][0]["suppliers"][:3])
```

- **Example output:**

```text
EcoRI GAATTC G^AATT_C 5_overhang 'AATT' True 4096.0 True
SmaI CCCGGG CCC^_GGG blunt '' True 4096.0 True
KpnI GGTACC G_GTAC^C 3_overhang 'GTAC' True 4096.0 True
BsaI GGTCTC GGTCTCN^NNNN_N 5_overhang 'NNNN' False 4096.0 True
['Thermo Fisher Scientific', 'Minotech Biotechnology', 'SibEnzyme Ltd.']
```

## 2) `dnakit.sequence.biopython_p19_restriction_functions.restriction_related_enzymes()`

- **Purpose:** Find the enzymes related to one enzyme: compatible ends, neoschizomers, equischizomers or isoschizomers.
- **Method:** Calls `Bio.Restriction.<Enzyme>.<relation>`. `compatible_end` lists enzymes whose ends can be ligated to this enzyme's ends, `neoschizomers` have the same site but a different cut position, `equischizomers` the same site and the same cut, and `isoschizomers` are both. The search runs over all REBASE enzymes unless `among` restricts it.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.restriction_related_enzymes(enzyme, relation='isoschizomers', *, among=None)`.
- **Input:** An enzyme name, `relation` and an optional list of enzyme names `among`.
- **Output:** `site`, `count` and the sorted `enzymes`.
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_related_enzymes

for enzyme, relation in (("KpnI", "isoschizomers"), ("KpnI", "equischizomers"), ("SmaI", "neoschizomers"), ("EcoRI", "compatible_end")):
    result = restriction_related_enzymes(enzyme, relation)
    print(enzyme, relation, result["count"], result["enzymes"][:4])
print(restriction_related_enzymes("SmaI", "neoschizomers", among=["XmaI", "SmaI", "KpnI"])["enzymes"])
```

- **Example output:**

```text
KpnI isoschizomers 2 ['Acc65I', 'Asp718I']
KpnI equischizomers 0 []
SmaI neoschizomers 3 ['Cfr9I', 'TspMI', 'XmaI']
EcoRI compatible_end 51 ['AarI', 'Acc36I', 'AceIII', 'AcsI']
['XmaI']
```

## 3) `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzymes_by_supplier()`

- **Purpose:** List the enzymes sold by one supplier, or list the suppliers.
- **Method:** Uses `Bio.Restriction.RestrictionBatch.add_supplier` and `suppl_codes`. A supplier is given by its one-letter REBASE code (`N` = New England Biolabs) or its full name. Without `supplier` the table of codes and names with the number of enzymes each supplier sells is returned.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.restriction_enzymes_by_supplier(supplier=None)`.
- **Input:** Optional `supplier` code or name.
- **Output:** With a supplier: `supplier_code`, `supplier`, `enzyme_count` and `enzymes`; without: the `suppliers` table.
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_enzymes_by_supplier

table = restriction_enzymes_by_supplier()["suppliers"]
print(len(table), table[0])
result = restriction_enzymes_by_supplier("N")
print(result["supplier"], result["enzyme_count"], result["enzymes"][:4])
print(restriction_enzymes_by_supplier("New England Biolabs")["enzyme_count"])
```

- **Example output:**

```text
15 {'code': 'B', 'name': 'Thermo Fisher Scientific', 'enzyme_count': 211}
New England Biolabs 238 ['AatII', 'AbaSI', 'Acc65I', 'AccI']
238
```

## 4) `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_enzymes()`

- **Purpose:** Select enzymes by property: end type, availability, supplier, palindromic site or site size.
- **Method:** Applies the predicates `is_blunt`, `is_5overhang`, `is_3overhang`, `is_comm`, `suppl`, `is_palindromic`, `size` and `is_defined` of the Biopython enzyme classes to `enzymes` (default: every enzyme). All given criteria must hold. Biopython 1.88's own `RestrictionBatch.lambdasplit` and `split` return empty batches, so the predicates are applied directly.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_enzymes(enzymes=None, *, end_type=None, commercial=None, supplier=None, palindromic=None, min_site_size=None, max_site_size=None, defined_site=None)`.
- **Input:** Optional `enzymes`, `end_type` (`blunt`, `5_overhang`, `3_overhang`), `commercial`, `supplier` code, `palindromic`, `min_site_size`, `max_site_size` and `defined_site`.
- **Output:** `candidate_count`, `count` and the sorted `enzymes`.
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import filter_restriction_enzymes

blunt = filter_restriction_enzymes(end_type="blunt")
print(blunt["count"], "of", blunt["candidate_count"])
print(filter_restriction_enzymes(["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], end_type="5_overhang", palindromic=True, min_site_size=6)["enzymes"])
sold = filter_restriction_enzymes(commercial=True, supplier="N", max_site_size=4, defined_site=True)
print(sold["count"], sold["enzymes"][:5])
```

- **Example output:**

```text
137 of 1088
['BamHI', 'EcoRI', 'HindIII']
24 ['AciI', 'AluI', 'BfaI', 'BstUI', 'CviAII']
```

## 5) `dnakit.sequence.biopython_p19_restriction_functions.restriction_analysis_report()`

- **Purpose:** Search a sequence for restriction sites and print the result as a list, a count table or a map.
- **Method:** Runs `Bio.Restriction.Analysis` on the chosen enzymes and prints with `print_as` (`list`, `number` or `map`) and `print_that`. Positions are 1-based: the first base after the top-strand cut. Use `linear=False` for a circular molecule.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.restriction_analysis_report(sequence, enzymes, *, linear=True, view='list')`.
- **Input:** A DNA `sequence` (IUPAC letters), a list of enzyme names, `linear` and `view`.
- **Output:** `cuts` (enzyme to positions), `no_cut` and the text `report`.
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import restriction_analysis_report

result = restriction_analysis_report("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], view="map")
print(result["cuts"], result["no_cut"])
tidy = lambda text: "\n".join(line.rstrip() for line in text.splitlines())
print(tidy(result["report"]))
print(tidy(restriction_analysis_report("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], view="number")["report"]))
```

- **Example output:**

```text
{'BamHI': [11], 'EcoRI': [2, 20], 'KpnI': [38], 'SmaI': [27]} ['HindIII']
 2 EcoRI
 |
 |        11 BamHI
 |        |
 |        |        20 EcoRI
 |        |        |
 |        |        |      27 SmaI
 |        |        |      |
 |        |        |      |          38 KpnI
 |        |        |      |          |
GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC
||||||||||||||||||||||||||||||||||||||
CTTAAGTTTCCTAGGTTTCTTAAGGGCCCAAACCATGG
1                                   38
   Enzymes which do not cut the sequence.
HindIII
enzymes which cut 1 times :
BamHI      :  11.
KpnI       :  38.
SmaI       :  27.
enzymes which cut 2 times :
EcoRI      :  2, 20.
   Enzymes which do not cut the sequence.
HindIII
```

## 6) `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_results()`

- **Purpose:** Run a restriction search and keep only the results that match end-type, cut-count, name, site-size or position criteria.
- **Method:** Uses the `Bio.Restriction.Analysis` filters and applies the given criteria one after another: `blunt`, `overhang5`, `overhang3`; `with_N_sites` or `with_number_list`; `with_sites` or `without_site`; the enzyme names and site sizes (applied directly, because the two Biopython 1.88 methods `with_name` and `with_site_size` do not accept them together with an earlier result); and the position filters `only_between`, `between`, `show_only_between`, `only_outside`, `outside` and `do_not_cut` with 1-based `start` and `end`. As in Biopython, enzymes of the chosen end type that do not cut stay in the result with an empty list; `cuts="cutting"` drops them. A step that leaves nothing ends the chain (Biopython would restart from the full result).
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.filter_restriction_results(sequence, enzymes, *, linear=True, end_type=None, cut_count=None, cuts=None, names=None, site_sizes=None, position_filter=None, start=None, end=None)`.
- **Input:** `sequence`, enzyme names, `linear`, and at least one of `end_type`, `cut_count` (a number or a list), `cuts` (`cutting` or `non_cutting`), `names`, `site_sizes`, and `position_filter` with `start` and `end`.
- **Output:** `filters` applied, `match_count` and `results` (enzyme to positions).
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import filter_restriction_results

def show(**criteria):
    print(criteria, filter_restriction_results("GAATTCAAAGGATCCAAAGAATTCCCGGGTTTGGTACC", ["EcoRI", "BamHI", "SmaI", "KpnI", "HindIII"], **criteria)["results"])

show(end_type="blunt")
show(end_type="5_overhang")
show(end_type="5_overhang", cuts="cutting")
show(cut_count=2)
show(cuts="non_cutting")
show(names=["EcoRI", "KpnI"])
show(site_sizes=[6], position_filter="between", start=1, end=20)
show(position_filter="show_only_between", start=15, end=38)
```

- **Example output:**

```text
{'end_type': 'blunt'} {'SmaI': [27]}
{'end_type': '5_overhang'} {'BamHI': [11], 'EcoRI': [2, 20], 'HindIII': []}
{'end_type': '5_overhang', 'cuts': 'cutting'} {'BamHI': [11], 'EcoRI': [2, 20]}
{'cut_count': 2} {'EcoRI': [2, 20]}
{'cuts': 'non_cutting'} {'HindIII': []}
{'names': ['EcoRI', 'KpnI']} {'EcoRI': [2, 20], 'KpnI': [38]}
{'site_sizes': [6], 'position_filter': 'between', 'start': 1, 'end': 20} {'BamHI': [11], 'EcoRI': [2, 20]}
{'position_filter': 'show_only_between', 'start': 15, 'end': 38} {'EcoRI': [20], 'SmaI': [27]}
```

## 7) `dnakit.sequence.biopython_p19_restriction_functions.caps_differential_cut_sites()`

- **Purpose:** Find CAPS marker candidates: restriction sites that cut some aligned sequences but not others.
- **Method:** Calls `Bio.CAPS.CAPSMap` (Konieczny and Ausubel 1993). For every enzyme a site is reported when at least one sequence is cut at that alignment position and at least one is not. Gap characters are not accepted because `Bio.Restriction` cannot search them.
- **API:** `dnakit.sequence.biopython_p19_restriction_functions.caps_differential_cut_sites(alignment, enzymes)`.
- **Input:** A mapping of sequence name to aligned DNA of equal length (IUPAC letters) and a list of enzyme names.
- **Output:** `site_count` and the `sites`: enzyme, 0-based alignment `start` of the cut, and the sequences it cuts (`cuts_in`) or does not cut (`blocked_in`).
- **Example code:**

```python
from dnakit.sequence.biopython_p19_restriction_functions import caps_differential_cut_sites

alignment = {"a": "AAGAATTCGGGAAGGATCCT", "b": "AAGAATTAGGGAAGGATCCT", "c": "AAGAATTCGGGAAGGATCGT"}
result = caps_differential_cut_sites(alignment, ["EcoRI", "BamHI"])
print(result["site_count"])
for site in result["sites"]:
    print(site)
```

- **Example output:**

```text
2
{'enzyme': 'EcoRI', 'start': 3, 'cuts_in': ['a', 'c'], 'blocked_in': ['b']}
{'enzyme': 'BamHI', 'start': 14, 'cuts_in': ['a', 'b'], 'blocked_in': ['c']}
```

## 8) `dnakit.databases.biopython_p19_biosql_functions.biosql_load_records()`

- **Purpose:** Load sequence records from a file into a BioSQL sub-database.
- **Method:** Reads the file with `Bio.SeqIO`, creates the sub-database (`BioSeqDatabase.new_database`) when it is new, loads the records with `BioSeqDatabase.load` (`BioSQL.Loader.DatabaseLoader.load_seqrecord`) and commits; a failure rolls back. No NCBI taxonomy look-up is made. For SQLite, `schema_sql` (the BioSQL file `biosqldb-sqlite.sql`, not shipped with Biopython) creates a new database file; an existing database is never recreated. MySQL (`mysql.connector`) and PostgreSQL (`psycopg2`) databases are reached with `driver`, `host`, `user`, `password` and `port`.
- **API:** `dnakit.databases.biopython_p19_biosql_functions.biosql_load_records(database, namespace, records_path, *, record_format='genbank', driver='sqlite3', schema_sql=None, host=None, user=None, password=None, port=None)`.
- **Input:** The `database` (SQLite file, or the database name for a server), the `namespace`, the `records_path` with `record_format` (default `genbank`), `driver` and optionally `schema_sql`. The examples read the schema path from the environment variable `DNAKIT_BIOSQL_SQLITE_SCHEMA` (download https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-sqlite.sql).
- **Output:** `records_loaded`, `records_in_namespace` and whether the namespace was newly created.
- **Example code:**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
result = biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
print(result["records_loaded"], result["records_in_namespace"], result["created_namespace"])
```

- **Example output:**

```text
2 2 True
```

## 9) `dnakit.databases.biopython_p19_biosql_functions.biosql_get_record()`

- **Purpose:** Retrieve one record from a BioSQL sub-database by identifier.
- **Method:** `id_type` selects the identifier kind: `display_id` (`BioSeqDatabase.get_Seq_by_id`), `accession` (`get_Seq_by_acc`), `version` (`get_Seq_by_ver`), or `primary_id`, `name` and `gi` (`BioSeqDatabase.lookup`). With `output_path` (a new file) the record is also written with `Bio.SeqIO` in `output_format`.
- **API:** `dnakit.databases.biopython_p19_biosql_functions.biosql_get_record(database, namespace, identifier, *, id_type='display_id', include_sequence=True, max_sequence_length=1000000, output_path=None, output_format='genbank', driver='sqlite3', host=None, user=None, password=None, port=None)`.
- **Input:** `database`, `namespace`, `identifier`, `id_type`, `include_sequence`, `max_sequence_length`, an optional `output_path` with `output_format`, and the connection arguments.
- **Output:** `record` with id, name, description, length, sequence, feature summary, annotations and cross-references.
- **Example code:**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_get_record

record = biosql_get_record(database, "demo", "REC1")["record"]
print(record["id"], record["name"], record["length"], record["sequence"], record["feature_types"])
print(biosql_get_record(database, "demo", "ACC1", id_type="accession")["record"]["name"])
print(biosql_get_record(database, "demo", "ACC0.1", id_type="version", include_sequence=False)["record"]["length"])
```

- **Example output:**

```text
ACC0.1 REC1 30 ATGGCTAAAGGTTTCCCAGAACTGGTTTAA ['CDS']
REC2
30
```

## 10) `dnakit.databases.biopython_p19_biosql_functions.biosql_get_subsequence()`

- **Purpose:** Return a slice of a stored sequence without loading the whole record.
- **Method:** Cuts the slice `[start, end)` (0-based, half-open) inside the database with SQL `SUBSTR` (`BioSQL.BioSeqDatabase.Adaptor.get_subseq_as_string`). The record is chosen by `identifier` and `id_type` as in `biosql_get_record`.
- **API:** `dnakit.databases.biopython_p19_biosql_functions.biosql_get_subsequence(database, namespace, identifier, start, end, *, id_type='display_id', driver='sqlite3', host=None, user=None, password=None, port=None)`.
- **Input:** `database`, `namespace`, `identifier`, `start`, `end`, `id_type` and the connection arguments.
- **Output:** The `sequence` slice and the full `sequence_length`.
- **Example code:**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_get_subsequence

result = biosql_get_subsequence(database, "demo", "REC1", 3, 12)
print(result["sequence"], result["sequence_length"])
print(biosql_get_subsequence(database, "demo", "ACC1.1", 16, 20, id_type="version")["sequence"])
```

- **Example output:**

```text
GCTAAAGGT 30
ACGT
```

## 11) `dnakit.databases.biopython_p19_biosql_functions.biosql_list_entries()`

- **Purpose:** List the sub-databases of a BioSQL database, or the entries of one sub-database.
- **Method:** Without `namespace` it lists the sub-database names with their entry counts (`Adaptor.list_biodatabase_names`). With `namespace` it lists the internal entry ids and the display ids (`Adaptor.list_bioentry_ids`, `list_bioentry_display_ids`), at most `max_entries`.
- **API:** `dnakit.databases.biopython_p19_biosql_functions.biosql_list_entries(database, *, namespace=None, max_entries=10000, driver='sqlite3', host=None, user=None, password=None, port=None)`.
- **Input:** `database`, optional `namespace`, `max_entries` and the connection arguments.
- **Output:** `namespaces` with counts, or `entry_count`, `entry_ids` and `display_ids`.
- **Example code:**

```python
import os, tempfile
from dnakit.databases.biopython_p19_biosql_functions import biosql_load_records

database = os.path.join(tempfile.mkdtemp(), "bio.db")
biosql_load_records(database, "demo", "records.gb", schema_sql=os.environ["DNAKIT_BIOSQL_SQLITE_SCHEMA"])
from dnakit.databases.biopython_p19_biosql_functions import biosql_list_entries

print(biosql_list_entries(database)["namespaces"])
entries = biosql_list_entries(database, namespace="demo")
print(entries["entry_count"], entries["entry_ids"], entries["display_ids"])
```

- **Example output:**

```text
[{'name': 'demo', 'entry_count': 2}]
2 [1, 2] ['REC1', 'REC2']
```

## 12) `dnakit.visualization.biopython_p19_graphics_functions.genome_diagram_crosslinks()`

- **Purpose:** Draw genome tracks as stacked bands joined by cross-link ribbons.
- **Method:** Builds a `Bio.Graphics.GenomeDiagram.Diagram` with one track per genome (odd track levels leave room for the ribbons), adds the features as boxes or arrows (by strand), joins regions of different tracks with `CrossLink` (`flip=True` draws an inverted ribbon for a reversed match) and writes the figure. All tracks use the scale of the longest. The output is `pdf`, `svg` or `eps`; `png` needs the ReportLab Cairo backend `rlPyCairo`. The output file must not exist.
- **API:** `dnakit.visualization.biopython_p19_graphics_functions.genome_diagram_crosslinks(tracks, links, output_path, *, image_format='pdf', title=None, width_cm=24.0, track_height_cm=3.0)`.
- **Input:** `tracks` (`name`, `length`, `features` with `start`, `end`, `strand`, `label`, `color`), `links` (`track_a`, `start_a`, `end_a`, `track_b`, `start_b`, `end_b`, `color`, `flip`; tracks by name or index), a new `output_path`, `image_format`, `title`, `width_cm` and `track_height_cm`. Coordinates are 0-based and half-open; colours are names such as `red` or `#RRGGBB`.
- **Output:** The `output` path, `image_format`, `bytes` and the drawn `tracks` and `links`.
- **Example code:**

```python
import os, tempfile
from dnakit.visualization.biopython_p19_graphics_functions import genome_diagram_crosslinks

tracks = [
    {"name": "genomeA", "length": 1000, "features": [
        {"start": 100, "end": 300, "strand": 1, "label": "geneA", "color": "steelblue"},
        {"start": 500, "end": 700, "strand": -1, "label": "geneB", "color": "#cc3333"}]},
    {"name": "genomeB", "length": 800, "features": [{"start": 50, "end": 250, "strand": -1, "label": "geneA2"}]},
    {"name": "genomeC", "length": 900},
]
links = [
    {"track_a": "genomeA", "start_a": 100, "end_a": 300, "track_b": "genomeB", "start_b": 50, "end_b": 250, "flip": True, "color": "red"},
    {"track_a": 0, "start_a": 500, "end_a": 700, "track_b": 2, "start_b": 100, "end_b": 300},
]
target = os.path.join(tempfile.mkdtemp(), "links.svg")
result = genome_diagram_crosslinks(tracks, links, target, image_format="svg", title="demo")
print(result["image_format"], result["bytes"] > 1000, [t["name"] for t in result["tracks"]], len(result["links"]))
text = open(target).read()
print("geneA" in text, "genomeB" in text)
```

- **Example output:**

```text
svg True ['genomeA', 'genomeB', 'genomeC'] 2
True True
```

## 13) `dnakit.visualization.biopython_p19_graphics_functions.chromosome_ideogram()`

- **Purpose:** Draw chromosomes side by side as an ideogram with telomere caps and marked features.
- **Method:** Builds a `Bio.Graphics.BasicChromosome.Organism` with one `Chromosome` per entry: two `TelomereSegment` caps around a `ChromosomeSegment`, or an `AnnotatedChromosomeSegment` when the chromosome has features (strand +1 is drawn on the right, -1 on the left, 0 across). All chromosomes share one scale (the longest). `telomere_length` defaults to 2 percent of the longest chromosome. The output is `pdf`, `svg` or `eps` (`png` needs `rlPyCairo`) and must not exist.
- **API:** `dnakit.visualization.biopython_p19_graphics_functions.chromosome_ideogram(chromosomes, output_path, *, title='Chromosomes', image_format='pdf', telomere_length=None, page_width_cm=24.0, page_height_cm=18.0)`.
- **Input:** `chromosomes` (`name`, `length`, `features` with `start`, `end`, `strand`, `label`, `color`), a new `output_path`, `title`, `image_format`, `telomere_length`, `page_width_cm` and `page_height_cm`.
- **Output:** The `output` path, `image_format`, `bytes`, `telomere_length` and the drawn `chromosomes`.
- **Example code:**

```python
import os, tempfile
from dnakit.visualization.biopython_p19_graphics_functions import chromosome_ideogram

chromosomes = [
    {"name": "chr1", "length": 1000, "features": [
        {"start": 100, "end": 200, "strand": 1, "label": "a", "color": "red"},
        {"start": 400, "end": 420, "strand": -1, "label": "b"}]},
    {"name": "chr2", "length": 700},
]
target = os.path.join(tempfile.mkdtemp(), "ideogram.svg")
result = chromosome_ideogram(chromosomes, target, title="Demo", image_format="svg")
print(result["image_format"], result["bytes"] > 1000, result["telomere_length"], [c["name"] for c in result["chromosomes"]])
text = open(target).read()
print("chr1" in text, "Demo" in text)
```

- **Example output:**

```text
svg True 20 ['chr1', 'chr2']
True True
```
