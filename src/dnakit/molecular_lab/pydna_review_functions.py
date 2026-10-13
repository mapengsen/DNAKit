"""Review wrappers for pydna sequence, primer, assembly and recombination functions (G18n).

Each function imports pydna (and Biopython) only when it runs. A missing package
raises ``BackendUnavailableError`` with the ``dnakit[bio-assembly]`` install hint.
Positions are 0-based with ``start`` inclusive and ``end`` exclusive. Every function
returns a JSON-safe dictionary with ``status``, ``operation``, the result fields and
``provenance``. Inputs are DNA strings of A, C, G and T; they are upper-cased.
"""

from __future__ import annotations

import importlib.metadata
import random
import re
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-assembly"
_MAX_DNA = 200_000
_MAX_GENOME = 20_000
_MAX_FRAGMENTS = 20
_MAX_PRIMER = 200
_MAX_PRODUCTS = 200
_MAX_CODONS = 10_000
_MAX_TEXT_BYTES = 1_000_000
_MAX_TAIL = 200
_SITE = re.compile(r"[A-Z]+[a-z]+[A-Z]+")
_POLYMERASE_FUNCTIONS = {"taq": "program", "dbd": "dbd_program", "pfu_sso7d": "pfu_sso7d_program"}
_GENOME_EDITS = {
    ("homologous_recombination", "integration"): "homologous_recombination_integration",
    ("homologous_recombination", "excision"): "homologous_recombination_excision",
    (
        "homologous_recombination",
        "excision_or_inversion",
    ): "homologous_recombination_excision_or_inversion",
    ("cre_lox", "integration"): "cre_lox_integration",
    ("cre_lox", "excision"): "cre_lox_excision",
    ("cre_lox", "excision_or_inversion"): "cre_lox_excision_or_inversion",
    ("recombinase", "integration"): "recombinase_integration",
    ("recombinase", "excision"): "recombinase_excision",
    ("recombinase", "excision_or_inversion"): "recombinase_excision_or_inversion",
    ("crispr", "integration"): "crispr_integration",
}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        options = ", ".join(allowed)
        raise ConfigurationError(f"{name} must be one of: {options}.")
    return value


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _dna_list(
    value: object, name: str, *, minimum: int = 1, maximum: int = _MAX_FRAGMENTS
) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of DNA strings.")
    if not minimum <= len(value) <= maximum:
        raise ConfigurationError(f"{name} must contain between {minimum} and {maximum} items.")
    texts = [_dna(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if sum(len(item) for item in texts) > _MAX_DNA:
        raise ConfigurationError(f"{name} exceeds {_MAX_DNA} bases in total.")
    return texts


def _reverse_complement(text: str) -> str:
    return text.translate(_COMPLEMENT)[::-1]


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("pydna")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "pydna", "version": installed, **notes}


def _product_rows(products: Sequence[Any], *, listed: int) -> tuple[list[dict[str, Any]], int]:
    rows = sorted(
        (
            {
                "sequence": str(product.seq).upper(),
                "length": len(product),
                "circular": bool(product.circular),
            }
            for product in products
        ),
        key=lambda row: (row["length"], row["sequence"], row["circular"]),
    )
    return rows[:listed], len(rows)


def _run(action: Any, message: str) -> Any:
    try:
        return action()
    except (AssertionError, AttributeError, IndexError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"{message}: {exc}") from exc


def _dseq(text: str, name: str = "", *, circular: bool = False) -> Any:
    dseq_class = require("pydna.dseqrecord", _EXTRA).Dseqrecord
    return dseq_class(text, name=name, circular=circular)


def dna_molecules_equal(
    sequences: Sequence[str],
    *,
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether DNA strings describe the same double-stranded molecule.

    A sequence equals its reverse complement. With ``circular=True`` every
    rotation of a sequence also counts as the same molecule. At least two
    sequences are needed; the result is ``True`` only when all of them match.
    """

    eq = require("pydna.utils", _EXTRA).eq
    texts = _dna_list(sequences, "sequences", minimum=2)
    flag = _flag(circular, "circular")
    equal = _run(lambda: bool(eq(*texts, circular=flag)), "pydna could not compare the sequences")
    return {
        "status": "completed",
        "operation": "dna_molecules_equal",
        "sequence_count": len(texts),
        "lengths": [len(text) for text in texts],
        "circular": flag,
        "equal": equal,
        "provenance": _provenance(method="pydna.utils.eq"),
    }


def rare_codon_positions(
    sequence: str,
    *,
    organism: str = "sce",
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate the codons that are rare in an organism's codon usage table.

    ``organism`` is ``"sce"`` (yeast) or ``"eco"`` (E. coli), the keys of pydna's
    bundled rare-codon tables. Each hit is ``{"start", "end", "codon"}`` with 0-based
    coordinates. Bases that do not form a complete codon at the end are ignored.
    """

    codon_module = require("pydna.codon", _EXTRA)
    text = _dna(sequence, "sequence", maximum=_MAX_DNA)
    name = _choice(organism, "organism", sorted(codon_module.rare_codons))
    raw = _run(
        lambda: require("pydna.utils", _EXTRA).rarecodons(text, name),
        "pydna could not scan the codons",
    )
    positions = [
        {"start": int(part.start), "end": int(part.stop), "codon": text[part.start : part.stop]}
        for part in raw
    ]
    return {
        "status": "completed",
        "operation": "rare_codon_positions",
        "organism": name,
        "codon_count": len(text) // 3,
        "rare_codon_count": len(positions),
        "positions": positions,
        "provenance": _provenance(method="pydna.utils.rarecodons"),
    }


def smallest_circular_rotation(sequence: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the lexicographically smallest rotation of a circular DNA sequence.

    The rotation is a canonical form for circular molecules: two plasmid maps
    that differ only by their starting base give the same result.
    """

    text = _dna(sequence, "sequence")
    rotated = _run(
        lambda: require("pydna.utils", _EXTRA).smallest_rotation(text),
        "pydna could not compute the rotation",
    )
    return {
        "status": "completed",
        "operation": "smallest_circular_rotation",
        "sequence": str(rotated).upper(),
        "length": len(text),
        "provenance": _provenance(method="pydna.utils.smallest_rotation"),
    }


def generate_random_orf(
    codons: int,
    *,
    max_codons: int | None = None,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a random open reading frame: ATG, random sense codons, then a stop.

    ``codons`` counts every codon including the start and the stop, so the result
    has ``3 * codons`` bases. With ``max_codons`` the length is drawn from a
    triangular distribution between the two values. The random seed is set for
    the call only; the global random state is restored afterwards.
    """

    count = _integer(codons, "codons", minimum=2, maximum=_MAX_CODONS)
    upper = None
    if max_codons is not None:
        upper = _integer(max_codons, "max_codons", minimum=count, maximum=_MAX_CODONS)
    value = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    module = require("pydna.utils", _EXTRA)
    saved = random.getstate()
    try:
        random.seed(value)
        text = _run(lambda: module.randomORF(count, upper), "pydna could not build the ORF")
    finally:
        random.setstate(saved)
    text = str(text).upper()
    return {
        "status": "completed",
        "operation": "generate_random_orf",
        "sequence": text,
        "length": len(text),
        "codon_count": len(text) // 3,
        "seed": value,
        "max_codons": upper,
        "provenance": _provenance(method="pydna.utils.randomORF"),
    }


def annotate_primer_binding_sites(
    product: str,
    forward_primer: str,
    reverse_primer: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate where two primers bind the ends of an assembled or amplified product.

    The forward primer must match the start of ``product`` and the reverse primer
    (written 5' to 3') must match its end as a reverse complement. Each binding site
    becomes a ``primer_bind`` feature: the forward one on the top strand, the
    reverse one on the bottom strand. Features are sorted by start, then end.
    """

    assembly = require("pydna.assembly2", _EXTRA)
    text = _dna(product, "product")
    fwd = _dna(forward_primer, "forward_primer", maximum=_MAX_PRIMER)
    rev = _dna(reverse_primer, "reverse_primer", maximum=_MAX_PRIMER)
    if not text.startswith(fwd):
        raise ConfigurationError("forward_primer must match the start of product.")
    if not text.endswith(_reverse_complement(rev)):
        raise ConfigurationError("reverse_primer must match the end of product.")
    annotated = _run(
        lambda: assembly.annotate_primer_binding_sites(
            _dseq(text, "product"),
            [
                _dseq(fwd, "forward_primer"),
                _dseq(text, "product"),
                _dseq(rev, "reverse_primer"),
            ],
        ),
        "pydna could not annotate the primer sites",
    )
    features = []
    for feature in annotated.features:
        labels = feature.qualifiers.get("label", [""])
        features.append(
            {
                "type": str(feature.type),
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "strand": int(feature.location.strand or 0),
                "label": str(labels[0]),
            }
        )
    features.sort(key=lambda row: (row["start"], row["end"], row["strand"], row["label"]))
    return {
        "status": "completed",
        "operation": "annotate_primer_binding_sites",
        "product_length": len(text),
        "forward_length": len(fwd),
        "reverse_length": len(rev),
        "feature_count": len(features),
        "features": features,
        "provenance": _provenance(method="pydna.assembly2.annotate_primer_binding_sites"),
    }


def suggest_pcr_program(
    template: str,
    forward_primer: str,
    reverse_primer: str,
    *,
    polymerase: str = "taq",
    min_annealing: int = 14,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate PCR with a primer pair and return pydna's suggested cycling programme.

    ``polymerase`` picks the Tm and annealing rules: ``"taq"`` (default Taq or
    similar), ``"dbd"`` (DBD-type polymerases) or ``"pfu_sso7d"``. Each primer must
    anneal over at least ``min_annealing`` bases. The programme is a text block;
    ``program_lines`` holds the same block split by line.
    """

    choice = _choice(polymerase, "polymerase", sorted(_POLYMERASE_FUNCTIONS))
    text = _dna(template, "template")
    fwd = _dna(forward_primer, "forward_primer", maximum=_MAX_PRIMER)
    rev = _dna(reverse_primer, "reverse_primer", maximum=_MAX_PRIMER)
    limit = _integer(min_annealing, "min_annealing", minimum=8, maximum=60)
    amplify = require("pydna.amplify", _EXTRA)
    primer_class = require("pydna.primer", _EXTRA).Primer
    tm_module = require("pydna.tm", _EXTRA)
    amplicon = _run(
        lambda: amplify.pcr(primer_class(fwd), primer_class(rev), _dseq(text), limit=limit),
        "pydna could not amplify the template",
    )
    function = _POLYMERASE_FUNCTIONS[choice]
    program = _run(
        lambda: getattr(tm_module, function)(amplicon),
        "pydna could not write the programme",
    )
    return {
        "status": "completed",
        "operation": "suggest_pcr_program",
        "polymerase": choice,
        "forward_primer": fwd,
        "reverse_primer": rev,
        "min_annealing": limit,
        "product_length": len(amplicon),
        "program": program,
        "program_lines": program.splitlines(),
        "provenance": _provenance(method=f"pydna.tm.{function}"),
    }


def ligate_blunt_fragments(
    fragments: Sequence[str],
    *,
    allow_partial_overlap: bool = False,
    circular_only: bool = False,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate the products of ligating blunt, double-stranded DNA fragments.

    Each entry of ``fragments`` is a blunt linear DNA string; sticky ends from a
    restriction digest are not represented. ``allow_partial_overlap`` also joins
    fragments whose ends overlap partially. Products are sorted by length and
    sequence; at most ``max_products`` are listed and ``product_count`` gives the total.
    """

    assembly = require("pydna.assembly2", _EXTRA)
    texts = _dna_list(fragments, "fragments", minimum=2)
    partial = _flag(allow_partial_overlap, "allow_partial_overlap")
    circular = _flag(circular_only, "circular_only")
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    products = _run(
        lambda: assembly.ligation_assembly(
            [_dseq(text) for text in texts],
            allow_blunt=True,
            allow_partial_overlap=partial,
            circular_only=circular,
        ),
        "pydna could not ligate the fragments",
    )
    rows, total = _product_rows(products, listed=listed)
    return {
        "status": "completed",
        "operation": "ligate_blunt_fragments",
        "fragment_count": len(texts),
        "allow_partial_overlap": partial,
        "circular_only": circular,
        "product_count": total,
        "products": rows,
        "provenance": _provenance(method="pydna.assembly2.ligation_assembly"),
    }


def gateway_recombination_products(
    fragments: Sequence[str],
    *,
    reaction_type: str = "BP",
    greedy: bool = False,
    circular_only: bool = False,
    multi_site_only: bool = False,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate Gateway BP or LR recombination products of DNA fragments.

    ``reaction_type`` is ``"BP"`` (attB with attP) or ``"LR"`` (attL with attR).
    Fragments must carry the Gateway sites in their DNA; ``greedy`` uses the wider
    site definitions. Products are sorted by length and sequence.
    """

    module = require("pydna.assembly2", _EXTRA)
    texts = _dna_list(fragments, "fragments", minimum=2)
    reaction = _choice(reaction_type, "reaction_type", ("BP", "LR"))
    flags = (
        _flag(greedy, "greedy"),
        _flag(circular_only, "circular_only"),
        _flag(multi_site_only, "multi_site_only"),
    )
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    products = _run(
        lambda: module.gateway_assembly(
            [_dseq(text) for text in texts],
            reaction,
            greedy=flags[0],
            circular_only=flags[1],
            multi_site_only=flags[2],
        ),
        "pydna could not run the Gateway reaction",
    )
    rows, total = _product_rows(products, listed=listed)
    return {
        "status": "completed",
        "operation": "gateway_recombination_products",
        "reaction_type": reaction,
        "site_set": "greedy" if flags[0] else "conservative",
        "fragment_count": len(texts),
        "product_count": total,
        "products": rows[:listed],
        "provenance": _provenance(method="pydna.assembly2.gateway_assembly"),
    }


def _recombinase(sites: object) -> Any:
    if isinstance(sites, (str, bytes)) or not isinstance(sites, Sequence) or len(sites) != 2:
        raise ConfigurationError("recombinase_sites must be a pair of site strings.")
    pair = []
    for index, site in enumerate(sites):
        if not isinstance(site, str) or not _SITE.fullmatch(site.strip()):
            raise ConfigurationError(
                f"recombinase_sites[{index}] must look like ATGCCCTAAaaTT: upper-case flanks "
                "around a lower-case homology core."
            )
        pair.append(site.strip())
    recombinase = require("pydna.recombinase", _EXTRA).Recombinase
    return _run(lambda: recombinase(pair[0], pair[1]), "pydna rejected the recombinase sites")


def simulate_genome_edit(
    genome: str,
    *,
    mechanism: str,
    operation: str,
    inserts: Sequence[str] = (),
    guides: Sequence[str] = (),
    recombinase_sites: Sequence[str] | None = None,
    genome_circular: bool = False,
    inserts_circular: bool = False,
    limit: int = 40,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate one genome-editing event on a DNA sequence with pydna.

    ``mechanism`` is ``"homologous_recombination"``, ``"cre_lox"``, ``"recombinase"``
    or ``"crispr"``. ``operation`` is ``"integration"`` (needs ``inserts``),
    ``"excision"`` or ``"excision_or_inversion"`` (no inserts). CRISPR supports only
    integration and needs ``guides``, 20-base spacers that cut the genome. Recombinase
    needs ``recombinase_sites``: two sites such as ``("ATGCCCTAAaaTT", "AAaaTTTTTTTCCCT")``,
    with a lower-case homology core. CRE-lox integration needs circular inserts, so set
    ``inserts_circular=True``. ``limit`` is the minimum homology length and is used only
    for homologous recombination and CRISPR. Products are sorted by length and sequence.
    """

    mech = _choice(mechanism, "mechanism", sorted({key[0] for key in _GENOME_EDITS}))
    op = _choice(
        operation,
        "operation",
        sorted({key[1] for key in _GENOME_EDITS if key[0] == mech}),
    )
    function_name = _GENOME_EDITS[(mech, op)]
    text = _dna(genome, "genome", maximum=_MAX_GENOME)
    insert_texts = _dna_list(inserts, "inserts", minimum=0) if inserts else []
    guide_texts = _dna_list(guides, "guides", minimum=0) if guides else []
    needs_inserts = op == "integration"
    if needs_inserts and not insert_texts:
        raise ConfigurationError("integration needs at least one insert.")
    if not needs_inserts and insert_texts:
        raise ConfigurationError("excision and inversion take no inserts.")
    if mech == "crispr" and not guide_texts:
        raise ConfigurationError("crispr needs at least one guide spacer.")
    if mech != "crispr" and guide_texts:
        raise ConfigurationError("guides apply only to the crispr mechanism.")
    if (mech == "recombinase") != (recombinase_sites is not None):
        raise ConfigurationError("recombinase_sites is required for, and only for, recombinase.")
    genome_flag = _flag(genome_circular, "genome_circular")
    insert_flag = _flag(inserts_circular, "inserts_circular")
    uses_limit = mech in {"homologous_recombination", "crispr"}
    minimum = _integer(limit, "limit", minimum=1, maximum=2_000)
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)

    module = require("pydna.assembly2", _EXTRA)
    function = getattr(module, function_name)
    genome_record = _dseq(text, "genome", circular=genome_flag)
    records = [_dseq(item, "insert", circular=insert_flag) for item in insert_texts]
    args: list[Any] = [genome_record]
    if needs_inserts:
        args.append(records)
    if mech == "recombinase":
        args.append(_recombinase(recombinase_sites))
    if mech == "crispr":
        primer_class = require("pydna.primer", _EXTRA).Primer
        args.append([primer_class(item) for item in guide_texts])
    if uses_limit:
        args.append(minimum)
    products = _run(lambda: function(*args), "pydna could not simulate the edit")
    rows, total = _product_rows(products, listed=listed)
    return {
        "status": "completed",
        "operation": "simulate_genome_edit",
        "mechanism": mech,
        "edit": op,
        "genome_length": len(text),
        "insert_count": len(insert_texts),
        "guide_count": len(guide_texts),
        "limit": minimum if uses_limit else None,
        "product_count": total,
        "products": rows,
        "provenance": _provenance(method=f"pydna.assembly2.{function_name}"),
    }


def assemble_recombinase_fragments(
    fragments: Sequence[str],
    *,
    recombinase_sites: Sequence[str],
    circular_only: bool = False,
    max_products: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Enumerate recombinase assembly products of DNA fragments.

    ``recombinase_sites`` is the site pair described in ``simulate_genome_edit``.
    Each fragment carries the site in its DNA; the sites join fragments into one
    product. Products are sorted by length and sequence.
    """

    module = require("pydna.assembly2", _EXTRA)
    texts = _dna_list(fragments, "fragments", minimum=2)
    recombinase = _recombinase(recombinase_sites)
    circular = _flag(circular_only, "circular_only")
    listed = _integer(max_products, "max_products", minimum=1, maximum=_MAX_PRODUCTS)
    products = _run(
        lambda: module.recombinase_assembly(
            [_dseq(text) for text in texts], recombinase, circular_only=circular
        ),
        "pydna could not run the recombinase assembly",
    )
    rows, total = _product_rows(products, listed=listed)
    return {
        "status": "completed",
        "operation": "assemble_recombinase_fragments",
        "fragment_count": len(texts),
        "circular_only": circular,
        "product_count": total,
        "products": rows,
        "provenance": _provenance(method="pydna.assembly2.recombinase_assembly"),
    }


def extract_protospacers(
    guide_construct: str,
    *,
    cas: str = "cas9",
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract Cas9 protospacers that sit directly before the guide scaffold.

    The construct is searched on both strands. A circular construct is searched
    across its origin. ``cas`` is ``"cas9"``, the only nuclease pydna defines.
    """

    crispr = require("pydna.crispr", _EXTRA)
    text = _dna(guide_construct, "guide_construct")
    name = _choice(cas, "cas", ("cas9",))
    flag = _flag(circular, "circular")
    nuclease = crispr.cas9
    found = _run(
        lambda: crispr.protospacer(_dseq(text, "guide_construct", circular=flag), cas=nuclease),
        "pydna could not scan the construct",
    )
    protospacers = sorted(str(item).upper() for item in found)
    return {
        "status": "completed",
        "operation": "extract_protospacers",
        "cas": name,
        "protospacer_length": nuclease.size,
        "scaffold": nuclease.scaffold,
        "circular": flag,
        "protospacer_count": len(protospacers),
        "protospacers": protospacers,
        "provenance": _provenance(method="pydna.crispr.protospacer"),
    }


def parse_sequence_text(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse EMBL, GenBank or FASTA text into DNA records with pydna.

    The text is read in memory, so no file is needed. Every record keeps its
    input order and reports its name, length, topology and upper-case sequence.
    """

    parser = require("pydna.parsers", _EXTRA)
    if not isinstance(text, str) or not text.strip():
        raise ConfigurationError("text must be non-empty EMBL, GenBank or FASTA text.")
    if len(text.encode("utf-8")) > _MAX_TEXT_BYTES:
        raise ConfigurationError(f"text exceeds {_MAX_TEXT_BYTES} bytes.")
    parsed = _run(lambda: parser.embl_gb_fasta(text), "pydna could not parse the text")
    if not parsed:
        raise ConfigurationError("no sequence records were found in the text.")
    records = [
        {
            "name": str(getattr(item, "name", "")),
            "length": len(item),
            "topology": "circular" if bool(getattr(item, "circular", False)) else "linear",
            "sequence": str(item.seq).upper(),
        }
        for item in parsed
    ]
    return {
        "status": "completed",
        "operation": "parse_sequence_text",
        "record_count": len(records),
        "records": records,
        "provenance": _provenance(method="pydna.parsers.embl_gb_fasta"),
    }


def design_user_cloning_primers(
    fragments: Sequence[Mapping[str, str]],
    *,
    max_overlap: int = 15,
    min_overlap: int = 4,
    max_tail: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design USER-cloning primer tails that give each fragment overlaps with its neighbours.

    Each fragment is a mapping with ``template``, ``forward_primer`` and
    ``reverse_primer`` DNA strings. The primers must amplify the template. The
    result gives the tailed primers in the input order; uracil in the tails is kept
    as ``U``, as pydna writes it.
    """

    amplify = require("pydna.amplify", _EXTRA)
    primer_class = require("pydna.primer", _EXTRA).Primer
    design = require("pydna.design", _EXTRA).user_assembly_design
    if isinstance(fragments, (str, bytes)) or not isinstance(fragments, Sequence):
        raise ConfigurationError("fragments must be a list of mappings.")
    if not 2 <= len(fragments) <= _MAX_FRAGMENTS:
        raise ConfigurationError(f"fragments must contain between 2 and {_MAX_FRAGMENTS} items.")
    amplicons = []
    for index, item in enumerate(fragments):
        if not isinstance(item, Mapping) or set(item) != {
            "template",
            "forward_primer",
            "reverse_primer",
        }:
            raise ConfigurationError(
                f"fragments[{index}] needs template, forward_primer and reverse_primer."
            )
        template = _dna(item["template"], f"fragments[{index}].template")
        forward = _dna(item["forward_primer"], f"fragments[{index}].forward_primer", maximum=200)
        reverse = _dna(item["reverse_primer"], f"fragments[{index}].reverse_primer", maximum=200)
        amplicons.append(
            _run(
                lambda t=template, f=forward, r=reverse: amplify.pcr(
                    primer_class(f), primer_class(r), _dseq(t), limit=14
                ),
                f"pydna could not amplify fragments[{index}]",
            )
        )
    overlap_high = _integer(max_overlap, "max_overlap", minimum=4, maximum=60)
    overlap_low = _integer(min_overlap, "min_overlap", minimum=1, maximum=overlap_high)
    tail = _integer(max_tail, "max_tail", minimum=0, maximum=_MAX_TAIL)
    designed = _run(
        lambda: design(
            amplicons,
            max_overlap=overlap_high,
            min_overlap=overlap_low,
            max_tail=tail,
        ),
        "pydna could not design the tails",
    )
    rows = [
        {
            "forward_primer": str(item.forward_primer.seq),
            "reverse_primer": str(item.reverse_primer.seq),
            "product_length": len(item),
        }
        for item in designed
    ]
    return {
        "status": "completed",
        "operation": "design_user_cloning_primers",
        "fragment_count": len(rows),
        "max_overlap": overlap_high,
        "min_overlap": overlap_low,
        "max_tail": tail,
        "fragments": rows,
        "provenance": _provenance(method="pydna.design.user_assembly_design"),
    }


__all__ = [
    "annotate_primer_binding_sites",
    "assemble_recombinase_fragments",
    "design_user_cloning_primers",
    "dna_molecules_equal",
    "extract_protospacers",
    "gateway_recombination_products",
    "generate_random_orf",
    "ligate_blunt_fragments",
    "parse_sequence_text",
    "rare_codon_positions",
    "simulate_genome_edit",
    "smallest_circular_rotation",
    "suggest_pcr_program",
]

_register_public_verbose(globals(), tuple(__all__))
