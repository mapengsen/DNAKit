"""Contract tests for the public third-party coverage registry."""

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.integrations import (
    beta_diversity_scikit_bio,
    describe_third_party_tool,
    describe_unified_capability,
    integration_inventory,
    list_third_party_tools,
    list_unified_capabilities,
    run_unified_capability,
    third_party_coverage,
    third_party_route,
    unified_capability_summary,
)


def test_coverage_preserves_declared_and_recoverable_counts() -> None:
    summary = third_party_coverage(verbose=False)["summary"]
    assert summary["declared_record_count"] == 282
    assert summary["raw_mentions"] == 425
    assert summary["unrecovered_declared_mentions"] == 0
    assert summary["candidate_capabilities"] == 837
    assert summary["capability_matches"] == 32
    assert summary["integration_status"]["pip_adapter"] == 26
    assert summary["integration_status"]["native_facade"] == 62
    assert summary["unique_records"] >= 200


def test_integration_inventory_separates_surface_counts() -> None:
    inventory = integration_inventory(verbose=False)
    assert inventory["python_functions"] >= inventory["agent_compatible_python_functions"]
    assert inventory["mcp_server_recipes"] > 0
    assert inventory["registered_skills"] > 0
    assert inventory["tool_surface_total"] == 3648
    assert inventory["survey_unique_projects"] == 320
    assert inventory["survey_records_by_kind"]["mcp"] == 103
    assert inventory["survey_records_by_integration_status"]["project_evidence_only"] == 44
    assert inventory["survey_records_by_integration_status"]["no_local_match"] == 171
    assert inventory["unified_candidate_capabilities"] == 837
    assert sum(inventory["unified_candidate_statuses"].values()) == 837


def test_tool_query_contains_provenance_and_citation_status() -> None:
    records = list_third_party_tools("biopython", verbose=False)
    assert records
    record = next(item for item in records if item["name"].casefold() == "biopython")
    assert record["source_lines"]
    assert "citation_status" in record
    assert describe_third_party_tool("Biopython", verbose=False)["name"] == record["name"]


def test_route_points_to_existing_analysis_bridge_without_running_it() -> None:
    route = third_party_route("Cas-OFFinder", verbose=False)
    assert route["api"] == "dnakit.integrations.run_analysis"
    assert "cas-offinder" in route["analysis_ids"]
    assert route["provenance"]["survey_lines"]


def test_biopython_route_lists_related_native_facades() -> None:
    route = third_party_route("Biopython", verbose=False)
    assert route["status"] == "native_facade"
    assert "dnakit.bioinformatics.reads.read_fastx" in route["native_functions"]
    assert "not a claim" in route["boundary"]


def test_existing_pip_dependencies_are_recorded_as_bounded_adapters() -> None:
    summary = third_party_coverage(verbose=False)["summary"]
    assert summary["integration_status"]["pip_adapter"] == 26
    pydna = third_party_route("pydna", verbose=False)
    assert pydna["status"] == "native_facade"
    assert "dnakit.bioinformatics.cloning.assemble_duplexes" in pydna["native_functions"]
    pysam = third_party_route("pysam", verbose=False)
    assert pysam["status"] == "native_facade"
    assert "dnakit.bioinformatics.reads.coverage_consensus" in pysam["native_functions"]
    for name, api, doi in (
        ("BWA", "dnakit.ngs.map_reads", "10.1093/bioinformatics/btp324"),
        ("BWA-MEM2", "dnakit.ngs.map_reads", "10.1109/IPDPS.2019.00041"),
        ("IQ-TREE2", "dnakit.phylogeny.infer_ml_tree", "10.1093/molbev/msaa015"),
        ("Bowtie2", "dnakit.ngs.map_reads", "10.1038/nmeth.1923"),
        ("SAMtools", "dnakit.ngs.prepare_bam", "10.1093/gigascience/giab008"),
        ("BCFtools", "dnakit.ngs.call_variants", "10.1093/gigascience/giab008"),
        ("SPAdes", "dnakit.ngs.assemble_reads", "10.1089/cmb.2012.0021"),
        (
            "NCBI BLAST+",
            "dnakit.alignment.build_blast_database / dnakit.alignment.local_blast",
            "10.1186/1471-2105-10-421",
        ),
        ("Minimap2", "dnakit.alignment.map_long_reads", "10.1093/bioinformatics/bty191"),
        ("BEDTools", "dnakit.alignment.intersect_bed", "10.1093/bioinformatics/btq033"),
        ("DIAMOND", "dnakit.alignment.search_diamond_blastp / dnakit.alignment.search_diamond_blastx", "10.1038/s41592-021-01101-x"),
        ("MMseqs2", "dnakit.alignment.search_mmseqs_protein / dnakit.alignment.cluster_mmseqs_protein", "10.1038/nbt.3988"),
        ("Kraken2", "dnakit.microbiome.classify_kraken2", "10.1186/s13059-019-1891-0"),
        ("Prokka", "dnakit.annotation.annotate_prokka", "10.1093/bioinformatics/btu153"),
        ("Bakta", "dnakit.annotation.annotate_bakta", "10.1099/mgen.0.000685"),
        ("CheckM2", "dnakit.annotation.predict_checkm2", "10.1038/s41592-023-01940-w"),
        ("Bracken", "dnakit.microbiome.estimate_bracken_abundance", "10.7717/peerj-cs.104"),
        ("Centrifuge", "dnakit.microbiome.classify_centrifuge", "10.1101/gr.210641.116"),
        ("Kaiju", "dnakit.microbiome.classify_kaiju", "10.1038/ncomms11257"),
        ("Unicycler", "dnakit.assembly.assemble_unicycler", "10.1371/journal.pcbi.1005595"),
        ("GTDB-Tk", "dnakit.microbiome.classify_gtdbtk", "10.1093/bioinformatics/btz848"),
        ("HUMAnN", "dnakit.microbiome.profile_humann", "10.7554/eLife.65088"),
        ("MethylDackel", "dnakit.epigenome.extract_methylation_methyldackel", ""),
        (
            "GEMINI",
            "dnakit.variants.load_gemini_database / dnakit.variants.query_gemini_database",
            "10.1371/journal.pcbi.1003153",
        ),
        (
            "Ensembl VEP",
            "dnakit.variants.annotate_vcf_ensembl_vep",
            "10.1186/s13059-016-0974-4",
        ),
    ):
        route = third_party_route(name, verbose=False)
        assert route["status"] == "native_facade"
        assert route["api"] == api
        if doi:
            assert doi in route["citation_dois"]
    macs3 = third_party_route("MACS3", verbose=False)
    assert macs3["status"] == "native_facade"
    assert macs3["api"] == "dnakit.ngs.call_chip_peaks"
    assert "BSD-3-Clause" in macs3["license_evidence"]
    freebayes = third_party_route("FreeBayes", verbose=False)
    assert freebayes["status"] == "native_facade"
    assert freebayes["api"] == "dnakit.variants.call_freebayes"
    assert any("arXiv:1207.3907" in item for item in freebayes["citation_references"])
    gatk = third_party_route("GATK4（现行源码）", verbose=False)
    assert gatk["status"] == "native_facade"
    assert gatk["api"] == "dnakit.variants.call_gatk_haplotypecaller"
    assert "Apache-2.0" in gatk["license_evidence"][0]
    vcftools = third_party_route("VCFtools", verbose=False)
    assert vcftools["status"] == "native_facade"
    assert vcftools["api"] == "dnakit.variants.filter_vcf_vcftools"
    assert "10.1093/bioinformatics/btr330" in vcftools["citation_dois"]
    fasttree = third_party_route("FastTree", verbose=False)
    assert fasttree["status"] == "native_facade"
    assert fasttree["api"] == "dnakit.phylogeny.infer_fast_tree"
    assert "10.1371/journal.pone.0009490" in fasttree["citation_dois"]
    shasta = third_party_route("Shasta", verbose=False)
    assert shasta["status"] == "native_facade"
    assert shasta["api"] == "dnakit.assembly.assemble_shasta"
    assert "10.1038/s41587-020-0503-6" in shasta["citation_dois"]
    merqury = third_party_route("Merqury", verbose=False)
    assert merqury["status"] == "native_facade"
    assert merqury["api"] == "dnakit.assembly.evaluate_merqury"
    assert "10.1186/s13059-020-02134-9" in merqury["citation_dois"]
    busco = third_party_route("BUSCO", verbose=False)
    assert busco["status"] == "native_facade"
    assert busco["api"] == "dnakit.assembly.run_busco"
    assert "10.1093/molbev/msab199" in busco["citation_dois"]
    deepvariant = third_party_route("DeepVariant", verbose=False)
    assert deepvariant["status"] == "native_facade"
    assert deepvariant["api"] == "dnakit.variants.call_deepvariant"
    assert "10.1038/nbt.4235" in deepvariant["citation_dois"]
    bedops = third_party_route("BEDOPS", verbose=False)
    assert bedops["status"] == "native_facade"
    assert bedops["api"] == "dnakit.alignment.sort_bed_bedops"
    assert "10.1093/bioinformatics/bts277" in bedops["citation_dois"]
    for name, api, doi in (
        ("Flye", "dnakit.assembly.assemble_flye", "10.1038/s41587-019-0072-8"),
        ("hifiasm", "dnakit.assembly.assemble_hifiasm", "10.1038/s41592-020-01056-5"),
        ("MEGAHIT", "dnakit.assembly.assemble_megahit", "10.1093/bioinformatics/btv033"),
        ("cuteSV", "dnakit.variants.call_cutesv", "10.1186/s13059-020-02107-y"),
        ("Sniffles2", "dnakit.variants.call_sniffles2", "10.1038/s41587-023-02024-y"),
        ("DELLY", "dnakit.variants.call_delly", "10.1093/bioinformatics/bts378"),
        ("HTSlib", "dnakit.ngs.compress_index_vcf_htslib", "10.1093/gigascience/giab007"),
        ("seqtk", "dnakit.ngs.sample_sequences_seqtk", ""),
        ("WhatsHap", "dnakit.ngs.phase_variants_whatshap", "10.1101/085050"),
        ("BamTools", "dnakit.ngs.sort_bam_bamtools", "10.1093/bioinformatics/btr174"),
        ("VSEARCH", "dnakit.fasta_fastq.cluster_vsearch", "10.7717/peerj.2584"),
    ):
        route = third_party_route(name, verbose=False)
        assert route["status"] == "native_facade"
        assert route["api"] == api
        if doi:
            assert doi in route["citation_dois"]
    vcfanno = third_party_route("vcfanno", verbose=False)
    assert vcfanno["status"] == "native_facade"
    assert vcfanno["api"] == "dnakit.variants.annotate_vcf_vcfanno"
    assert "10.1186/s13059-016-0973-5" in vcfanno["citation_dois"]
    survivor = third_party_route("SURVIVOR", verbose=False)
    assert survivor["status"] == "native_facade"
    assert survivor["api"] == "dnakit.variants.merge_structural_variants_survivor"
    assert "10.1038/ncomms14061" in survivor["citation_dois"]
    picard = third_party_route("Picard", verbose=False)
    assert picard["status"] == "native_facade"
    assert picard["api"] == "dnakit.ngs.mark_duplicates_picard"
    spoa = third_party_route("SPOA", verbose=False)
    assert spoa["status"] == "native_facade"
    assert spoa["api"] == "dnakit.alignment.consensus_spoa"
    assert "10.1093/bioinformatics/btg109" in spoa["citation_dois"]


def test_optional_pip_tools_have_explicit_adapters_or_external_boundaries() -> None:
    skbio = third_party_route("scikit-bio", verbose=False)
    assert skbio["api"].startswith("dnakit.integrations.beta_diversity_scikit_bio")
    assert skbio["status"] == "optional_backend"
    assert "full namespace" in skbio["boundary"]
    seqme = third_party_route("seqme", verbose=False)
    assert seqme["api"] == "dnakit.integrations.evaluate_seqme"
    assert seqme["status"] == "optional_backend"
    assert "arXiv:2511.04239" in seqme["citation_references"]
    assert "model extras" in seqme["boundary"]
    dnachisel = third_party_route("DNA Chisel", verbose=False)
    assert dnachisel["api"] == "dnakit.integrations.optimize_with_dnachisel"
    assert dnachisel["status"] == "optional_backend"
    assert "optional backend" in dnachisel["boundary"]
    ugene = third_party_route("UGENE", verbose=False)
    assert ugene["api"] == "dnakit.external_adapters.run"
    assert ugene["status"] == "external_facade"
    ifeatureomega = third_party_route("iFeatureOmega", verbose=False)
    assert ifeatureomega["status"] == "external_facade"
    assert third_party_route("ga", verbose=False)["status"] == "external_facade"
    assert third_party_route("dnachisel", verbose=False)["status"] == "optional_backend"
    assert third_party_route("pyfaidx", verbose=False)["status"] == "optional_backend"
    assert third_party_route("vcfpy", verbose=False)["status"] == "optional_backend"
    assert third_party_route("cyvcf2", verbose=False)["status"] == "optional_backend"
    cutadapt = third_party_route("Cutadapt", verbose=False)
    assert cutadapt["api"] == "dnakit.sequencing_qc.run_cutadapt"
    assert cutadapt["status"] == "native_facade"
    assert "10.14806/ej.17.1.200" in cutadapt["citation_dois"]
    fastp = third_party_route("fastp", verbose=False)
    assert fastp["api"] == "dnakit.sequencing_qc.run_fastp"
    assert fastp["status"] == "native_facade"
    assert "10.1002/imt2.70078" in fastp["citation_dois"]
    fastqc = third_party_route("FastQC", verbose=False)
    assert fastqc["api"] == "dnakit.sequencing_qc.run_fastqc"
    assert fastqc["status"] == "native_facade"
    multiqc = third_party_route("MultiQC", verbose=False)
    assert multiqc["api"] == "dnakit.sequencing_qc.run_multiqc"
    assert multiqc["status"] == "native_facade"
    assert "10.1093/bioinformatics/btw354" in multiqc["citation_dois"]
    nanofilt = third_party_route("NanoFilt", verbose=False)
    assert nanofilt["api"] == "dnakit.sequencing_qc.filter_nanopore_reads"
    assert nanofilt["status"] == "native_facade"
    assert "10.1093/bioinformatics/bty149" in nanofilt["citation_dois"]
    filtlong = third_party_route("Filtlong", verbose=False)
    assert filtlong["api"] == "dnakit.sequencing_qc.filter_long_reads"
    assert filtlong["status"] == "native_facade"
    nanoplot = third_party_route("NanoPlot", verbose=False)
    assert nanoplot["api"] == "dnakit.sequencing_qc.summarize_long_reads"
    assert nanoplot["status"] == "native_facade"
    assert "10.1093/bioinformatics/btad311" in nanoplot["citation_dois"]
    porechop = third_party_route("Porechop", verbose=False)
    assert porechop["api"] == "dnakit.sequencing_qc.trim_porechop"
    assert porechop["status"] == "native_facade"
    assert "10.1099/mgen.0.000132" in porechop["citation_dois"]
    trimmomatic = third_party_route("Trimmomatic", verbose=False)
    assert trimmomatic["api"] == "dnakit.sequencing_qc.trim_trimmomatic"
    assert trimmomatic["status"] == "native_facade"
    assert "10.1093/bioinformatics/btu170" in trimmomatic["citation_dois"]
    afterqc = third_party_route("AfterQC", verbose=False)
    assert afterqc["api"] == "dnakit.sequencing_qc.run_afterqc"
    assert afterqc["status"] == "native_facade"
    assert "10.1186/s12859-017-1469-3" in afterqc["citation_dois"]
    seqmagick = third_party_route("SeqMagick", verbose=False)
    assert seqmagick["api"] == "dnakit.fasta_fastq.convert_seqmagick"
    assert seqmagick["status"] == "native_facade"
    assert "GPL-3.0" in seqmagick["license_evidence"]
    seqkit = third_party_route("SeqKit", verbose=False)
    assert "10.1002/imt2.191" in seqkit["citation_dois"]
    biotite = third_party_route("Biotite", verbose=False)
    assert biotite["status"] == "optional_backend"
    assert "10.1186/s12859-018-2367-z" in biotite["citation_dois"]
    edlib = third_party_route("edlib", verbose=False)
    assert edlib["status"] == "optional_backend"
    assert "10.1093/bioinformatics/btw753" in edlib["citation_dois"]
    pyranges = third_party_route("PyRanges", verbose=False)
    assert pyranges["status"] == "optional_backend"
    assert "10.1093/bioinformatics/btz615" in pyranges["citation_dois"]
    pyfastx = third_party_route("pyfastx", verbose=False)
    assert pyfastx["status"] == "optional_backend"
    assert pyfastx["api"] == "dnakit.integrations.summarize_pyfastx"
    assert "10.1093/bib/bbaa368" in pyfastx["citation_dois"]
    dnaio = third_party_route("dnaio", verbose=False)
    assert dnaio["status"] == "optional_backend"
    assert dnaio["api"] == "dnakit.integrations.summarize_dnaio_fastx"
    assert "MIT" in dnaio["license_evidence"]
    assert "10.5281/zenodo.10548864" in dnaio["citation_dois"]
    mappy = third_party_route("mappy", verbose=False)
    assert mappy["status"] == "optional_backend"
    assert mappy["api"] == "dnakit.integrations.map_mappy_sequences"
    assert "10.1093/bioinformatics/bty191" in mappy["citation_dois"]
    pyfaidx = third_party_route("pyfaidx", verbose=False)
    assert "10.7287/peerj.preprints.970v1" in pyfaidx["citation_dois"]
    pybigwig = third_party_route("pyBigWig", verbose=False)
    assert pybigwig["status"] == "optional_backend"
    assert pybigwig["api"] == "dnakit.integrations.summarize_pybigwig_region"
    assert "10.5281/zenodo.45238" in pybigwig["citation_dois"]
    pybedtools = third_party_route("pybedtools", verbose=False)
    assert pybedtools["status"] == "optional_backend"
    assert pybedtools["api"] == "dnakit.integrations.intersect_pybedtools_bed"
    assert "10.1093/bioinformatics/btr539" in pybedtools["citation_dois"]
    gffutils = third_party_route("gffutils", verbose=False)
    assert gffutils["status"] == "optional_backend"
    assert gffutils["api"] == "dnakit.integrations.query_gffutils_database"
    assert "MIT" in gffutils["license_evidence"]
    sbol3 = third_party_route("PySBOL3", verbose=False)
    assert sbol3["status"] == "optional_backend"
    assert sbol3["api"] == "dnakit.integrations.summarize_sbol3_document"
    assert "10.1021/acssynbio.2c00164" in sbol3["citation_dois"]
    scadnano = third_party_route("scadnano", verbose=False)
    assert scadnano["status"] == "optional_backend"
    assert scadnano["api"] == "dnakit.integrations.summarize_scadnano_design"
    assert "MIT" in scadnano["license_evidence"]
    polars_bio = third_party_route("polars-bio", verbose=False)
    assert polars_bio["status"] == "optional_backend"
    assert polars_bio["api"] == "dnakit.integrations.overlap_polars_bio_bed"
    assert "10.1093/bioinformatics/btaf640" in polars_bio["citation_dois"]
    dnacauldron = third_party_route("DnaCauldron", verbose=False)
    assert dnacauldron["status"] == "optional_backend"
    assert dnacauldron["api"] == "dnakit.integrations.simulate_dnacauldron_assembly"
    assert "MIT" in dnacauldron["license_evidence"]
    dna_features_viewer = third_party_route("DnaFeaturesViewer", verbose=False)
    assert dna_features_viewer["status"] == "optional_backend"
    assert dna_features_viewer["api"] == "dnakit.integrations.summarize_dna_features_viewer"
    assert "MIT" in dna_features_viewer["license_evidence"]
    gtars = third_party_route("gtars", verbose=False)
    assert gtars["status"] == "optional_backend"
    assert gtars["api"] == "dnakit.integrations.summarize_gtars_regions"
    assert "BSD-2-Clause" in gtars["license_evidence"]
    geneblocks = third_party_route("GeneBlocks", verbose=False)
    assert geneblocks["status"] == "optional_backend"
    assert geneblocks["api"] == "dnakit.integrations.diff_geneblocks_sequences"
    assert "MIT" in geneblocks["license_evidence"]
    sourmash = third_party_route("sourmash", verbose=False)
    assert sourmash["status"] == "optional_backend"
    assert sourmash["api"] == "dnakit.integrations.compare_sourmash_sequences"
    assert "10.21105/joss.06830" in sourmash["citation_dois"]


def test_unified_external_capability_requires_explicit_side_effect_authorization() -> None:
    capability = describe_unified_capability("CAP-0683", verbose=False)
    assert capability["api"] == "dnakit.sequencing_qc.run_cutadapt"
    assert capability["status"] == "native_facade"
    assert capability["effect"] == "external_process"
    with pytest.raises(ConfigurationError, match="may write files"):
        run_unified_capability(
            "CAP-0683",
            {"input": "reads.fastq", "options": {"output": "trimmed.fastq"}},
            verbose=False,
        )


def test_pyfaidx_and_vcf_adapters_keep_bounded_json_contract(monkeypatch, tmp_path) -> None:
    import dnakit.integrations.pip_tools as pip_tools

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nAACCGGTT\n", encoding="ascii")

    class FakeRecord:
        def __getitem__(self, bounds):
            assert bounds == slice(1, 5)
            return "ACCG"

    class FakeFasta:
        def __init__(self, path, **kwargs):
            assert path == str(reference.resolve())
            assert kwargs["as_raw"] is True

        def __getitem__(self, key):
            assert key == "chr1"
            return FakeRecord()

        def close(self):
            return None

    class FakePyfaidx:
        Fasta = FakeFasta

    class FakeSamples:
        names: tuple[str, ...] = ("sample1",)

    class FakeHeader:
        samples = FakeSamples()

    class FakeVcfRecord:
        CHROM = "chr1"
        POS = 2
        REF = "A"
        ALT: tuple[str, ...] = ("G",)

    class FakeReader:
        header = FakeHeader()

        def __iter__(self):
            return iter([FakeVcfRecord()])

        def close(self):
            return None

    class FakeReaderFactory:
        @staticmethod
        def from_path(path):
            assert path == str(reference.resolve())
            return FakeReader()

    class FakeVcfpy:
        class Reader:
            from_path = FakeReaderFactory.from_path

    class FakeCyvcf2:
        @staticmethod
        def VCF(path):
            assert path == str(reference.resolve())
            return FakeReader()

    def fake_require(module, extra):
        if module == "pyfaidx":
            assert extra == "bio-pyfaidx"
            return FakePyfaidx
        if module == "vcfpy":
            assert extra == "bio-vcfpy"
            return FakeVcfpy
        assert module == "cyvcf2"
        assert extra == "bio-cyvcf2"
        return FakeCyvcf2

    monkeypatch.setattr(pip_tools, "require", fake_require)
    fetched = pip_tools.fetch_pyfaidx_region(reference, "chr1", 1, 5, strand="-", verbose=False)
    assert fetched["sequence"] == "CGGT"
    summary = pip_tools.summarize_vcf_vcfpy(reference, verbose=False)
    assert summary["records"] == 1
    assert summary["variant_kinds"] == {"snp": 1}
    cyvcf2_summary = pip_tools.summarize_vcf_cyvcf2(reference, verbose=False)
    assert cyvcf2_summary["records"] == 1


def test_biotite_edlib_and_pyranges_adapters_keep_json_contract(monkeypatch, tmp_path) -> None:
    import dnakit.integrations.pip_tools as pip_tools

    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t1\t5\nchr1\t4\t8\n", encoding="ascii")

    class FakeAlignment:
        score = 3

        class trace:
            @staticmethod
            def tolist():
                return [[0, 0], [1, -1], [-1, 1]]

    class FakeAlign:
        class SubstitutionMatrix:
            @staticmethod
            def std_nucleotide_matrix():
                return object()

        @staticmethod
        def align_optimal(*args, **kwargs):
            assert kwargs["local"] is False
            assert kwargs["max_number"] == 2
            return [FakeAlignment()]

        @staticmethod
        def get_symbols(alignment):
            assert isinstance(alignment, FakeAlignment)
            return [["A", "C", None], ["A", None, "C"]]

        @staticmethod
        def get_sequence_identity(alignment, mode):
            assert mode == "all"
            return 1 / 3

        @staticmethod
        def write_alignment_to_cigar(alignment, **kwargs):
            assert kwargs["reference_index"] == 1
            assert kwargs["segment_index"] == 0
            assert kwargs["include_terminal_gaps"] is True
            return "1=1I1D"

    class FakeSequence:
        NucleotideSequence = str

    class FakeEdlib:
        @staticmethod
        def align(query, target, **kwargs):
            assert query == "AC"
            assert target == "AG"
            assert kwargs["task"] == "path"
            if kwargs["k"] == 0:
                return {"editDistance": -1, "locations": [], "cigar": None}
            return {"editDistance": 1, "locations": [(0, 1)], "cigar": "1=1X"}

        @staticmethod
        def getNiceAlignment(result, query, target):
            assert result["editDistance"] == 1
            return {"query_aligned": query, "target_aligned": target}

    class FakeFrame:
        def __init__(self, rows):
            self.rows = rows

        def to_dict(self, orient):
            assert orient == "records"
            return self.rows

    class FakeRanges:
        def __init__(self, rows):
            self.df = FakeFrame(rows)

        def merge(self, **kwargs):
            assert kwargs == {"slack": 1, "count": True, "strand": None}
            return FakeRanges([{"Chromosome": "chr1", "Start": 1, "End": 8, "Count": 2}])

    class FakePyRanges:
        @staticmethod
        def read_bed(path, **kwargs):
            assert path == str(bed.resolve())
            assert kwargs["as_df"] is False
            assert kwargs["nrows"] in {2, 11}
            return FakeRanges(
                [
                    {"Chromosome": "chr1", "Start": 1, "End": 5},
                    {"Chromosome": "chr1", "Start": 4, "End": 8},
                ]
            )

    def fake_require(module, extra):
        expected = {
            "biotite.sequence": ("bio-biotite", FakeSequence),
            "biotite.sequence.align": ("bio-biotite", FakeAlign),
            "edlib": ("bio-edlib", FakeEdlib),
            "pyranges": ("bio-pyranges", FakePyRanges),
        }
        required_extra, backend = expected[module]
        assert extra == required_extra
        return backend

    monkeypatch.setattr(pip_tools, "require", fake_require)
    biotite = pip_tools.align_biotite("AC", "AC", max_alignments=2, verbose=False)
    assert biotite["alignments"][0]["cigar"] == "1=1I1D"
    assert biotite["alignments"][0]["spans"] == {"query": [0, 2], "target": [0, 2]}
    edlib = pip_tools.align_edlib("AC", "AG", verbose=False)
    assert edlib["result"]["nice_alignment"]["target_aligned"] == "AG"
    no_edlib_match = pip_tools.align_edlib("AC", "AG", k=0, verbose=False)
    assert no_edlib_match["found"] is False
    assert "nice_alignment" not in no_edlib_match["result"]
    summary = pip_tools.summarize_pyranges_bed(bed, max_rows=10, verbose=False)
    assert summary["row_count"] == 2
    merged = pip_tools.merge_pyranges_bed(bed, slack=1, max_rows=10, verbose=False)
    assert merged["intervals"][0]["Count"] == 2
    from dnakit.exceptions import ConfigurationError

    try:
        pip_tools.summarize_pyranges_bed(bed, max_rows=1, verbose=False)
    except ConfigurationError as exc:
        assert "max_rows=1" in str(exc)
    else:
        raise AssertionError("BED input over max_rows must not be silently truncated")


def test_pyfastx_and_pybigwig_adapters_keep_json_contract(monkeypatch, tmp_path) -> None:
    import dnakit.integrations.pip_tools as pip_tools

    fastx = tmp_path / "reads.fa"
    fastx.write_text(">r1\nACGT\n>r2\nGG\n", encoding="ascii")
    bigwig = tmp_path / "signal.bw"
    bigwig.write_bytes(b"placeholder")
    other_bed = tmp_path / "other.bed"
    other_bed.write_text("chr1\t2\t6\n", encoding="ascii")
    gff_db = tmp_path / "annotations.db"
    gff_db.write_bytes(b"sqlite-placeholder")
    sbol = tmp_path / "design.nt"
    sbol.write_text("placeholder", encoding="ascii")
    scadnano_file = tmp_path / "design.sc"
    scadnano_file.write_text("{}", encoding="ascii")
    polars_source = tmp_path / "source.bed"
    polars_source.write_text("chr1\t1\t5\n", encoding="ascii")
    polars_other = tmp_path / "other.csv"
    polars_other.write_text("chrom,start,end\nchr1,2,6\n", encoding="ascii")
    dna_features_source = tmp_path / "features.gb"
    dna_features_source.write_text("placeholder", encoding="ascii")

    class FakeFastx:
        def __init__(self, path, **kwargs):
            assert path == str(fastx.resolve())
            assert kwargs == {"uppercase": True}

        def __iter__(self):
            return iter([("r1", "ACGT"), ("r2", "GG")])

        def close(self):
            return None

    class FakePyfastx:
        Fastx = FakeFastx

    class FakeBigWigReader:
        def chroms(self):
            return {"chr1": 100}

        def stats(self, chrom, start, end, **kwargs):
            assert (chrom, start, end) == ("chr1", 10, 20)
            assert kwargs == {"type": "mean", "nBins": 2}
            return [0.25, None]

        def close(self):
            return None

    class FakePyBigWig:
        @staticmethod
        def open(path):
            assert path == str(bigwig.resolve())
            return FakeBigWigReader()

    class FakeBedTool:
        def __init__(self, path):
            assert path in {str(fastx.resolve()), str(other_bed.resolve())}

        def intersect(self, other, **kwargs):
            assert isinstance(other, FakeBedTool)
            assert kwargs == {"u": True, "wa": False, "wb": True, "stream": True}
            return iter(["chr1\t1\t5\tchr1\t2\t6\n"])

    class FakePybedtools:
        BedTool = FakeBedTool

    class FakeFeature:
        def __init__(self):
            self.id = "gene1"
            self.seqid = "chr1"
            self.source = "test"
            self.featuretype = "gene"
            self.start = 1
            self.end = 10
            self.score = None
            self.strand = "+"
            self.frame = None
            self.attributes = {"Name": ["demo"], "Alias": ("g1",)}

    class FakeFeatureNotFoundError(Exception):
        pass

    class FakeGffutils:
        class exceptions:
            FeatureNotFoundError = FakeFeatureNotFoundError

        class FeatureDB:
            def __init__(self, path, **kwargs):
                assert path == str(gff_db.resolve())
                assert kwargs == {"keep_order": True}

            def __getitem__(self, identifier):
                if identifier == "gene1":
                    return FakeFeature()
                raise FakeFeatureNotFoundError(identifier)

            def region(self, **kwargs):
                assert kwargs == {"region": ("chr1", 1, 10), "featuretype": "gene"}
                return iter([FakeFeature()])

            def close(self):
                return None

    class FakeSbolObject:
        identity = "https://example.org/gene1"
        display_id = "gene1"

    class FakeSbolDocument:
        def read(self, path):
            assert path == str(sbol.resolve())

        def __iter__(self):
            return iter([FakeSbolObject()])

    class FakeSbol3:
        Document = FakeSbolDocument

    class FakeStrand:
        def __init__(self):
            self.domains = [object(), object()]
            self.sequence = "ACGT"

    class FakeDesign:
        def __init__(self):
            self.strands = [FakeStrand()]
            self.helices = {0: object(), 1: object()}
            self.groups = {"default": object()}

        @staticmethod
        def from_scadnano_file(path):
            assert path == str(scadnano_file.resolve())
            return FakeDesign()

    class FakeScadnano:
        Design = FakeDesign

    class FakePolarsResult:
        def limit(self, value):
            assert value == 21
            return self

        def collect(self):
            return self

        def to_dicts(self):
            return [{"chrom": "chr1", "start": 1, "end": 5}]

    class FakePolarsBio:
        @staticmethod
        def overlap(source, other, **kwargs):
            assert source == str(polars_source.resolve())
            assert other == str(polars_other.resolve())
            assert kwargs == {
                "overlap_output": "left",
                "distinct_output": True,
                "output_type": "polars.LazyFrame",
            }
            return FakePolarsResult()

    class FakeGraphicFeature:
        start = 1
        end = 5
        strand = 1
        label = "gene1"

    class FakeGraphicRecord:
        sequence_length = 12

        def __init__(self):
            self.features = [FakeGraphicFeature()]

    class FakeDnaFeaturesTranslator:
        def translate_record(self, path, *, filetype):
            assert path == str(dna_features_source.resolve())
            assert filetype == "genbank"
            return FakeGraphicRecord()

    class FakeDnaFeaturesViewer:
        BiopythonTranslator = FakeDnaFeaturesTranslator

    class FakeGeneBlocksLocation:
        def __init__(self, start, end):
            self.start = start
            self.end = end

    class FakeGeneBlocksBlock:
        operation = "replace"
        s1_location = FakeGeneBlocksLocation(1, 2)
        s2_location = FakeGeneBlocksLocation(1, 2)

    class FakeGeneBlocksDiff:
        def __init__(self):
            self.blocks = [FakeGeneBlocksBlock()]

    class FakeGeneBlocksDiffBlocks:
        @staticmethod
        def from_sequences(reference, target, **kwargs):
            assert (reference, target) == ("AAC", "AGC")
            assert kwargs == {"blast_over": None, "max_complexity": 4_000_000}
            return FakeGeneBlocksDiff()

    class FakeGeneBlocks:
        DiffBlocks = FakeGeneBlocksDiffBlocks

    class FakeSourmashSketch:
        def __init__(self, **kwargs):
            assert kwargs == {
                "n": 0,
                "ksize": 3,
                "scaled": 1_000,
                "track_abundance": False,
            }
            self.hashes = {11: 1, 7: 1}

        def add_sequence(self, sequence, *, force):
            assert sequence in {"AACCGG", "AATCGG"}
            assert force is True

        def similarity(self, other, **kwargs):
            assert kwargs == {"ignore_abundance": True, "downsample": False}
            return 0.5

        def contained_by(self, other):
            return 0.75

    class FakeSourmash:
        MinHash = FakeSourmashSketch

    class FakeGtarsRegion:
        chr = "chr1"
        start = 1
        end = 5
        rest = "peak1"

    class FakeGtarsRegionSet:
        def __init__(self, path):
            assert path == str(polars_source.resolve())

        def __iter__(self):
            return iter([FakeGtarsRegion()])

        def reduce(self):
            return self

    class FakeGtarsModels:
        RegionSet = FakeGtarsRegionSet

    class FakeConstruct:
        id = "construct_1"
        seq = "AACCGGTT"

    class FakeSimulation:
        def __init__(self):
            self.construct_records = [FakeConstruct()]
            self.warnings = []
            self.errors = []

    class FakeRepository:
        def __init__(self):
            self.records = []

        def add_record(self, record):
            self.records.append(record)

    class FakeAssembly:
        def __init__(self, **kwargs):
            assert kwargs["parts"] == ["left", "right"]
            assert kwargs["max_constructs"] == 3
            assert kwargs["expected_constructs"] == "any_number"

        def simulate(self, sequence_repository):
            assert len(sequence_repository.records) == 2
            return FakeSimulation()

    class FakeDnaCauldron:
        SequenceRepository = FakeRepository
        GibsonAssembly = FakeAssembly

    def fake_require(module, extra):
        if module == "pyfastx":
            assert extra == "bio-pyfastx"
            return FakePyfastx
        if module == "pybedtools":
            assert extra == "bio-pybedtools"
            return FakePybedtools
        if module == "gffutils":
            assert extra == "bio-gffutils"
            return FakeGffutils
        if module == "sbol3":
            assert extra == "bio-pysbol3"
            return FakeSbol3
        if module == "scadnano":
            assert extra == "bio-scadnano"
            return FakeScadnano
        if module == "polars_bio":
            assert extra == "bio-polarsbio"
            return FakePolarsBio
        if module == "dna_features_viewer":
            assert extra == "bio-dna-features-viewer"
            return FakeDnaFeaturesViewer
        if module == "geneblocks":
            assert extra == "bio-geneblocks"
            return FakeGeneBlocks
        if module == "sourmash":
            assert extra == "bio-sourmash"
            return FakeSourmash
        if module == "gtars.models":
            assert extra == "bio-gtars"
            return FakeGtarsModels
        if module == "dnacauldron":
            assert extra == "bio-dnacauldron"
            return FakeDnaCauldron
        assert module == "pyBigWig"
        assert extra == "bio-pybigwig"
        return FakePyBigWig

    monkeypatch.setattr(pip_tools, "require", fake_require)
    summary = pip_tools.summarize_pyfastx(fastx, verbose=False)
    assert summary["records"] == 2
    assert summary["total_bases"] == 6
    assert summary["gc_bases"] == 4
    assert summary["format"] == "fasta_or_empty"
    signal = pip_tools.summarize_pybigwig_region(
        bigwig,
        "chr1",
        10,
        20,
        bins=2,
        verbose=False,
    )
    assert signal["values"] == [0.25, None]
    assert signal["provenance"]["citation"] == "10.5281/zenodo.45238"
    intersection = pip_tools.intersect_pybedtools_bed(
        fastx,
        other_bed,
        unique=True,
        write_b=True,
        verbose=False,
    )
    assert intersection["rows"] == [["chr1", "1", "5", "chr1", "2", "6"]]
    assert intersection["provenance"]["bedtools_citation"] == "10.1093/bioinformatics/btq033"
    queried = pip_tools.query_gffutils_database(
        gff_db,
        feature_ids=["gene1"],
        verbose=False,
    )
    assert queried["features"][0]["attributes"]["Name"] == ["demo"]
    region = pip_tools.query_gffutils_database(
        gff_db,
        chromosome="chr1",
        start=1,
        end=10,
        featuretype="gene",
        verbose=False,
    )
    assert region["feature_count"] == 1
    sbol_summary = pip_tools.summarize_sbol3_document(sbol, verbose=False)
    assert sbol_summary["objects"] == 1
    assert sbol_summary["provenance"]["license"] == "MIT"
    scadnano_summary = pip_tools.summarize_scadnano_design(
        scadnano_file,
        verbose=False,
    )
    assert scadnano_summary["helices"] == 2
    assert scadnano_summary["preview"][0]["sequence_length"] == 4
    intervals = pip_tools.overlap_polars_bio_bed(
        polars_source,
        polars_other,
        distinct=True,
        max_rows=20,
        verbose=False,
    )
    assert intervals["row_count"] == 1
    assert intervals["provenance"]["citation"] == "10.1093/bioinformatics/btaf640"
    feature_summary = pip_tools.summarize_dna_features_viewer(
        dna_features_source,
        filetype="genbank",
        verbose=False,
    )
    assert feature_summary["features"] == 1
    assert feature_summary["preview"][0]["label"] == "gene1"
    assert feature_summary["provenance"]["license"] == "MIT"
    diff_summary = pip_tools.diff_geneblocks_sequences(
        "AAC",
        "AGC",
        verbose=False,
    )
    assert diff_summary["block_count"] == 1
    assert diff_summary["blocks"][0]["operation"] == "replace"
    assert diff_summary["provenance"]["license"] == "MIT"
    sourmash_summary = pip_tools.compare_sourmash_sequences(
        "AACCGG",
        "AATCGG",
        ksize=3,
        verbose=False,
    )
    assert sourmash_summary["similarity"] == 0.5
    assert sourmash_summary["reference_hash_preview"] == [7, 11]
    assert sourmash_summary["provenance"]["citation"] == "10.21105/joss.06830"
    gtars_summary = pip_tools.summarize_gtars_regions(
        polars_source,
        verbose=False,
    )
    assert gtars_summary["regions"] == 1
    assert gtars_summary["total_interval_bases"] == 4
    assert gtars_summary["covered_bases"] == 4
    assert gtars_summary["provenance"]["license"] == "BSD-2-Clause"
    assembly = pip_tools.simulate_dnacauldron_assembly(
        {"left": "AACCGG", "right": "CCGGTT"},
        max_constructs=3,
        verbose=False,
    )
    assert assembly["constructs"][0]["length"] == 8
    assert assembly["provenance"]["license"] == "MIT"


def test_scikit_bio_adapters_keep_json_contract(monkeypatch) -> None:
    class FakeArray:
        def tolist(self):
            return [[0.0, 0.5], [0.5, 0.0]]

    class FakeDistance:
        def __init__(self):
            self.data = FakeArray()

    class FakeBeta:
        @staticmethod
        def beta_diversity(metric, counts, ids):
            assert metric == "braycurtis"
            assert counts == [[1.0, 0.0], [0.0, 1.0]]
            assert ids == ["a", "b"]
            return FakeDistance()

    class FakeDistanceModule:
        @staticmethod
        def DistanceMatrix(matrix, ids):
            return matrix, ids

        @staticmethod
        def permanova(distance, grouping, permutations):
            assert grouping == ["case", "case"]
            assert permutations == 0
            return {"test statistic": 0.0, "p-value": 1.0}

    import dnakit.integrations.pip_tools as pip_tools

    def fake_require(module, extra):
        assert extra == "bio-skibio"
        return FakeBeta if module == "skbio.diversity" else FakeDistanceModule

    monkeypatch.setattr(pip_tools, "require", fake_require)
    beta = beta_diversity_scikit_bio([[1, 0], [0, 1]], ["a", "b"], verbose=False)
    assert beta["distance_matrix"] == [[0.0, 0.5], [0.5, 0.0]]
    permanova = pip_tools.permanova_scikit_bio(
        [[0, 0.5], [0.5, 0]], ["a", "b"], ["case", "case"], permutations=0, verbose=False
    )
    assert permanova["result"]["p-value"] == 1.0


def test_unified_capability_schemas_are_shipped() -> None:
    summary = unified_capability_summary(verbose=False)
    assert summary["count"] == 837
    assert summary["modules"]["sequence"] > 0
    capability = describe_unified_capability("CAP-0001", verbose=False)
    assert capability["api"].startswith("dnakit.")
    assert capability["input_schema"]
    assert list_unified_capabilities(module="sequence", verbose=False)


def test_unified_capability_executes_through_agent_registry() -> None:
    result = run_unified_capability(
        "CAP-0001",
        {"input": {"reads": ["ACGT"], "reference": "ACGT"}},
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["api"] == "dnakit.alignment.align_long_reads"
    assert isinstance(result["result"], list)


def test_unified_capability_maps_native_single_input_parameter() -> None:
    result = run_unified_capability(
        "CAP-0187",
        {"input": "@r1\nACGT\n+\nIIII\n"},
        verbose=False,
    )
    assert result["api"] == "dnakit.fasta_fastq.count_reads"
    assert result["result"] == 1
