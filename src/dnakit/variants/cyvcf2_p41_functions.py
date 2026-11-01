"""Native cyvcf2 IBD, QC-site selection and textual export (P41)."""
from __future__ import annotations
from itertools import islice
from typing import Any
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.exceptions import ConfigurationError

def cyvcf2_ibd_bins(path: str, *, max_variants: int = 100000, verbose: bool = True) -> dict[str, Any]:
    """Return the 16 native IBD run-length bins for every sample pair (gts012=True)."""
    source = h.existing(path, 'path')
    limit = h.integer(max_variants, 'max_variants', minimum=1, maximum=1000000)
    reader = require('cyvcf2', 'bio-cyvcf2').VCF(str(source), gts012=True)
    try:
        if len(reader.samples) > 100:
            raise ConfigurationError('At most 100 samples are supported for bounded pairwise output.')
        values = reader.ibd(nmax=limit)
    finally:
        reader.close()
    return {'status': 'completed', 'backend': 'cyvcf2', 'max_variants': limit, 'pairs': [{'sample_a': a, 'sample_b': b, 'bins': x.tolist()} for (a,b), x in values.items()]}

def cyvcf2_informative_sites(path: str, *, sites: list[list[Any]] | None = None, offset: int = 0, each: int = 1, call_rate: float = 0.8, max_records: int = 1000, verbose: bool = True) -> dict[str, Any]:
    """Select informative sites via native gen_variants, including its native no-sites QC rules.

    Explicit sites use [chromosome, 1-based position, optional REF, ALT] rows and require an index.
    Without sites upstream applies its fixed depth, frequency and 5kb spacing filters.
    """
    source = h.existing(path, 'path')
    limit = h.integer(max_records, 'max_records', minimum=1, maximum=100000)
    if not 0 <= call_rate <= 1:
        raise ConfigurationError('call_rate must lie between zero and one.')
    reader = require('cyvcf2', 'bio-cyvcf2').VCF(str(source), gts012=True)
    try:
        records = list(islice(reader.gen_variants(sites, offset=h.integer(offset, 'offset', minimum=0), each=h.integer(each, 'each', minimum=1), call_rate=call_rate), limit + 1))
        rows = [{'index': int(i), 'chrom': v.CHROM, 'position': v.POS, 'ref': v.REF, 'alt': list(v.ALT), 'call_rate': float(v.call_rate)} for i,v in records[:limit]]
    finally:
        reader.close()
    return {'status': 'completed', 'backend': 'cyvcf2', 'records': rows, 'truncated': len(records) > limit}

def cyvcf2_export_text(path: str, output: str, *, max_records: int = 1000000, verbose: bool = True) -> dict[str, Any]:
    """Write the actual VCF header and record text from a VCF/BCF input."""
    source, target = h.existing(path, 'path'), h.new_file(output, 'output')
    limit = h.integer(max_records, 'max_records', minimum=1, maximum=10000000)
    reader = require('cyvcf2', 'bio-cyvcf2').VCF(str(source))
    count = 0
    try:
        with target.open('x') as handle:
            handle.write(reader.raw_header)
            for record in reader:
                if count == limit:
                    raise ConfigurationError('Input exceeds max_records; no complete export was written.')
                handle.write(str(record))
                count += 1
    except Exception:
        target.unlink(missing_ok=True)
        raise
    finally:
        reader.close()
    return {'status': 'completed', 'backend': 'cyvcf2', 'record_count': count, 'outputs': [str(target)]}

def vcf_call_is_filtered(genotype: str, *, filters: list[str] | None = None, verbose: bool = True) -> dict[str, Any]:
    """Return vcfpy.Call.is_filtered for an explicit genotype call and FT list."""
    data = {'GT': genotype}
    if filters is not None:
        data['FT'] = filters
    value = require('vcfpy', 'bio-vcfpy').Call('sample', data).is_filtered()
    return {'status': 'completed', 'backend': 'vcfpy', 'is_filtered': bool(value)}

__all__ = ['cyvcf2_ibd_bins', 'cyvcf2_informative_sites', 'cyvcf2_export_text', 'vcf_call_is_filtered']
_register_public_verbose(globals(), tuple(__all__))
