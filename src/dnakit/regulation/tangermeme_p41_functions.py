"""JSON interfaces for tangermeme numeric utilities (P41)."""
from __future__ import annotations
from typing import Any
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.exceptions import ConfigurationError

def _tensor(value: Any) -> Any:
    torch = require('torch', 'bio-tangermeme')
    try:
        result = torch.tensor(value, dtype=torch.float32)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError('Expected rectangular numeric input.') from exc
    if result.numel() > 1000000 or not torch.isfinite(result).all():
        raise ConfigurationError('Input must contain at most one million finite values.')
    return result

def tangermeme_hypothetical_attributions(multipliers: list[Any], inputs: list[Any], references: list[Any], *, verbose: bool = True) -> dict[str, Any]:
    """Compute hypothetical DeepLIFT attributions; arrays have native batch/channel/length shape."""
    tensors = [_tensor(x) for x in (multipliers, inputs, references)]
    if any(x.ndim != 3 for x in tensors) or len({tuple(x.shape) for x in tensors}) != 1:
        raise ConfigurationError('Require equal-shaped three-dimensional arrays.')
    values = require('tangermeme.deep_lift_shap', 'bio-tangermeme').hypothetical_attributions(*[(x,) for x in tensors])
    return {'status': 'completed', 'backend': 'tangermeme', 'attributions': values[0].detach().cpu().tolist()}

def tangermeme_unchunk(chunks: list[Any], lengths: list[int], *, overlap: int = 0, verbose: bool = True) -> dict[str, Any]:
    """Restore native sequence outputs from chunks; upstream overlap trimming is preserved."""
    data = _tensor(chunks)
    amount = h.integer(overlap, 'overlap', minimum=0)
    if data.ndim != 3 or not lengths or amount >= data.shape[-1]:
        raise ConfigurationError('Require 3D chunks, lengths and overlap smaller than chunk length.')
    size = data.shape[-1]
    values = [h.integer(x, 'length', minimum=size) for x in lengths]
    if sum((x-size)//(size-amount)+1 for x in values) != data.shape[0]:
        raise ConfigurationError('Lengths do not account for exactly the supplied chunks.')
    out = require('tangermeme.utils', 'bio-tangermeme').unchunk(data, values, overlap=amount)
    return {'status': 'completed', 'backend': 'tangermeme', 'sequences': [x.tolist() for x in out]}

def tangermeme_extract_signal(signal: list[Any], loci: list[list[int]], *, verbose: bool = True) -> dict[str, Any]:
    """Sum per-base signal for rows [example_index, start, end], half-open coordinates."""
    data = _tensor(signal)
    if data.ndim != 3:
        raise ConfigurationError('Signal must have example/channel/position dimensions.')
    for row in loci:
        if len(row) != 3 or any(type(x) is not int for x in row) or not (0 <= row[0] < data.shape[0] and 0 <= row[1] < row[2] <= data.shape[2]):
            raise ConfigurationError('Invalid locus coordinates.')
    pd = require('pandas', 'bio-tangermeme')
    out = require('tangermeme.utils', 'bio-tangermeme').extract_signal(pd.DataFrame(loci), data, verbose=verbose)
    return {'status': 'completed', 'backend': 'tangermeme', 'signal_sums': out.tolist()}

def tangermeme_position_entropy(probabilities: list[Any], *, eps: float = 1e-9, verbose: bool = True) -> dict[str, Any]:
    """Compute native per-position entropy on batch/alphabet/position probability tensors."""
    data = _tensor(probabilities)
    if data.ndim != 3 or (data < 0).any() or (data > 1).any() or not 0 < eps < 1:
        raise ConfigurationError('Use 3D probabilities between zero and one and 0 < eps < 1.')
    out = require('tangermeme.utils', 'bio-tangermeme').entropy(data, eps=eps)
    return {'status': 'completed', 'backend': 'tangermeme', 'entropy': out.tolist()}

__all__ = ['tangermeme_hypothetical_attributions', 'tangermeme_unchunk', 'tangermeme_extract_signal', 'tangermeme_position_entropy']
_register_public_verbose(globals(), tuple(__all__))
