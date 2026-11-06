"""Tabulate imported DNA Cauldron backbone-choice results (pack P35)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.molecular_lab.dnacauldron_p35_functions import _provenance, _record_from_spec


def tabulate_dnacauldron_backbone_choices(
    choices: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return BackboneChoice.list_to_infos_spreadsheet as JSON records.

    Each imported choice needs ``record`` (id, sequence, optional topology/features),
    and may have ``already_on_backbone`` (bool), ``backbone_record``, ``final_record``
    in the same record format, and ``error`` (text). This reports caller-supplied
    insertion results; it does not run backbone selection or infer successful cloning.
    """
    if (
        isinstance(choices, (str, bytes))
        or not isinstance(choices, Sequence)
        or len(choices) > 100000
    ):
        raise ConfigurationError(
            "choices must be a list of at most 100000 backbone-choice records."
        )
    dc = require("dnacauldron.utils", "bio-dnacauldron")
    built = []
    for index, item in enumerate(choices):
        if not isinstance(item, Mapping) or "record" not in item:
            raise ConfigurationError(f"choices[{index}] needs a record mapping.")
        already = item.get("already_on_backbone", False)
        if not isinstance(already, bool):
            raise ConfigurationError("already_on_backbone must be a boolean.")
        error = item.get("error")
        if error is not None and not isinstance(error, str):
            raise ConfigurationError("error must be text or null.")
        record = _record_from_spec(item["record"])
        backbone = (
            None
            if item.get("backbone_record") is None
            else _record_from_spec(item["backbone_record"])
        )
        final = (
            None if item.get("final_record") is None else _record_from_spec(item["final_record"])
        )
        built.append(
            dc.BackboneChoice(
                record,
                already_on_backbone=already,
                backbone_record=backbone,
                final_record=final,
                error=None if error is None else ValueError(error),
            )
        )
    table = dc.BackboneChoice.list_to_infos_spreadsheet(built)
    return {
        "status": "completed",
        "count": len(table),
        "columns": list(table.columns),
        "choices": table.to_dict(orient="records"),
        "provenance": _provenance(method="BackboneChoice.list_to_infos_spreadsheet"),
    }


__all__ = ["tabulate_dnacauldron_backbone_choices"]
_register_public_verbose(globals())
