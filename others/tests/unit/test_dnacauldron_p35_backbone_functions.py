"""Backbone result tables compared with direct DNA Cauldron calls."""

import pytest

dc = pytest.importorskip("dnacauldron")
from dnacauldron.utils import BackboneChoice

from dnakit.molecular_lab.dnacauldron_p35_backbone_functions import (
    tabulate_dnacauldron_backbone_choices,
)
from dnakit.molecular_lab.dnacauldron_p35_functions import _record_from_spec


@pytest.mark.parametrize(
    "choices",
    [
        [],
        [{"record": {"id": "x", "sequence": "ACGT"}}],
        [
            {
                "record": {"id": "x", "sequence": "A"},
                "already_on_backbone": True,
                "backbone_record": {"id": "backbone", "sequence": "ACGT"},
                "final_record": {"id": "assembled", "sequence": "ACGTA"},
            },
            {"record": {"id": "y", "sequence": "CC"}, "error": "No compatible backbone"},
        ],
    ],
)
def test_backbone_summary_matches_upstream(choices):
    direct = []
    for item in choices:
        direct.append(
            BackboneChoice(
                _record_from_spec(item["record"]),
                already_on_backbone=item.get("already_on_backbone", False),
                backbone_record=None
                if "backbone_record" not in item
                else _record_from_spec(item["backbone_record"]),
                final_record=None
                if "final_record" not in item
                else _record_from_spec(item["final_record"]),
                error=None if "error" not in item else ValueError(item["error"]),
            )
        )
    expected = BackboneChoice.list_to_infos_spreadsheet(direct)
    result = tabulate_dnacauldron_backbone_choices(choices, verbose=False)
    assert result["choices"] == expected.to_dict(orient="records")
    assert result["columns"] == list(expected.columns)
