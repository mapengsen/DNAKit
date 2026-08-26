"""Tests for the single internal coordinate convention."""

import pytest

from dnakit.core import (
    CompoundLocation,
    CoordinateSystem,
    ExternalInterval,
    Interval,
    Strand,
    UnresolvedLocation,
    export_location,
    import_location,
    reverse_strand_location,
)
from dnakit.exceptions import CoordinateError


@pytest.mark.parametrize(
    ("external", "expected"),
    [
        (ExternalInterval(0, 4, CoordinateSystem.ZERO_BASED_HALF_OPEN), Interval(0, 4)),
        (ExternalInterval(0, 3, CoordinateSystem.ZERO_BASED_CLOSED), Interval(0, 4)),
        (ExternalInterval(1, 4, CoordinateSystem.ONE_BASED_CLOSED), Interval(0, 4)),
        (ExternalInterval(1, 5, CoordinateSystem.ONE_BASED_HALF_OPEN), Interval(0, 4)),
    ],
)
def test_external_systems_import_to_zero_based_half_open(
    external: ExternalInterval,
    expected: Interval,
) -> None:
    assert import_location(external) == expected


@pytest.mark.parametrize(
    ("system", "expected"),
    [
        (CoordinateSystem.ZERO_BASED_HALF_OPEN, (0, 4)),
        (CoordinateSystem.ZERO_BASED_CLOSED, (0, 3)),
        (CoordinateSystem.ONE_BASED_CLOSED, (1, 4)),
        (CoordinateSystem.ONE_BASED_HALF_OPEN, (1, 5)),
    ],
)
def test_internal_interval_exports_without_hidden_coordinate_state(
    system: CoordinateSystem,
    expected: tuple[int, int],
) -> None:
    (external,) = export_location(Interval(0, 4), target_system=system)

    assert (external.start, external.end) == expected
    assert external.system is system
    assert external.strand is Strand.UNKNOWN


def test_circular_origin_wrap_becomes_ordered_compound_location() -> None:
    wrapped = import_location(
        ExternalInterval(9, 2, CoordinateSystem.ONE_BASED_CLOSED, Strand.FORWARD),
        sequence_length=10,
    )

    assert wrapped == CompoundLocation((Interval(8, 10), Interval(0, 2)))
    assert export_location(
        wrapped,
        target_system=CoordinateSystem.ONE_BASED_CLOSED,
        sequence_length=10,
    ) == (
        ExternalInterval(9, 10, CoordinateSystem.ONE_BASED_CLOSED),
        ExternalInterval(1, 2, CoordinateSystem.ONE_BASED_CLOSED),
    )


def test_adjacent_closed_endpoints_across_origin_represent_full_circle() -> None:
    wrapped = import_location(
        ExternalInterval(5, 4, CoordinateSystem.ONE_BASED_CLOSED),
        sequence_length=10,
    )

    assert wrapped == CompoundLocation((Interval(4, 10), Interval(0, 4)))
    assert len(wrapped) == 10


def test_wrapped_interval_requires_sequence_length() -> None:
    with pytest.raises(CoordinateError, match="requires sequence_length"):
        import_location(ExternalInterval(9, 2, CoordinateSystem.ONE_BASED_CLOSED))


def test_reverse_strand_location_uses_sequence_length() -> None:
    assert reverse_strand_location(Interval(2, 5), sequence_length=10) == Interval(5, 8)

    compound = CompoundLocation((Interval(2, 4), Interval(8, 10)))
    assert reverse_strand_location(compound, sequence_length=12) == CompoundLocation(
        (Interval(2, 4), Interval(8, 10))
    )


def test_unresolved_location_keeps_reason_and_cannot_be_exported() -> None:
    location = UnresolvedLocation("unknown gap", anchors=[Interval(2, 3)])
    reversed_location = reverse_strand_location(location, sequence_length=10)

    assert reversed_location == UnresolvedLocation("unknown gap", [Interval(7, 8)])
    with pytest.raises(CoordinateError, match="cannot be exported"):
        export_location(location, target_system=CoordinateSystem.ZERO_BASED_HALF_OPEN)


def test_coordinate_invariants_reject_negative_reversed_and_out_of_bounds_values() -> None:
    with pytest.raises(CoordinateError):
        Interval(-1, 2)
    with pytest.raises(CoordinateError):
        Interval(3, 2)
    with pytest.raises(CoordinateError, match="exceeds"):
        import_location(
            ExternalInterval(0, 11, CoordinateSystem.ZERO_BASED_HALF_OPEN),
            sequence_length=10,
        )
    with pytest.raises(CoordinateError, match="closed"):
        export_location(Interval(2, 2), target_system=CoordinateSystem.ONE_BASED_CLOSED)
