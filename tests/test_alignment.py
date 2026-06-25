"""Tests for optimal rigid-body molecular alignment."""

from math import dist

import pytest

from conformer_analyzer.alignment import align_coordinates, align_geometry, centroid
from conformer_analyzer.models import Atom, MoleculeGeometry
from conformer_analyzer.rmsd import calculate_rmsd


def _geometry(
    coordinates: tuple[tuple[float, float, float], ...],
    atomic_numbers: tuple[int, ...] = (6, 7, 8, 16),
) -> MoleculeGeometry:
    return MoleculeGeometry(
        atoms=tuple(
            Atom(atomic_number, coordinate)
            for atomic_number, coordinate in zip(atomic_numbers, coordinates)
        )
    )


def _rotate_and_translate(
    coordinates: tuple[tuple[float, float, float], ...],
) -> tuple[tuple[float, float, float], ...]:
    return tuple((-y + 4.0, x - 3.0, z + 2.5) for x, y, z in coordinates)


def test_centroid_is_arithmetic_coordinate_mean() -> None:
    coordinates = ((0.0, 1.0, 2.0), (2.0, 3.0, 4.0), (4.0, 5.0, 6.0))

    assert centroid(coordinates) == pytest.approx((2.0, 3.0, 4.0))


def test_rotated_translated_copy_has_near_zero_aligned_rmsd() -> None:
    reference_coordinates = (
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 2.5),
    )
    reference = _geometry(reference_coordinates)
    candidate = _geometry(_rotate_and_translate(reference_coordinates))

    assert calculate_rmsd(reference, candidate) == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_alignment_preserves_internal_distances_and_atom_order() -> None:
    reference_coordinates = (
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (0.2, 1.8, 0.0),
        (0.1, 0.3, 2.2),
    )
    reference = _geometry(reference_coordinates)
    candidate = _geometry(_rotate_and_translate(reference_coordinates))

    aligned = align_geometry(reference, candidate)

    assert tuple(atom.atomic_number for atom in aligned.atoms) == (6, 7, 8, 16)
    for first_index in range(len(candidate.atoms)):
        for second_index in range(first_index + 1, len(candidate.atoms)):
            before = dist(
                candidate.atoms[first_index].coordinates,
                candidate.atoms[second_index].coordinates,
            )
            after = dist(
                aligned.atoms[first_index].coordinates,
                aligned.atoms[second_index].coordinates,
            )
            assert after == pytest.approx(before, abs=1e-12)


def test_heavy_atom_alignment_ignores_hydrogen_displacement() -> None:
    reference_coordinates = (
        (0.0, 0.0, 0.0),
        (1.5, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.4, 0.5, 1.0),
    )
    candidate_coordinates = list(_rotate_and_translate(reference_coordinates))
    hydrogen_x, hydrogen_y, hydrogen_z = candidate_coordinates[-1]
    candidate_coordinates[-1] = (
        hydrogen_x + 3.0,
        hydrogen_y - 2.0,
        hydrogen_z + 1.0,
    )
    reference = _geometry(reference_coordinates, (6, 7, 8, 1))
    candidate = _geometry(tuple(candidate_coordinates), (6, 7, 8, 1))

    assert calculate_rmsd(
        reference,
        candidate,
        exclude_hydrogens=True,
    ) == pytest.approx(0.0, abs=1e-12)
    assert calculate_rmsd(reference, candidate) > 0.1


def test_align_coordinates_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="at least one coordinate pair"):
        align_coordinates((), ())


def test_reflection_is_not_treated_as_rotation() -> None:
    reference_coordinates = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    reflected_coordinates = tuple(
        (-x, y, z) for x, y, z in reference_coordinates
    )
    reference = _geometry(reference_coordinates)
    reflected = _geometry(reflected_coordinates)

    assert calculate_rmsd(reference, reflected) > 0.1
