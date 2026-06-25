"""Tests for fitted pairwise RMSD calculations and their numeric kernel."""

from pathlib import Path

import pytest

from conformer_analyzer.models import Atom, MoleculeGeometry
from conformer_analyzer.parsing import parse_gaussian_log
from conformer_analyzer.rmsd import (
    _coordinate_rmsd,
    calculate_rmsd,
    pairwise_rmsd_matrix,
)
from conformer_analyzer.validation import GeometryComparisonError


def _geometry(
    atoms: tuple[tuple[int, tuple[float, float, float]], ...],
) -> MoleculeGeometry:
    return MoleculeGeometry(
        atoms=tuple(
            Atom(atomic_number=atomic_number, coordinates=coordinates)
            for atomic_number, coordinates in atoms
        )
    )


def test_identical_geometries_have_zero_rmsd() -> None:
    geometry = _geometry(
        (
            (8, (0.0, 0.0, 0.0)),
            (1, (1.0, 0.0, 0.0)),
            (1, (-1.0, 0.0, 0.0)),
        )
    )

    assert calculate_rmsd(geometry, geometry) == pytest.approx(0.0)


def test_internal_coordinate_kernel_handles_a_known_translation() -> None:
    reference = ((0.0, 0.0, 0.0), (2.0, -1.0, 4.0))
    translation = (1.0, 2.0, 2.0)
    candidate = tuple(
        (x + translation[0], y + translation[1], z + translation[2])
        for x, y, z in reference
    )

    assert _coordinate_rmsd(reference, candidate) == pytest.approx(3.0)


def test_hydrogen_exclusion_ignores_hydrogen_displacement() -> None:
    reference = _geometry(
        (
            (6, (0.0, 0.0, 0.0)),
            (1, (1.0, 0.0, 0.0)),
        )
    )
    candidate = _geometry(
        (
            (6, (0.0, 0.0, 0.0)),
            (1, (11.0, 0.0, 0.0)),
        )
    )

    assert calculate_rmsd(reference, candidate) > 0.0
    assert calculate_rmsd(
        reference,
        candidate,
        exclude_hydrogens=True,
    ) == pytest.approx(0.0)


def test_geometry_rmsd_rejects_noncomparable_atom_order() -> None:
    reference = _geometry(((6, (0.0, 0.0, 0.0)), (8, (1.0, 0.0, 0.0))))
    candidate = _geometry(((8, (0.0, 0.0, 0.0)), (6, (1.0, 0.0, 0.0))))

    with pytest.raises(GeometryComparisonError, match="atomic-number sequences"):
        calculate_rmsd(reference, candidate)


def test_internal_coordinate_kernel_rejects_different_point_counts() -> None:
    with pytest.raises(ValueError, match="same number of points"):
        _coordinate_rmsd(((0.0, 0.0, 0.0),), ())


def test_internal_coordinate_kernel_rejects_empty_coordinate_sets() -> None:
    with pytest.raises(ValueError, match="at least one coordinate pair"):
        _coordinate_rmsd((), ())


def test_hydrogen_exclusion_rejects_geometry_without_heavy_atoms() -> None:
    hydrogen = _geometry(((1, (0.0, 0.0, 0.0)),))

    with pytest.raises(ValueError, match="at least one non-hydrogen atom"):
        calculate_rmsd(hydrogen, hydrogen, exclude_hydrogens=True)


def test_pairwise_rmsd_matrix_is_square_symmetric() -> None:
    geometries = (
        _geometry(
            (
                (6, (0.0, 0.0, 0.0)),
                (1, (1.0, 0.0, 0.0)),
            )
        ),
        _geometry(
            (
                (6, (0.0, 0.0, 0.0)),
                (1, (2.0, 0.0, 0.0)),
            )
        ),
        _geometry(
            (
                (6, (0.0, 0.0, 0.0)),
                (1, (3.0, 0.0, 0.0)),
            )
        ),
    )

    matrix = pairwise_rmsd_matrix(geometries)

    assert matrix[0][0] == pytest.approx(0.0)
    assert matrix[0][1] == pytest.approx(0.5)
    assert matrix[0][2] == pytest.approx(1.0)
    assert matrix[1][2] == pytest.approx(0.5)
    assert all(matrix[index][index] == 0.0 for index in range(len(matrix)))
    assert all(len(row) == len(matrix) for row in matrix)
    assert all(
        matrix[row][column] == pytest.approx(matrix[column][row])
        for row in range(len(matrix))
        for column in range(len(matrix))
    )


def test_pairwise_rmsd_matrix_always_uses_alignment() -> None:
    base = _geometry(
        (
            (6, (0.0, 0.0, 0.0)),
            (1, (1.0, 0.0, 0.0)),
        )
    )
    rotated = _geometry(
        (
            (6, (0.0, 0.0, 0.0)),
            (1, (0.0, 1.0, 0.0)),
        )
    )

    matrix = pairwise_rmsd_matrix((base, rotated))

    assert matrix == (
        (0.0, pytest.approx(0.0)),
        (pytest.approx(0.0), 0.0),
    )


def test_pairwise_rmsd_matrix_accepts_empty_sequence() -> None:
    assert pairwise_rmsd_matrix(()) == ()


@pytest.mark.parametrize("exclude_hydrogens", (False, True))
def test_pairwise_atom_modes_accept_lists_and_preserve_matrix_invariants(
    exclude_hydrogens: bool,
) -> None:
    geometries = [
        _geometry(
            (
                (6, (0.0, 0.0, 0.0)),
                (8, (1.0, 0.0, 0.0)),
                (7, (0.0, 1.0, 0.0)),
                (1, (0.0, 0.0, 1.0)),
            )
        ),
        _geometry(
            (
                (6, (0.1, 0.0, 0.0)),
                (8, (1.0, 0.1, 0.0)),
                (7, (0.0, 1.0, 0.1)),
                (1, (0.0, 0.0, 1.2)),
            )
        ),
        _geometry(
            (
                (6, (-0.1, 0.0, 0.0)),
                (8, (1.0, -0.1, 0.0)),
                (7, (0.0, 1.0, -0.1)),
                (1, (0.0, 0.0, 0.8)),
            )
        ),
    ]

    matrix = pairwise_rmsd_matrix(
        geometries,
        exclude_hydrogens=exclude_hydrogens,
    )

    assert len(matrix) == len(geometries)
    assert all(len(row) == len(geometries) for row in matrix)
    assert all(matrix[index][index] == 0.0 for index in range(len(matrix)))
    assert all(
        matrix[row][column] == pytest.approx(matrix[column][row])
        for row in range(len(matrix))
        for column in range(len(matrix))
    )


def test_pairwise_matrix_rejects_an_invalid_conformer_set_with_clear_message() -> None:
    reference = _geometry(
        (
            (6, (0.0, 0.0, 0.0)),
            (8, (1.0, 0.0, 0.0)),
        )
    )
    comparable = _geometry(
        (
            (6, (0.0, 1.0, 0.0)),
            (8, (1.0, 1.0, 0.0)),
        )
    )
    invalid = _geometry(
        (
            (8, (0.0, 0.0, 0.0)),
            (6, (1.0, 0.0, 0.0)),
        )
    )

    with pytest.raises(
        GeometryComparisonError,
        match="atomic-number sequences differ at atom 1",
    ):
        pairwise_rmsd_matrix([reference, comparable, invalid])


@pytest.mark.parametrize(
    ("exclude_hydrogens", "expected"),
    (
        (False, 2.1199556774594157),
        (True, 1.5507099839046241),
    ),
)
def test_current_cooh_conformers_have_expected_pairwise_rmsd(
    exclude_hydrogens: bool,
    expected: float,
) -> None:
    conformer_directory = Path(__file__).resolve().parents[1] / "conformers" / "COOH"
    conformers = [
        parse_gaussian_log(conformer_directory / filename).conformers[0]
        for filename in ("conf1.log", "conf2.log")
    ]

    matrix = pairwise_rmsd_matrix(
        conformers,
        exclude_hydrogens=exclude_hydrogens,
    )

    assert matrix == (
        (0.0, pytest.approx(expected)),
        (pytest.approx(expected), 0.0),
    )
