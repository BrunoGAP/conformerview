"""Tests for atom-by-atom geometry comparability validation."""

from pathlib import Path

import pytest

from conformer_analyzer.elements import ATOMIC_NUMBER_TO_SYMBOL
from conformer_analyzer.models import Atom, MoleculeGeometry
from conformer_analyzer.parsing import parse_gaussian_log, parse_xyz_text
from conformer_analyzer.validation import (
    GeometryComparisonError,
    atomic_number_sequence,
    heavy_atom_geometry,
    validate_comparable_geometries,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFORMERS = PROJECT_ROOT / "conformers" / "COOH"


def _geometry(*atomic_numbers: int) -> MoleculeGeometry:
    return MoleculeGeometry(
        atoms=tuple(
            Atom(
                atomic_number=atomic_number,
                coordinates=(float(index), 0.0, 0.0),
            )
            for index, atomic_number in enumerate(atomic_numbers)
        )
    )


def test_matching_atomic_number_sequences_are_comparable() -> None:
    reference = _geometry(8, 6, 1, 1)
    candidate = MoleculeGeometry(
        atoms=tuple(
            Atom(atom.atomic_number, (atom.x + 5.0, atom.y, atom.z))
            for atom in reference.atoms
        )
    )

    assert validate_comparable_geometries(reference, candidate) is None
    assert atomic_number_sequence(reference) == (8, 6, 1, 1)


def test_different_atom_counts_fail_with_clear_message() -> None:
    with pytest.raises(GeometryComparisonError) as error:
        validate_comparable_geometries(_geometry(8, 1, 1), _geometry(8, 1))

    assert "atom counts differ" in str(error.value)
    assert "reference=3, candidate=2" in str(error.value)


def test_different_atom_order_fails_without_reordering() -> None:
    reference = _geometry(6, 8, 1)
    candidate = _geometry(8, 6, 1)

    with pytest.raises(GeometryComparisonError) as error:
        validate_comparable_geometries(reference, candidate)

    assert "atomic-number sequences differ at atom 1" in str(error.value)
    assert "reference=6, candidate=8" in str(error.value)
    assert "Atom reordering is not performed" in str(error.value)
    assert atomic_number_sequence(candidate) == (8, 6, 1)


def test_heavy_atom_filter_preserves_relative_order_and_objects() -> None:
    geometry = _geometry(1, 6, 1, 8, 7, 1, 16)

    filtered = heavy_atom_geometry(geometry)

    assert atomic_number_sequence(filtered) == (6, 8, 7, 16)
    assert filtered.atoms == (
        geometry.atoms[1],
        geometry.atoms[3],
        geometry.atoms[4],
        geometry.atoms[6],
    )
    assert atomic_number_sequence(geometry) == (1, 6, 1, 8, 7, 1, 16)


def test_extracted_example_geometries_are_comparable() -> None:
    conf1 = parse_gaussian_log(CONFORMERS / "conf1.log").conformers[0]
    conf2 = parse_gaussian_log(CONFORMERS / "conf2.log").conformers[0]

    validate_comparable_geometries(conf1.geometry, conf2.geometry)
    conf1_heavy = heavy_atom_geometry(conf1.geometry)
    conf2_heavy = heavy_atom_geometry(conf2.geometry)

    assert len(conf1.geometry.atoms) == len(conf2.geometry.atoms) == 58
    assert len(conf1_heavy.atoms) == len(conf2_heavy.atoms) == 37
    assert atomic_number_sequence(conf1_heavy) == atomic_number_sequence(conf2_heavy)


def test_xyz_geometry_can_validate_against_gaussian_geometry() -> None:
    conf2 = parse_gaussian_log(CONFORMERS / "conf2.log").conformers[0]
    xyz_rows = "\n".join(
        f"{ATOMIC_NUMBER_TO_SYMBOL[atom.atomic_number]} "
        f"{atom.x:.8f} {atom.y:.8f} {atom.z:.8f}"
        for atom in conf2.geometry.atoms
    )
    xyz = parse_xyz_text(
        f"{len(conf2.geometry.atoms)}\nconverted from Gaussian\n{xyz_rows}\n",
        "conf2.xyz",
    ).conformers[0]

    validate_comparable_geometries(conf2.geometry, xyz.geometry)
    assert atomic_number_sequence(xyz.geometry) == atomic_number_sequence(conf2.geometry)
