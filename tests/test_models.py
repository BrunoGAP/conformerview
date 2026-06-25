"""Tests for the molecular data models."""

from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from conformer_analyzer.models import (
    Atom,
    Conformer,
    MoleculeGeometry,
    ParseResult,
    ParseWarning,
)


def test_small_molecule_preserves_atom_order() -> None:
    oxygen = Atom(atomic_number=8, symbol="O", coordinates=(0.0, 0.0, 0.0))
    first_hydrogen = Atom(
        atomic_number=1,
        symbol="H",
        coordinates=(0.7586, 0.0, 0.5043),
    )
    second_hydrogen = Atom(
        atomic_number=1,
        symbol="H",
        coordinates=(-0.7586, 0.0, 0.5043),
    )

    geometry = MoleculeGeometry(
        atoms=(oxygen, first_hydrogen, second_hydrogen),
    )

    assert geometry.atoms == (oxygen, first_hydrogen, second_hydrogen)
    assert tuple(atom.atomic_number for atom in geometry.atoms) == (8, 1, 1)
    assert first_hydrogen.x == pytest.approx(0.7586)
    assert first_hydrogen.y == pytest.approx(0.0)
    assert first_hydrogen.z == pytest.approx(0.5043)


def test_conformer_stores_optional_calculation_metadata() -> None:
    geometry = MoleculeGeometry(
        atoms=(Atom(atomic_number=6, symbol="C", coordinates=(0.0, 0.0, 0.0)),),
    )

    conformer = Conformer(
        name="conf1",
        source_file=Path("conformers/COOH/conf1.log"),
        geometry=geometry,
        charge=0,
        multiplicity=1,
        energy_hartree=-40.123456,
    )

    assert conformer.geometry is geometry
    assert conformer.charge == 0
    assert conformer.multiplicity == 1
    assert conformer.energy_hartree == pytest.approx(-40.123456)


def test_parse_result_keeps_structured_warnings() -> None:
    warning = ParseWarning(
        message="Multiple eligible optimization jobs found.",
        job_index=2,
        line_number=145,
    )
    result = ParseResult(
        source_file=Path("conformers/COOH/conf2.log"),
        conformers=(),
        warnings=(warning,),
    )

    assert result.warnings == (warning,)
    assert result.warnings[0].job_index == 2


def test_models_are_immutable() -> None:
    atom = Atom(atomic_number=1, coordinates=(0.0, 0.0, 0.0))

    with pytest.raises(FrozenInstanceError):
        atom.atomic_number = 8  # type: ignore[misc]


@pytest.mark.parametrize(
    ("atomic_number", "coordinates", "message"),
    (
        (0, (0.0, 0.0, 0.0), "atomic_number"),
        (1, (0.0, 0.0), "coordinates"),
    ),
)
def test_atom_rejects_invalid_core_values(
    atomic_number: int,
    coordinates: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        Atom(
            atomic_number=atomic_number,
            coordinates=coordinates,  # type: ignore[arg-type]
        )


def test_conformer_rejects_invalid_multiplicity() -> None:
    geometry = MoleculeGeometry(atoms=())

    with pytest.raises(ValueError, match="multiplicity"):
        Conformer(
            name="invalid",
            source_file=Path("invalid.log"),
            geometry=geometry,
            multiplicity=0,
        )
