"""Tests for aligned molecular overlay rendering and PNG export."""

from pathlib import Path

from PIL import Image
import pytest

from conformer_analyzer.models import Atom, Conformer, MoleculeGeometry
from conformer_analyzer.visualization import (
    align_conformers_for_overlay,
    export_overlay_png,
    infer_bonds,
)


def _geometry(
    coordinates: tuple[tuple[float, float, float], ...],
    atomic_numbers: tuple[int, ...] | None = None,
) -> MoleculeGeometry:
    numbers = atomic_numbers or (6,) * len(coordinates)
    return MoleculeGeometry(
        tuple(
            Atom(atomic_number, coordinate)
            for atomic_number, coordinate in zip(numbers, coordinates)
        )
    )


def _conformer(name: str, geometry: MoleculeGeometry) -> Conformer:
    return Conformer(name, Path(f"{name}.log"), geometry)


def test_infer_bonds_uses_conservative_covalent_radius_cutoffs() -> None:
    geometry = _geometry(
        (
            (0.0, 0.0, 0.0),
            (1.09, 0.0, 0.0),
            (4.0, 0.0, 0.0),
        ),
        (6, 1, 8),
    )

    assert infer_bonds(geometry) == ((0, 1),)


def test_infer_bonds_reports_unsupported_elements() -> None:
    geometry = _geometry(((0.0, 0.0, 0.0),), (92,))

    with pytest.raises(ValueError, match=r"atomic number\(s\) 92"):
        infer_bonds(geometry)


def test_overlay_alignment_is_rigid_and_uses_first_conformer_as_reference() -> None:
    reference_coordinates = (
        (0.0, 0.0, 0.0),
        (1.4, 0.0, 0.0),
        (1.4, 1.4, 0.0),
        (0.0, 1.4, 0.2),
    )
    candidate_coordinates = tuple(
        (-y + 4.0, x - 3.0, z + 2.0)
        for x, y, z in reference_coordinates
    )
    reference = _geometry(reference_coordinates)
    candidate = _geometry(candidate_coordinates)

    aligned = align_conformers_for_overlay(
        (_conformer("C1", reference), _conformer("C2", candidate))
    )

    assert aligned[0] is reference
    for reference_atom, aligned_atom in zip(reference.atoms, aligned[1].atoms):
        assert aligned_atom.coordinates == pytest.approx(
            reference_atom.coordinates,
            abs=1e-12,
        )


def test_export_overlay_png_has_white_background_colors_and_legend(
    tmp_path: Path,
) -> None:
    reference = _geometry(
        (
            (0.0, 0.0, 0.0),
            (1.4, 0.0, 0.0),
            (1.4, 1.4, 0.0),
            (0.0, 1.4, 0.0),
        )
    )
    candidate = _geometry(
        (
            (0.0, 0.0, 0.0),
            (1.5, 0.0, 0.0),
            (1.4, 1.5, 0.0),
            (0.0, 1.4, 0.1),
        )
    )
    output_file = tmp_path / "overlay.png"

    result = export_overlay_png(
        (_conformer("C1", reference), _conformer("C2", candidate)),
        output_file,
        labels=("C1", "C2"),
        percentages=(80.0, 20.0),
        width=400,
        height=300,
        dpi=100,
    )

    assert result.output_file == output_file
    assert output_file.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert "inferred from covalent radii" in result.warnings[0]
    with Image.open(output_file) as image:
        assert image.size == (400, 300)
        rgb = image.convert("RGB")
        assert rgb.getpixel((0, 0)) == (255, 255, 255)
        colors = {
            color
            for _, color in (rgb.getcolors(maxcolors=400 * 300) or ())
        }
        assert (239, 32, 32) in colors
        assert (36, 54, 245) in colors


def test_export_reports_visible_atoms_without_inferred_connectivity(
    tmp_path: Path,
) -> None:
    geometry = _geometry(
        (
            (0.0, 0.0, 0.0),
            (1.4, 0.0, 0.0),
            (0.7, 1.2, 0.0),
            (8.0, 8.0, 0.0),
        )
    )

    result = export_overlay_png(
        (_conformer("C1", geometry),),
        tmp_path / "isolated.png",
        width=400,
        height=300,
    )

    assert any("atom position(s): 4" in warning for warning in result.warnings)
