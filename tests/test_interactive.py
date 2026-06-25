"""Tests for the local py3Dmol conformer viewer export."""

from pathlib import Path

from conformer_analyzer.interactive import export_interactive_overlay_html
from conformer_analyzer.models import Atom, Conformer, MoleculeGeometry


def _conformer(
    name: str,
    coordinates: tuple[tuple[float, float, float], ...],
) -> Conformer:
    geometry = MoleculeGeometry(
        tuple(Atom(6, coordinate) for coordinate in coordinates)
    )
    return Conformer(name, Path(f"{name}.log"), geometry)


def test_interactive_export_embeds_aligned_models_and_controls(
    tmp_path: Path,
) -> None:
    reference = _conformer(
        "conf1",
        (
            (0.0, 0.0, 0.0),
            (1.4, 0.0, 0.0),
            (1.4, 1.4, 0.0),
            (0.0, 1.4, 0.2),
        ),
    )
    candidate = _conformer(
        "conf2",
        (
            (4.0, -3.0, 2.0),
            (4.0, -1.6, 2.0),
            (2.6, -1.6, 2.0),
            (2.6, -3.0, 2.2),
        ),
    )
    output_file = tmp_path / "viewer.html"

    result = export_interactive_overlay_html(
        (reference, candidate),
        output_file,
        labels=("C1", "C2"),
    )

    html = output_file.read_text(encoding="utf-8")
    assert result.output_file == output_file
    assert len(result.bonds) == 4
    assert "3dmol@2.5.5" in html
    assert html.count("addModel") == 2
    assert "toggleModel(0" in html
    assert "toggleModel(1" in html
    assert html.count('class="color-picker" type="color"') == 2
    assert 'value="#ef2020"' in html
    assert 'value="#2436f5"' in html
    assert "setConformerColor(0" in html
    assert "setConformerColor(1" in html
    assert "function applyConformerStyle" in html
    assert "toggleHydrogens" in html
    assert "resetView" in html
    assert "saveSnapshot" in html
    assert "Drag to rotate" in html
    assert "Connectivity was inferred" in result.warnings[0]
