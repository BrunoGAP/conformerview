"""Tests for the command-line conformer analysis workflow."""

from pathlib import Path

import pytest

from conformer_analyzer.elements import ATOMIC_NUMBER_TO_SYMBOL
from conformer_analyzer.interface import analyze_log_files
from conformer_analyzer.parsing import parse_gaussian_log


PROJECT_ROOT = Path(__file__).resolve().parents[1]
COOH_LOGS = (
    PROJECT_ROOT / "conformers" / "COOH" / "conf1.log",
    PROJECT_ROOT / "conformers" / "COOH" / "conf2.log",
)


def test_analysis_workflow_writes_both_matrices_and_interactive_viewer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    opened_uris: list[str] = []
    monkeypatch.setattr(
        "conformer_analyzer.interface.webbrowser.open_new_tab",
        lambda uri: opened_uris.append(uri) or True,
    )

    result = analyze_log_files(
        COOH_LOGS,
        tmp_path,
        viewer_fit_atoms="heavy",
        show_hydrogens=True,
        open_viewer=True,
    )

    assert [conformer.name for conformer in result.conformers] == ["conf1", "conf2"]
    assert result.all_atom_matrix[0][1] == pytest.approx(2.1199556774594157)
    assert result.heavy_atom_matrix[0][1] == pytest.approx(1.5507099839046241)
    assert result.matrix_file == tmp_path / "rmsd_matrix.txt"
    assert result.viewer_file == tmp_path / "interactive_overlay.html"
    assert result.viewer_opened is True
    assert opened_uris == [result.viewer_file.resolve().as_uri()]
    assert any(
        "TD optimization job 1" in warning
        for warning in result.warnings
    )

    report = result.matrix_file.read_text(encoding="utf-8")
    assert "All atoms (fit and RMSD use all atoms)" in report
    assert "2.1199556775" in report
    assert "Heavy atoms (fit and RMSD exclude hydrogen)" in report
    assert "1.5507099839" in report
    assert "TD optimization job 1" in report
    assert "Connectivity was inferred" in report

    viewer_html = result.viewer_file.read_text(encoding="utf-8")
    assert "Aligned conformer overlay" in viewer_html
    assert 'id="hydrogens" type="checkbox" checked' in viewer_html


def test_analysis_requires_at_least_two_log_files(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="at least two Gaussian .log files"):
        analyze_log_files((COOH_LOGS[0],), tmp_path, open_viewer=False)


def test_analysis_rejects_unknown_viewer_fit_mode(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="must be 'all' or 'heavy'"):
        analyze_log_files(
            COOH_LOGS,
            tmp_path,
            viewer_fit_atoms="carbon",
            open_viewer=False,
        )


def test_analysis_accepts_mixed_gaussian_and_xyz_inputs(tmp_path: Path) -> None:
    conf2 = parse_gaussian_log(COOH_LOGS[1]).conformers[0]
    xyz_rows = "\n".join(
        f"{ATOMIC_NUMBER_TO_SYMBOL[atom.atomic_number]} "
        f"{atom.x:.8f} {atom.y:.8f} {atom.z:.8f}"
        for atom in conf2.geometry.atoms
    )
    xyz_file = tmp_path / "conf2_copy.xyz"
    xyz_file.write_text(
        f"{len(conf2.geometry.atoms)}\nconverted from Gaussian\n{xyz_rows}\n",
        encoding="utf-8",
    )

    result = analyze_log_files(
        (COOH_LOGS[1], xyz_file),
        tmp_path / "outputs",
        open_viewer=False,
    )

    assert [conformer.name for conformer in result.conformers] == [
        "conf2",
        "conf2_copy",
    ]
    assert result.all_atom_matrix[0][1] == pytest.approx(0.0, abs=1e-10)
    assert result.heavy_atom_matrix[0][1] == pytest.approx(0.0, abs=1e-10)
