"""Command-line workflow for conformer comparison."""

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import webbrowser

from conformer_analyzer.interactive import export_interactive_overlay_html
from conformer_analyzer.models import Conformer
from conformer_analyzer.parsing import parse_conformer_file
from conformer_analyzer.rmsd import pairwise_rmsd_matrix
from conformer_analyzer.validation import validate_comparable_geometries


RmsdMatrix = tuple[tuple[float, ...], ...]


@dataclass(frozen=True, slots=True)
class AnalysisRun:
    """Files, values, and warnings produced by one CLI analysis."""

    conformers: tuple[Conformer, ...]
    all_atom_matrix: RmsdMatrix
    heavy_atom_matrix: RmsdMatrix
    matrix_file: Path
    viewer_file: Path
    warnings: tuple[str, ...]
    viewer_opened: bool


SUPPORTED_INPUT_SUFFIXES = {".log", ".xyz"}


def analyze_log_files(
    input_files: Sequence[str | Path],
    output_directory: str | Path = "outputs",
    *,
    viewer_fit_atoms: str = "all",
    show_hydrogens: bool = False,
    open_viewer: bool = True,
) -> AnalysisRun:
    """Run the complete validated RMSD and interactive-viewer workflow."""

    paths = tuple(Path(path) for path in input_files)
    if len(paths) < 2:
        raise ValueError(
            "Analysis requires at least two Gaussian .log files or XYZ .xyz files."
        )
    if viewer_fit_atoms not in {"all", "heavy"}:
        raise ValueError("viewer_fit_atoms must be 'all' or 'heavy'.")

    conformers: list[Conformer] = []
    warnings: list[str] = []
    for path in paths:
        if path.suffix.lower() not in SUPPORTED_INPUT_SUFFIXES:
            raise ValueError(f"Input file must use the .log or .xyz extension: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"Input file was not found: {path}")
        parse_result = parse_conformer_file(path)
        warnings.extend(
            f"{path.name}: {warning.message}"
            for warning in parse_result.warnings
        )
        if not parse_result.conformers:
            raise ValueError(f"No usable conformer was extracted from {path}.")
        conformers.extend(parse_result.conformers)

    conformer_tuple = tuple(conformers)
    if len(conformer_tuple) < 2:
        raise ValueError("At least two usable conformers are required.")
    reference = conformer_tuple[0].geometry
    for conformer in conformer_tuple[1:]:
        validate_comparable_geometries(reference, conformer.geometry)

    all_atom_matrix = pairwise_rmsd_matrix(conformer_tuple)
    heavy_atom_matrix = pairwise_rmsd_matrix(
        conformer_tuple,
        exclude_hydrogens=True,
    )
    labels = tuple(f"C{index}" for index in range(1, len(conformer_tuple) + 1))

    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    matrix_file = output_path / "rmsd_matrix.txt"
    viewer_file = output_path / "interactive_overlay.html"
    interactive_result = export_interactive_overlay_html(
        conformer_tuple,
        viewer_file,
        labels=labels,
        exclude_hydrogens_from_fit=viewer_fit_atoms == "heavy",
        show_hydrogens=show_hydrogens,
    )
    warnings.extend(interactive_result.warnings)

    report = format_rmsd_report(
        conformer_tuple,
        labels,
        all_atom_matrix,
        heavy_atom_matrix,
        warnings=warnings,
    )
    matrix_file.write_text(report, encoding="utf-8")

    viewer_opened = False
    if open_viewer:
        viewer_opened = webbrowser.open_new_tab(viewer_file.resolve().as_uri())

    return AnalysisRun(
        conformers=conformer_tuple,
        all_atom_matrix=all_atom_matrix,
        heavy_atom_matrix=heavy_atom_matrix,
        matrix_file=matrix_file,
        viewer_file=viewer_file,
        warnings=tuple(warnings),
        viewer_opened=viewer_opened,
    )


def format_rmsd_report(
    conformers: Sequence[Conformer],
    labels: Sequence[str],
    all_atom_matrix: RmsdMatrix,
    heavy_atom_matrix: RmsdMatrix,
    *,
    warnings: Sequence[str] = (),
) -> str:
    """Format fitted RMSD matrices and conformer provenance as plain text."""

    conformer_tuple = tuple(conformers)
    label_tuple = tuple(labels)
    if len(conformer_tuple) != len(label_tuple):
        raise ValueError("The number of labels must match the conformer count.")
    lines = [
        "Fitted RMSD matrices",
        "====================",
        "",
        "Conformers (matrix order):",
    ]
    lines.extend(
        f"{label}: {conformer.name} - {conformer.source_file}"
        for label, conformer in zip(label_tuple, conformer_tuple)
    )
    lines.extend(
        (
            "",
            "Units: angstrom (A)",
            "",
            "All atoms (fit and RMSD use all atoms)",
            _format_matrix(all_atom_matrix, label_tuple),
            "",
            "Heavy atoms (fit and RMSD exclude hydrogen)",
            _format_matrix(heavy_atom_matrix, label_tuple),
            "",
        )
    )
    if warnings:
        lines.extend(("Warnings:", *(f"- {warning}" for warning in warnings), ""))
    return "\n".join(lines)


def _format_matrix(matrix: RmsdMatrix, labels: tuple[str, ...]) -> str:
    if len(matrix) != len(labels) or any(len(row) != len(labels) for row in matrix):
        raise ValueError("RMSD matrix dimensions must match the conformer labels.")
    label_width = max(8, *(len(label) for label in labels))
    value_width = 15
    header = " " * label_width + "".join(
        f"{label:>{value_width}}" for label in labels
    )
    rows = [header]
    rows.extend(
        f"{label:<{label_width}}"
        + "".join(f"{value:>{value_width}.10f}" for value in row)
        for label, row in zip(labels, matrix)
    )
    return "\n".join(rows)


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="conformer-analyzer",
        description=(
            "Align conformers, write fitted RMSD matrices, and open "
            "an interactive py3Dmol overlay."
        ),
    )
    parser.add_argument(
        "logs",
        nargs="+",
        type=Path,
        help="Two or more Gaussian .log or XYZ .xyz files in comparison order.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory for rmsd_matrix.txt and interactive_overlay.html.",
    )
    parser.add_argument(
        "--viewer-fit-atoms",
        choices=("all", "heavy"),
        default="all",
        help="Atoms used to fit conformers in the 3D viewer (default: all).",
    )
    parser.add_argument(
        "--show-hydrogens",
        action="store_true",
        help="Show hydrogens initially; the viewer toolbar can toggle them.",
    )
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Generate the viewer without opening a browser tab.",
    )
    return parser


def main() -> None:
    """Run the command-line interface."""

    parser = _argument_parser()
    arguments = parser.parse_args()
    try:
        result = analyze_log_files(
            arguments.logs,
            arguments.output_dir,
            viewer_fit_atoms=arguments.viewer_fit_atoms,
            show_hydrogens=arguments.show_hydrogens,
            open_viewer=not arguments.no_open,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(result.matrix_file.read_text(encoding="utf-8"))
    print(f"RMSD report: {result.matrix_file.resolve()}")
    print(f"Interactive viewer: {result.viewer_file.resolve()}")
    if not arguments.no_open and not result.viewer_opened:
        print("The browser could not be opened automatically; open the HTML file manually.")


if __name__ == "__main__":
    main()
