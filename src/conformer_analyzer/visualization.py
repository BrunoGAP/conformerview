"""Render aligned conformers as PyMOL-inspired static overlay images."""

from dataclasses import dataclass
from math import ceil
import os
from pathlib import Path
import tempfile
from typing import Sequence, TypeAlias

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "conformer-analyzer-matplotlib"),
)

import matplotlib

matplotlib.use("Agg")

from matplotlib import colors as matplotlib_colors
from matplotlib import pyplot as plt
import numpy as np
from numpy.typing import NDArray

from conformer_analyzer.alignment import align_geometry
from conformer_analyzer.models import Conformer, MoleculeGeometry
from conformer_analyzer.validation import validate_comparable_geometries


Bond: TypeAlias = tuple[int, int]
CoordinateArray = NDArray[np.float64]


DEFAULT_COLORS = (
    "#ef2020",
    "#2436f5",
    "#27833a",
    "#ff7a16",
    "#ed16d2",
    "#18ced1",
)

# Single-bond covalent radii in angstroms for elements currently encountered
# in the Gaussian examples. Unsupported elements fail explicitly so uncertain
# connectivity is never silently invented.
_COVALENT_RADII = {
    1: 0.31,
    6: 0.76,
    7: 0.71,
    8: 0.66,
    9: 0.57,
    15: 1.07,
    16: 1.05,
    17: 1.02,
    35: 1.20,
    53: 1.39,
}


@dataclass(frozen=True, slots=True)
class OverlayExportResult:
    """Metadata and warnings from one overlay image export."""

    output_file: Path
    bonds: tuple[Bond, ...]
    warnings: tuple[str, ...]


def infer_bonds(
    geometry: MoleculeGeometry,
    *,
    tolerance: float = 1.20,
    minimum_distance: float = 0.40,
) -> tuple[Bond, ...]:
    """Infer conservative connectivity from interatomic covalent radii."""

    if tolerance <= 0.0:
        raise ValueError("Bond tolerance must be positive.")
    unsupported = sorted(
        {atom.atomic_number for atom in geometry.atoms} - _COVALENT_RADII.keys()
    )
    if unsupported:
        numbers = ", ".join(str(number) for number in unsupported)
        raise ValueError(
            "Cannot infer connectivity: no covalent radius is configured for "
            f"atomic number(s) {numbers}."
        )

    bonds: list[Bond] = []
    for first_index, first_atom in enumerate(geometry.atoms):
        for second_index in range(first_index + 1, len(geometry.atoms)):
            second_atom = geometry.atoms[second_index]
            distance = float(
                np.linalg.norm(
                    np.asarray(first_atom.coordinates)
                    - np.asarray(second_atom.coordinates)
                )
            )
            cutoff = tolerance * (
                _COVALENT_RADII[first_atom.atomic_number]
                + _COVALENT_RADII[second_atom.atomic_number]
            )
            if minimum_distance <= distance <= cutoff:
                bonds.append((first_index, second_index))
    return tuple(bonds)


def align_conformers_for_overlay(
    conformers: Sequence[Conformer],
    *,
    exclude_hydrogens: bool = False,
) -> tuple[MoleculeGeometry, ...]:
    """Align every conformer onto the first conformer's geometry."""

    if not conformers:
        raise ValueError("Overlay rendering requires at least one conformer.")
    reference = conformers[0].geometry
    aligned = [reference]
    for conformer in conformers[1:]:
        validate_comparable_geometries(reference, conformer.geometry)
        aligned.append(
            align_geometry(
                reference,
                conformer.geometry,
                exclude_hydrogens=exclude_hydrogens,
            )
        )
    return tuple(aligned)


def export_overlay_png(
    conformers: Sequence[Conformer],
    output_file: str | Path,
    *,
    labels: Sequence[str] | None = None,
    percentages: Sequence[float | None] | None = None,
    colors: Sequence[str] | None = None,
    exclude_hydrogens_from_fit: bool = False,
    show_hydrogens: bool = False,
    width: int = 1600,
    height: int = 900,
    dpi: int = 100,
) -> OverlayExportResult:
    """Align conformers and export a white-background molecular overlay PNG."""

    if width <= 0 or height <= 0 or dpi <= 0:
        raise ValueError("Image width, height, and DPI must be positive.")

    conformer_tuple = tuple(conformers)
    aligned_geometries = align_conformers_for_overlay(
        conformer_tuple,
        exclude_hydrogens=exclude_hydrogens_from_fit,
    )
    display_labels = _validated_labels(conformer_tuple, labels)
    display_percentages = _validated_percentages(
        len(conformer_tuple),
        percentages,
    )
    display_colors = _validated_colors(len(conformer_tuple), colors)

    reference = aligned_geometries[0]
    bonds = infer_bonds(reference)
    warnings = [
        "Connectivity was inferred from covalent radii and should be verified "
        "when explicit bond data becomes available."
    ]
    bonded_indices = {atom_index for bond in bonds for atom_index in bond}
    visible_indices = tuple(
        index
        for index, atom in enumerate(reference.atoms)
        if show_hydrogens or atom.atomic_number != 1
    )
    isolated = [index + 1 for index in visible_indices if index not in bonded_indices]
    if isolated:
        warnings.append(
            "No bond was inferred for visible atom position(s): "
            + ", ".join(str(index) for index in isolated)
            + "."
        )

    visible_bonds = tuple(
        (first, second)
        for first, second in bonds
        if first in visible_indices and second in visible_indices
    )
    if not visible_bonds:
        raise ValueError("No visible bonds could be inferred for overlay rendering.")

    projected = _project_geometries(aligned_geometries, visible_indices)
    figure = _render_figure(
        projected,
        visible_indices,
        visible_bonds,
        display_labels,
        display_percentages,
        display_colors,
        width=width,
        height=height,
        dpi=dpi,
    )

    path = Path(output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        path,
        format="png",
        dpi=dpi,
        facecolor="white",
        edgecolor="white",
    )
    plt.close(figure)
    return OverlayExportResult(path, bonds, tuple(warnings))


def _validated_labels(
    conformers: tuple[Conformer, ...],
    labels: Sequence[str] | None,
) -> tuple[str, ...]:
    values = tuple(labels) if labels is not None else tuple(
        conformer.name for conformer in conformers
    )
    if len(values) != len(conformers):
        raise ValueError("The number of labels must match the conformer count.")
    if any(not label.strip() for label in values):
        raise ValueError("Overlay labels cannot be empty.")
    return values


def _validated_percentages(
    conformer_count: int,
    percentages: Sequence[float | None] | None,
) -> tuple[float | None, ...]:
    values = (
        tuple(percentages)
        if percentages is not None
        else (None,) * conformer_count
    )
    if len(values) != conformer_count:
        raise ValueError(
            "The number of percentages must match the conformer count."
        )
    if any(value is not None and not 0.0 <= value <= 100.0 for value in values):
        raise ValueError("Overlay percentages must be between 0 and 100.")
    return values


def _validated_colors(
    conformer_count: int,
    colors: Sequence[str] | None,
) -> tuple[str, ...]:
    if colors is None:
        repeats = ceil(conformer_count / len(DEFAULT_COLORS))
        values = (DEFAULT_COLORS * repeats)[:conformer_count]
    else:
        values = tuple(colors)
    if len(values) != conformer_count:
        raise ValueError("The number of colors must match the conformer count.")
    try:
        for color in values:
            matplotlib_colors.to_rgba(color)
    except ValueError as error:
        raise ValueError(f"Invalid overlay color: {error}") from error
    return values


def _project_geometries(
    geometries: tuple[MoleculeGeometry, ...],
    visible_indices: tuple[int, ...],
) -> tuple[CoordinateArray, ...]:
    visible_coordinates = np.asarray(
        [
            geometry.atoms[index].coordinates
            for geometry in geometries
            for index in visible_indices
        ],
        dtype=np.float64,
    )
    center = visible_coordinates.mean(axis=0)
    _, _, principal_axes = np.linalg.svd(
        visible_coordinates - center,
        full_matrices=False,
    )
    projection = principal_axes.T
    if projection[np.argmax(np.abs(projection[:, 0])), 0] < 0.0:
        projection[:, 0] *= -1.0
    if np.linalg.det(projection) < 0.0:
        projection[:, 2] *= -1.0

    return tuple(
        (np.asarray([atom.coordinates for atom in geometry.atoms]) - center)
        @ projection
        for geometry in geometries
    )


def _render_figure(
    projected: tuple[CoordinateArray, ...],
    visible_indices: tuple[int, ...],
    bonds: tuple[Bond, ...],
    labels: tuple[str, ...],
    percentages: tuple[float | None, ...],
    colors: tuple[str, ...],
    *,
    width: int,
    height: int,
    dpi: int,
):
    legend_rows = ceil(len(labels) / min(3, len(labels)))
    legend_fraction = min(0.28, 0.08 + 0.08 * legend_rows)
    figure = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    figure.patch.set_facecolor("white")
    molecule_axis = figure.add_axes(
        (0.025, legend_fraction, 0.95, 0.97 - legend_fraction)
    )
    legend_axis = figure.add_axes((0.02, 0.01, 0.96, legend_fraction - 0.02))
    molecule_axis.set_facecolor("white")
    legend_axis.set_facecolor("white")

    line_records = []
    for conformer_index, coordinates in enumerate(projected):
        for first, second in bonds:
            depth = float((coordinates[first, 2] + coordinates[second, 2]) / 2.0)
            line_records.append(
                (depth, conformer_index, coordinates[first], coordinates[second])
            )
    for order, (_, conformer_index, first, second) in enumerate(
        sorted(line_records, key=lambda record: record[0])
    ):
        color = colors[conformer_index]
        dark_color = _darken(color, 0.56)
        molecule_axis.plot(
            (first[0], second[0]),
            (first[1], second[1]),
            color=dark_color,
            linewidth=10.5,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2 * order,
        )
        molecule_axis.plot(
            (first[0], second[0]),
            (first[1], second[1]),
            color=color,
            linewidth=6.8,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2 * order + 1,
        )

    all_visible = np.concatenate(
        [coordinates[list(visible_indices), :2] for coordinates in projected]
    )
    minimum = all_visible.min(axis=0)
    maximum = all_visible.max(axis=0)
    span = np.maximum(maximum - minimum, 1.0)
    molecule_axis.set_xlim(
        minimum[0] - 0.05 * span[0],
        maximum[0] + 0.05 * span[0],
    )
    molecule_axis.set_ylim(
        minimum[1] - 0.08 * span[1],
        maximum[1] + 0.08 * span[1],
    )
    molecule_axis.set_aspect("equal", adjustable="box")
    molecule_axis.axis("off")

    columns = min(3, len(labels))
    for index, (label, percentage, color) in enumerate(
        zip(labels, percentages, colors)
    ):
        row = index // columns
        column = index % columns
        x = (column + 0.08) / columns
        y = 1.0 - (row + 0.5) / legend_rows
        legend_axis.plot(
            (x, x + 0.075),
            (y, y),
            color=color,
            linewidth=8.0,
            solid_capstyle="butt",
            transform=legend_axis.transAxes,
        )
        legend_text = label
        if percentage is not None:
            formatted = f"{percentage:g}"
            legend_text = f"{label} ({formatted}%)"
        legend_axis.text(
            x + 0.095,
            y,
            legend_text,
            transform=legend_axis.transAxes,
            va="center",
            ha="left",
            fontsize=22,
            color="black",
        )
    legend_axis.set_xlim(0.0, 1.0)
    legend_axis.set_ylim(0.0, 1.0)
    legend_axis.axis("off")
    return figure


def _darken(color: str, factor: float) -> tuple[float, float, float]:
    red, green, blue, _ = matplotlib_colors.to_rgba(color)
    return red * factor, green * factor, blue * factor
