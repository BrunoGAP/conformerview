"""Root-mean-square deviation calculations after rigid-body fitting."""

from math import fsum, sqrt
from typing import Sequence

from conformer_analyzer.alignment import align_geometry
from conformer_analyzer.models import Conformer, Coordinates, MoleculeGeometry
from conformer_analyzer.validation import (
    heavy_atom_geometry,
    validate_comparable_geometries,
)


def _coordinate_rmsd(
    reference: Sequence[Coordinates],
    candidate: Sequence[Coordinates],
) -> float:
    """Calculate RMSD between corresponding Cartesian coordinate triples.

    No translation or rotation is performed. Coordinates are compared exactly
    in the supplied order and reference frame.
    """

    if len(reference) != len(candidate):
        raise ValueError(
            "Coordinate sets must contain the same number of points "
            f"(reference={len(reference)}, candidate={len(candidate)})."
        )
    if not reference:
        raise ValueError("RMSD requires at least one coordinate pair.")

    squared_distances = (
        (reference_x - candidate_x) ** 2
        + (reference_y - candidate_y) ** 2
        + (reference_z - candidate_z) ** 2
        for (reference_x, reference_y, reference_z), (
            candidate_x,
            candidate_y,
            candidate_z,
        ) in zip(reference, candidate)
    )
    return sqrt(fsum(squared_distances) / len(reference))


def _calculate_direct_rmsd(
    reference: MoleculeGeometry,
    candidate: MoleculeGeometry,
    *,
    exclude_hydrogens: bool = False,
) -> float:
    """Calculate direct atom-by-atom RMSD for two molecular geometries.

    Comparability is checked before optional hydrogen removal, ensuring a
    mismatched full atom sequence cannot be hidden by filtering.
    """

    validate_comparable_geometries(reference, candidate)

    if exclude_hydrogens:
        reference = heavy_atom_geometry(reference)
        candidate = heavy_atom_geometry(candidate)
        if not reference.atoms:
            raise ValueError("RMSD requires at least one non-hydrogen atom.")

    return _coordinate_rmsd(
        tuple(atom.coordinates for atom in reference.atoms),
        tuple(atom.coordinates for atom in candidate.atoms),
    )


def calculate_rmsd(
    reference: MoleculeGeometry,
    candidate: MoleculeGeometry,
    *,
    exclude_hydrogens: bool = False,
) -> float:
    """Calculate RMSD after optimal rigid-body alignment of the candidate."""

    aligned_candidate = align_geometry(
        reference,
        candidate,
        exclude_hydrogens=exclude_hydrogens,
    )
    return _calculate_direct_rmsd(
        reference,
        aligned_candidate,
        exclude_hydrogens=exclude_hydrogens,
    )


def _validate_geometry_sequence(
    geometries: Sequence[MoleculeGeometry],
) -> None:
    if not geometries:
        return
    reference = geometries[0]
    for candidate in geometries[1:]:
        validate_comparable_geometries(reference, candidate)


def _normalize_geometries(
    items: Sequence[MoleculeGeometry | Conformer],
) -> tuple[MoleculeGeometry, ...]:
    """Return geometries from raw geometries or parsed conformers."""

    geometries: list[MoleculeGeometry] = []
    for index, item in enumerate(items):
        if isinstance(item, Conformer):
            geometries.append(item.geometry)
        elif isinstance(item, MoleculeGeometry):
            geometries.append(item)
        else:
            raise TypeError(
                "Pairwise RMSD inputs must be MoleculeGeometry or Conformer "
                f"objects (item {index + 1} is {type(item).__name__})."
            )
    return tuple(geometries)


def _pairwise_matrix(
    geometries: Sequence[MoleculeGeometry],
    rmsd_function,
) -> tuple[tuple[float, ...], ...]:
    if not geometries:
        return ()

    n = len(geometries)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = rmsd_function(geometries[i], geometries[j])
            matrix[i][j] = value
            matrix[j][i] = value

    return tuple(tuple(row) for row in matrix)


def pairwise_rmsd_matrix(
    geometries: Sequence[MoleculeGeometry | Conformer],
    *,
    exclude_hydrogens: bool = False,
) -> tuple[tuple[float, ...], ...]:
    """Return a fitted RMSD matrix for geometries or parsed conformers."""
    normalized_geometries = _normalize_geometries(geometries)
    _validate_geometry_sequence(normalized_geometries)

    if exclude_hydrogens:
        rmsd_function = lambda reference, candidate: calculate_rmsd(
            reference,
            candidate,
            exclude_hydrogens=True,
        )
    else:
        rmsd_function = calculate_rmsd

    return _pairwise_matrix(normalized_geometries, rmsd_function)
