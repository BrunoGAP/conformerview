"""Optimal rigid-body alignment of molecular geometries."""

from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from conformer_analyzer.models import Atom, Coordinates, MoleculeGeometry
from conformer_analyzer.validation import (
    heavy_atom_geometry,
    validate_comparable_geometries,
)


CoordinateArray = NDArray[np.float64]


def centroid(coordinates: Sequence[Coordinates]) -> Coordinates:
    """Return the arithmetic centroid of Cartesian coordinates."""

    coordinate_array = _as_coordinate_array(coordinates)
    values = coordinate_array.mean(axis=0)
    return float(values[0]), float(values[1]), float(values[2])


def align_coordinates(
    reference: Sequence[Coordinates],
    candidate: Sequence[Coordinates],
) -> tuple[Coordinates, ...]:
    """Rigidly align candidate coordinates onto corresponding reference points."""

    reference_array, candidate_array = _validated_coordinate_arrays(
        reference,
        candidate,
    )
    rotation, reference_centroid, candidate_centroid = _fit_rigid_transform(
        reference_array,
        candidate_array,
    )
    aligned = (
        (candidate_array - candidate_centroid) @ rotation + reference_centroid
    )
    return _to_coordinates(aligned)


def align_geometry(
    reference: MoleculeGeometry,
    candidate: MoleculeGeometry,
    *,
    exclude_hydrogens: bool = False,
) -> MoleculeGeometry:
    """Return a rigidly aligned copy of ``candidate``.

    When hydrogens are excluded, the transform is fitted using corresponding
    heavy atoms and then applied unchanged to every atom in the candidate.
    """

    validate_comparable_geometries(reference, candidate)
    fit_reference = reference
    fit_candidate = candidate
    if exclude_hydrogens:
        fit_reference = heavy_atom_geometry(reference)
        fit_candidate = heavy_atom_geometry(candidate)
        if not fit_reference.atoms:
            raise ValueError("Alignment requires at least one non-hydrogen atom.")

    reference_array = _geometry_coordinates(fit_reference)
    candidate_array = _geometry_coordinates(fit_candidate)
    rotation, reference_centroid, candidate_centroid = _fit_rigid_transform(
        reference_array,
        candidate_array,
    )

    all_candidate_coordinates = _geometry_coordinates(candidate)
    aligned_coordinates = (
        (all_candidate_coordinates - candidate_centroid) @ rotation
        + reference_centroid
    )
    return MoleculeGeometry(
        atoms=tuple(
            Atom(
                atomic_number=atom.atomic_number,
                symbol=atom.symbol,
                coordinates=coordinates,
            )
            for atom, coordinates in zip(
                candidate.atoms,
                _to_coordinates(aligned_coordinates),
            )
        )
    )


def _validated_coordinate_arrays(
    reference: Sequence[Coordinates],
    candidate: Sequence[Coordinates],
) -> tuple[CoordinateArray, CoordinateArray]:
    if len(reference) != len(candidate):
        raise ValueError(
            "Coordinate sets must contain the same number of points "
            f"(reference={len(reference)}, candidate={len(candidate)})."
        )
    return _as_coordinate_array(reference), _as_coordinate_array(candidate)


def _as_coordinate_array(coordinates: Sequence[Coordinates]) -> CoordinateArray:
    if not coordinates:
        raise ValueError("Alignment requires at least one coordinate pair.")
    coordinate_array = np.asarray(coordinates, dtype=np.float64)
    if coordinate_array.shape != (len(coordinates), 3):
        raise ValueError("Each coordinate must contain exactly three values.")
    if not np.isfinite(coordinate_array).all():
        raise ValueError("Alignment coordinates must contain only finite values.")
    return coordinate_array


def _geometry_coordinates(geometry: MoleculeGeometry) -> CoordinateArray:
    return _as_coordinate_array(tuple(atom.coordinates for atom in geometry.atoms))


def _fit_rigid_transform(
    reference: CoordinateArray,
    candidate: CoordinateArray,
) -> tuple[CoordinateArray, CoordinateArray, CoordinateArray]:
    reference_centroid = reference.mean(axis=0)
    candidate_centroid = candidate.mean(axis=0)
    reference_centered = reference - reference_centroid
    candidate_centered = candidate - candidate_centroid

    covariance = candidate_centered.T @ reference_centered
    left_vectors, _, right_vectors_transposed = np.linalg.svd(covariance)
    rotation = left_vectors @ right_vectors_transposed

    # A negative determinant represents a reflection, which is not a valid
    # rigid-body rotation and would invert molecular chirality.
    if np.linalg.det(rotation) < 0.0:
        left_vectors[:, -1] *= -1.0
        rotation = left_vectors @ right_vectors_transposed

    return rotation, reference_centroid, candidate_centroid


def _to_coordinates(array: CoordinateArray) -> tuple[Coordinates, ...]:
    return tuple(
        (float(row[0]), float(row[1]), float(row[2]))
        for row in array
    )
