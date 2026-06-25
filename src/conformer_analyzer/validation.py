"""Check whether molecular geometries can be compared safely."""

from conformer_analyzer.models import MoleculeGeometry


class GeometryComparisonError(ValueError):
    """Raised when atom-by-atom geometry correspondence is invalid."""


def atomic_number_sequence(geometry: MoleculeGeometry) -> tuple[int, ...]:
    """Return atomic numbers in the geometry's original atom order."""

    return tuple(atom.atomic_number for atom in geometry.atoms)


def validate_comparable_geometries(
    reference: MoleculeGeometry,
    candidate: MoleculeGeometry,
) -> None:
    """Ensure two geometries have one-to-one atom correspondence.

    Validation never reorders atoms. A successful return means each atom can be
    compared with the atom at the same tuple position in the other geometry.
    """

    reference_count = len(reference.atoms)
    candidate_count = len(candidate.atoms)
    if reference_count != candidate_count:
        raise GeometryComparisonError(
            "Geometries are not comparable: atom counts differ "
            f"(reference={reference_count}, candidate={candidate_count})."
        )

    for atom_index, (reference_atom, candidate_atom) in enumerate(
        zip(reference.atoms, candidate.atoms),
        start=1,
    ):
        if reference_atom.atomic_number != candidate_atom.atomic_number:
            raise GeometryComparisonError(
                "Geometries are not comparable: atomic-number sequences differ "
                f"at atom {atom_index} (reference={reference_atom.atomic_number}, "
                f"candidate={candidate_atom.atomic_number}). Atom reordering is "
                "not performed."
            )


def heavy_atom_geometry(geometry: MoleculeGeometry) -> MoleculeGeometry:
    """Return non-hydrogen atoms while preserving their relative order."""

    return MoleculeGeometry(
        atoms=tuple(atom for atom in geometry.atoms if atom.atomic_number != 1)
    )
