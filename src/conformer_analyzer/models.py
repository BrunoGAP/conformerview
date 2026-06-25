"""Typed data structures shared across the application.

Coordinates use angstroms and electronic energies use Hartree, matching the
units normally printed in the Gaussian output sections this project will read.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TypeAlias


Coordinates: TypeAlias = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class Atom:
    """One atom at a Cartesian position.

    ``atomic_number`` is the canonical atom identity. The optional symbol is
    retained for display and input formats but is not used to reorder atoms.
    """

    atomic_number: int
    coordinates: Coordinates
    symbol: str | None = None

    def __post_init__(self) -> None:
        if self.atomic_number < 1:
            raise ValueError("atomic_number must be a positive integer")
        if len(self.coordinates) != 3:
            raise ValueError("coordinates must contain exactly three values")

    @property
    def x(self) -> float:
        """Return the x coordinate in angstroms."""

        return self.coordinates[0]

    @property
    def y(self) -> float:
        """Return the y coordinate in angstroms."""

        return self.coordinates[1]

    @property
    def z(self) -> float:
        """Return the z coordinate in angstroms."""

        return self.coordinates[2]


@dataclass(frozen=True, slots=True)
class MoleculeGeometry:
    """A molecular geometry whose atom tuple defines correspondence order."""

    atoms: tuple[Atom, ...]


class GeometrySource(str, Enum):
    """How a conformer's geometry was selected from its input file."""

    OPTIMIZED_GROUND_STATE = "optimized_ground_state"
    TD_OPTIMIZATION_INPUT = "td_optimization_input"
    XYZ_FILE = "xyz_file"


@dataclass(frozen=True, slots=True)
class Conformer:
    """A molecular conformer and the calculation metadata available for it."""

    name: str
    source_file: Path
    geometry: MoleculeGeometry
    charge: int | None = None
    multiplicity: int | None = None
    energy_hartree: float | None = None
    source_job_index: int | None = None
    geometry_source: GeometrySource | None = None

    def __post_init__(self) -> None:
        if self.multiplicity is not None and self.multiplicity < 1:
            raise ValueError("multiplicity must be a positive integer")


class ParseWarningCode(str, Enum):
    """Stable identifiers that an interface can use for warning behavior."""

    MULTIPLE_GROUND_STATE_OPTIMIZATIONS = "multiple_ground_state_optimizations"
    TD_INPUT_GEOMETRY_FALLBACK = "td_input_geometry_fallback"
    NO_USABLE_GEOMETRY = "no_usable_geometry"


@dataclass(frozen=True, slots=True)
class ParseWarning:
    """A non-fatal parsing issue with optional Gaussian location context."""

    message: str
    job_index: int | None = None
    line_number: int | None = None
    code: ParseWarningCode | None = None


@dataclass(frozen=True, slots=True)
class ParseResult:
    """Conformers and non-fatal warnings produced while parsing one file."""

    source_file: Path
    conformers: tuple[Conformer, ...]
    warnings: tuple[ParseWarning, ...] = ()


class GaussianJobType(str, Enum):
    """Calculation types detected from a Gaussian route section."""

    OPTIMIZATION = "optimization"
    FREQUENCY = "frequency"
    SINGLE_POINT = "single_point"


@dataclass(frozen=True, slots=True)
class GaussianJobSegment:
    """One Gaussian job and the route-level facts known about it.

    Line numbers and ``index`` are one-based for readable warnings. Route-level
    candidacy does not imply that an optimization converged; convergence is
    evaluated separately when geometries are extracted.
    """

    index: int
    start_line: int
    end_line: int
    text: str
    route_text: str
    job_types: frozenset[GaussianJobType]
    is_excited_state: bool

    @property
    def is_ground_state_optimization_candidate(self) -> bool:
        """Return whether the route describes a non-excited optimization."""

        return (
            GaussianJobType.OPTIMIZATION in self.job_types
            and not self.is_excited_state
        )
