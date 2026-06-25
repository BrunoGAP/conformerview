"""Read supported conformer files and extract molecular geometries."""

from dataclasses import dataclass
from pathlib import Path
import re

from conformer_analyzer.elements import SYMBOL_TO_ATOMIC_NUMBER
from conformer_analyzer.models import (
    Atom,
    Conformer,
    GaussianJobSegment,
    GaussianJobType,
    GeometrySource,
    MoleculeGeometry,
    ParseResult,
    ParseWarning,
    ParseWarningCode,
)


_LINK_ONE_RE = re.compile(r"^\s*Entering Link 1\s*=", re.IGNORECASE)
_ROUTE_START_RE = re.compile(r"^\s*#")
_SEPARATOR_RE = re.compile(r"^\s*-{5,}\s*$")
_OPTIMIZATION_RE = re.compile(r"(?<![a-z0-9])opt(?=$|[^a-z0-9])", re.IGNORECASE)
_FREQUENCY_RE = re.compile(r"(?<![a-z0-9])freq(?=$|[^a-z0-9])", re.IGNORECASE)
_EXCITED_STATE_RE = re.compile(
    r"(?<![a-z0-9])(?:td|tda|cis|eomccsd|zindo)(?=$|[^a-z0-9])",
    re.IGNORECASE,
)
_CHARGE_MULTIPLICITY_RE = re.compile(
    r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)",
    re.IGNORECASE,
)
_SCF_ENERGY_RE = re.compile(
    r"SCF Done:.*?=\s*([-+]?\d+(?:\.\d*)?(?:[DE][+-]?\d+)?)",
    re.IGNORECASE,
)


def parse_conformer_file(source_file: str | Path) -> ParseResult:
    """Parse one supported conformer input file."""

    path = Path(source_file)
    suffix = path.suffix.lower()
    if suffix == ".log":
        return parse_gaussian_log(path)
    if suffix == ".xyz":
        return parse_xyz_file(path)
    raise ValueError(f"Unsupported input file extension for {path}: expected .log or .xyz.")


@dataclass(frozen=True, slots=True)
class _OrientationBlock:
    line_number: int
    kind: str
    geometry: MoleculeGeometry


def parse_gaussian_log(source_file: str | Path) -> ParseResult:
    """Parse conformers and warnings from one Gaussian log file."""

    path = Path(source_file)
    return parse_gaussian_text(path.read_text(encoding="utf-8"), source_file=path)


def parse_gaussian_text(
    log_text: str,
    source_file: str | Path = "<memory>",
) -> ParseResult:
    """Parse Gaussian text, applying the documented TD-input fallback."""

    path = Path(source_file)
    jobs = segment_gaussian_jobs(log_text)
    selected: list[tuple[GaussianJobSegment, MoleculeGeometry]] = []

    for job in jobs:
        if not job.is_ground_state_optimization_candidate:
            continue
        if not _is_successful_optimization(job):
            continue
        geometry = _final_optimized_geometry(job)
        if geometry is not None:
            selected.append((job, geometry))

    warnings: list[ParseWarning] = []
    if len(selected) > 1:
        warnings.append(
            ParseWarning(
                code=ParseWarningCode.MULTIPLE_GROUND_STATE_OPTIMIZATIONS,
                message=(
                    f"Multiple successfully optimized ground-state jobs were found "
                    f"in {path.name}; all eligible geometries were extracted."
                ),
            )
        )

    if selected:
        conformers = tuple(
            _build_conformer(
                path,
                job,
                geometry,
                GeometrySource.OPTIMIZED_GROUND_STATE,
                include_energy=True,
                use_job_suffix=len(selected) > 1,
            )
            for job, geometry in selected
        )
        return ParseResult(path, conformers, tuple(warnings))

    fallback = _find_td_input_fallback(jobs)
    if fallback is not None:
        job, geometry = fallback
        warnings.append(
            ParseWarning(
                code=ParseWarningCode.TD_INPUT_GEOMETRY_FALLBACK,
                job_index=job.index,
                line_number=job.start_line,
                message=(
                    "No explicit successfully optimized ground-state geometry "
                    f"was found in {path.name}. The input geometry from TD "
                    f"optimization job {job.index} will be used; this assumes it "
                    "was previously optimized at the ground state."
                ),
            )
        )
        conformer = _build_conformer(
            path,
            job,
            geometry,
            GeometrySource.TD_OPTIMIZATION_INPUT,
            include_energy=False,
            use_job_suffix=False,
        )
        return ParseResult(path, (conformer,), tuple(warnings))

    warnings.append(
        ParseWarning(
            code=ParseWarningCode.NO_USABLE_GEOMETRY,
            message=(
                f"No successfully optimized ground-state geometry or TD "
                f"optimization input geometry was found in {path.name}."
            ),
        )
    )
    return ParseResult(path, (), tuple(warnings))


def parse_xyz_file(source_file: str | Path) -> ParseResult:
    """Parse a simple single-geometry XYZ file."""

    path = Path(source_file)
    return parse_xyz_text(path.read_text(encoding="utf-8"), source_file=path)


def parse_xyz_text(
    xyz_text: str,
    source_file: str | Path = "<memory>",
) -> ParseResult:
    """Parse atom symbols and Cartesian coordinates from simple XYZ text."""

    path = Path(source_file)
    lines = xyz_text.splitlines()
    if not lines:
        raise ValueError(f"XYZ file {path.name} is empty.")

    try:
        atom_count = int(lines[0].strip())
    except ValueError as error:
        raise ValueError(
            f"XYZ file {path.name} must start with an integer atom count."
        ) from error
    if atom_count <= 0:
        raise ValueError(f"XYZ file {path.name} must contain at least one atom.")
    if len(lines) < atom_count + 2:
        raise ValueError(
            f"XYZ file {path.name} declares {atom_count} atom(s) but contains "
            f"only {max(0, len(lines) - 2)} coordinate row(s)."
        )

    atoms: list[Atom] = []
    for atom_index, row in enumerate(lines[2 : 2 + atom_count], start=1):
        columns = row.split()
        if len(columns) < 4:
            raise ValueError(
                f"Malformed XYZ coordinate row {atom_index} in {path.name}: "
                "expected an element symbol followed by x, y, and z coordinates."
            )
        symbol = _canonicalize_symbol(columns[0])
        try:
            atomic_number = SYMBOL_TO_ATOMIC_NUMBER[symbol.lower()]
        except KeyError as error:
            raise ValueError(
                f"Invalid element symbol '{columns[0]}' on XYZ row {atom_index} "
                f"in {path.name}."
            ) from error
        try:
            coordinates = (
                float(columns[1]),
                float(columns[2]),
                float(columns[3]),
            )
        except ValueError as error:
            raise ValueError(
                f"Malformed XYZ coordinate row {atom_index} in {path.name}: "
                "x, y, and z coordinates must be numeric."
            ) from error
        atoms.append(
            Atom(
                atomic_number=atomic_number,
                coordinates=coordinates,
                symbol=symbol,
            )
        )

    trailing_rows = [row for row in lines[2 + atom_count :] if row.strip()]
    if trailing_rows:
        raise ValueError(
            f"XYZ file {path.name} declares {atom_count} atom(s) but contains "
            "additional coordinate rows."
        )

    conformer = Conformer(
        name=path.stem,
        source_file=path,
        geometry=MoleculeGeometry(tuple(atoms)),
        geometry_source=GeometrySource.XYZ_FILE,
    )
    return ParseResult(path, (conformer,), ())


def read_gaussian_jobs(source_file: str | Path) -> tuple[GaussianJobSegment, ...]:
    """Read a Gaussian log file and return its jobs in source order."""

    path = Path(source_file)
    return segment_gaussian_jobs(path.read_text(encoding="utf-8"))


def _build_conformer(
    path: Path,
    job: GaussianJobSegment,
    geometry: MoleculeGeometry,
    geometry_source: GeometrySource,
    *,
    include_energy: bool,
    use_job_suffix: bool,
) -> Conformer:
    charge, multiplicity = _extract_charge_and_multiplicity(job.text)
    name = f"{path.stem}_job{job.index}" if use_job_suffix else path.stem
    return Conformer(
        name=name,
        source_file=path,
        geometry=geometry,
        charge=charge,
        multiplicity=multiplicity,
        energy_hartree=_extract_last_scf_energy(job.text) if include_energy else None,
        source_job_index=job.index,
        geometry_source=geometry_source,
    )


def _find_td_input_fallback(
    jobs: tuple[GaussianJobSegment, ...],
) -> tuple[GaussianJobSegment, MoleculeGeometry] | None:
    for job in jobs:
        if GaussianJobType.OPTIMIZATION not in job.job_types or not job.is_excited_state:
            continue
        blocks = _parse_orientation_blocks(job)
        if blocks:
            return job, blocks[0].geometry
    return None


def _is_successful_optimization(job: GaussianJobSegment) -> bool:
    lowered = job.text.lower()
    converged = (
        "optimization completed." in lowered
        or "-- stationary point found." in lowered
    )
    normally_terminated = "normal termination of gaussian" in lowered
    errored = "error termination" in lowered
    return converged and normally_terminated and not errored


def _final_optimized_geometry(job: GaussianJobSegment) -> MoleculeGeometry | None:
    blocks = _parse_orientation_blocks(job)
    if not blocks:
        return None
    standard_blocks = [block for block in blocks if block.kind == "standard"]
    return (standard_blocks or blocks)[-1].geometry


def _extract_charge_and_multiplicity(text: str) -> tuple[int | None, int | None]:
    match = _CHARGE_MULTIPLICITY_RE.search(text)
    if match is None:
        return None, None
    return int(match.group(1)), int(match.group(2))


def _extract_last_scf_energy(text: str) -> float | None:
    matches = _SCF_ENERGY_RE.findall(text)
    if not matches:
        return None
    return float(matches[-1].replace("D", "E").replace("d", "e"))


def _parse_orientation_blocks(job: GaussianJobSegment) -> tuple[_OrientationBlock, ...]:
    lines = job.text.splitlines()
    blocks: list[_OrientationBlock] = []

    for marker_index, line in enumerate(lines):
        normalized = line.strip().lower()
        if normalized not in {"standard orientation:", "input orientation:"}:
            continue

        separator_count = 0
        data_start: int | None = None
        for line_index in range(marker_index + 1, len(lines)):
            if _SEPARATOR_RE.match(lines[line_index]):
                separator_count += 1
                if separator_count == 2:
                    data_start = line_index + 1
                    break
        if data_start is None:
            continue

        atoms: list[Atom] = []
        valid_block = True
        for row in lines[data_start:]:
            if _SEPARATOR_RE.match(row):
                break
            columns = row.split()
            if len(columns) < 6:
                valid_block = False
                break
            try:
                atoms.append(
                    Atom(
                        atomic_number=int(columns[1]),
                        coordinates=(
                            _parse_gaussian_float(columns[3]),
                            _parse_gaussian_float(columns[4]),
                            _parse_gaussian_float(columns[5]),
                        ),
                    )
                )
            except (ValueError, IndexError):
                valid_block = False
                break

        if valid_block and atoms:
            blocks.append(
                _OrientationBlock(
                    line_number=job.start_line + marker_index,
                    kind="standard" if normalized.startswith("standard") else "input",
                    geometry=MoleculeGeometry(tuple(atoms)),
                )
            )

    return tuple(blocks)


def _parse_gaussian_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _canonicalize_symbol(symbol: str) -> str:
    stripped = symbol.strip()
    if not stripped:
        return stripped
    return stripped[0].upper() + stripped[1:].lower()


def segment_gaussian_jobs(log_text: str) -> tuple[GaussianJobSegment, ...]:
    """Split Gaussian output text and classify each job's route section.

    Gaussian prints ``Entering Link 1`` once for each top-level job. Appended
    jobs in the same log usually begin with a new ``Initial command`` block, so
    that line is included when it precedes a later Link 1 marker.
    """

    lines = log_text.splitlines(keepends=True)
    if not lines:
        return ()

    starts = _find_job_starts(lines)
    jobs: list[GaussianJobSegment] = []

    for index, start in enumerate(starts, start=1):
        end = starts[index] if index < len(starts) else len(lines)
        segment_lines = lines[start:end]
        route_text = _extract_route_text(segment_lines)
        job_types = _detect_job_types(route_text)

        jobs.append(
            GaussianJobSegment(
                index=index,
                start_line=start + 1,
                end_line=end,
                text="".join(segment_lines),
                route_text=route_text,
                job_types=job_types,
                is_excited_state=bool(_EXCITED_STATE_RE.search(route_text)),
            )
        )

    return tuple(jobs)


def _find_job_starts(lines: list[str]) -> list[int]:
    link_starts = [
        line_index
        for line_index, line in enumerate(lines)
        if _LINK_ONE_RE.match(line)
    ]
    if len(link_starts) < 2:
        return [0]

    starts = [0]
    previous_link = link_starts[0]

    for link_start in link_starts[1:]:
        job_start = link_start
        for line_index in range(link_start - 1, previous_link, -1):
            if lines[line_index].strip().lower() == "initial command:":
                job_start = line_index
                break
        starts.append(job_start)
        previous_link = link_start

    return starts


def _extract_route_text(lines: list[str]) -> str:
    route_start = next(
        (
            line_index
            for line_index, line in enumerate(lines)
            if _ROUTE_START_RE.match(line)
        ),
        None,
    )
    if route_start is None:
        return ""

    route_lines: list[str] = []
    for line in lines[route_start:]:
        line_without_ending = line.rstrip("\r\n")
        if route_lines and (
            _SEPARATOR_RE.match(line_without_ending) or not line_without_ending.strip()
        ):
            break
        route_lines.append(line_without_ending)

    # Gaussian prefixes each displayed route line with one formatting space.
    # Removing only that space preserves a real separator while reconnecting a
    # keyword split at the fixed output width, such as "fre" + "q".
    unwrapped = "".join(
        line[1:] if line.startswith(" ") else line for line in route_lines
    )
    return " ".join(unwrapped.split())


def _detect_job_types(route_text: str) -> frozenset[GaussianJobType]:
    if not route_text:
        return frozenset()

    job_types: set[GaussianJobType] = set()
    if _OPTIMIZATION_RE.search(route_text):
        job_types.add(GaussianJobType.OPTIMIZATION)
    if _FREQUENCY_RE.search(route_text):
        job_types.add(GaussianJobType.FREQUENCY)
    if not job_types:
        job_types.add(GaussianJobType.SINGLE_POINT)

    return frozenset(job_types)
