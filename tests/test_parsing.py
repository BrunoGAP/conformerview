"""Tests for Gaussian job segmentation and route classification."""

from pathlib import Path

import pytest

from conformer_analyzer.models import (
    GaussianJobType,
    GeometrySource,
    ParseWarningCode,
)
from conformer_analyzer.parsing import (
    parse_gaussian_log,
    parse_gaussian_text,
    parse_xyz_text,
    read_gaussian_jobs,
    segment_gaussian_jobs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFORMERS = PROJECT_ROOT / "conformers" / "COOH"


def _orientation_block(kind: str, x_coordinate: float) -> str:
    return f""" {kind} orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0       {x_coordinate:10.6f}    0.000000    0.000000
 ---------------------------------------------------------------------
"""


def test_conf1_contains_excited_optimization_and_frequency_jobs() -> None:
    jobs = read_gaussian_jobs(CONFORMERS / "conf1.log")

    assert len(jobs) == 2
    assert jobs[0].job_types == frozenset({GaussianJobType.OPTIMIZATION})
    assert jobs[1].job_types == frozenset({GaussianJobType.FREQUENCY})
    assert all(job.is_excited_state for job in jobs)
    assert not any(job.is_ground_state_optimization_candidate for job in jobs)
    assert "opt td" in jobs[0].route_text.lower()
    assert "freq td" in jobs[1].route_text.lower()


def test_conf2_contains_ground_and_excited_state_jobs() -> None:
    jobs = read_gaussian_jobs(CONFORMERS / "conf2.log")

    assert len(jobs) == 3
    assert [job.job_types for job in jobs] == [
        frozenset({GaussianJobType.OPTIMIZATION}),
        frozenset({GaussianJobType.OPTIMIZATION}),
        frozenset({GaussianJobType.FREQUENCY}),
    ]
    assert [job.is_excited_state for job in jobs] == [False, True, True]
    assert [job.is_ground_state_optimization_candidate for job in jobs] == [
        True,
        False,
        False,
    ]


def test_conf1_uses_td_input_geometry_with_user_warning() -> None:
    result = parse_gaussian_log(CONFORMERS / "conf1.log")

    assert len(result.conformers) == 1
    conformer = result.conformers[0]
    assert len(conformer.geometry.atoms) == 58
    assert conformer.charge == 0
    assert conformer.multiplicity == 1
    assert conformer.energy_hartree is None
    assert conformer.source_job_index == 1
    assert conformer.geometry_source is GeometrySource.TD_OPTIMIZATION_INPUT
    assert conformer.geometry.atoms[0].atomic_number == 8
    assert conformer.geometry.atoms[0].coordinates == pytest.approx(
        (-9.978252, -1.743074, -0.700694)
    )

    assert len(result.warnings) == 1
    warning = result.warnings[0]
    assert warning.code is ParseWarningCode.TD_INPUT_GEOMETRY_FALLBACK
    assert warning.job_index == 1
    assert "No explicit successfully optimized ground-state geometry" in warning.message
    assert "input geometry from TD optimization job 1" in warning.message
    assert "assumes it was previously optimized at the ground state" in warning.message


def test_conf2_prefers_explicit_optimized_ground_state_geometry() -> None:
    result = parse_gaussian_log(CONFORMERS / "conf2.log")

    assert result.warnings == ()
    assert len(result.conformers) == 1
    conformer = result.conformers[0]
    assert len(conformer.geometry.atoms) == 58
    assert conformer.charge == 0
    assert conformer.multiplicity == 1
    assert conformer.energy_hartree == pytest.approx(-1905.38861799)
    assert conformer.source_job_index == 1
    assert conformer.geometry_source is GeometrySource.OPTIMIZED_GROUND_STATE
    assert conformer.geometry.atoms[0].coordinates == pytest.approx(
        (10.142352, -0.664438, 0.338399)
    )


def test_td_fallback_prefers_first_input_orientation() -> None:
    log_text = (
        " Entering Link 1 = /gaussian/l1.exe\n"
        " -----\n"
        " # b3lyp/6-31g(d) opt td\n"
        " -----\n"
        " Charge = 0 Multiplicity = 1\n"
        + _orientation_block("Input", 1.25)
        + _orientation_block("Standard", 9.75)
    )

    result = parse_gaussian_text(log_text, "td-only.log")

    assert result.conformers[0].geometry.atoms[0].x == pytest.approx(1.25)
    assert result.warnings[0].code is ParseWarningCode.TD_INPUT_GEOMETRY_FALLBACK


def test_multiple_ground_state_optimizations_are_all_reported() -> None:
    first_job = (
        " Entering Link 1 = /gaussian/l1.exe\n"
        " -----\n"
        " # b3lyp/6-31g(d) opt\n"
        " -----\n"
        " Charge = 0 Multiplicity = 1\n"
        + _orientation_block("Standard", 1.0)
        + " SCF Done: E(RB3LYP) = -10.0 A.U.\n"
        " Optimization completed.\n"
        " Normal termination of Gaussian 16\n"
    )
    second_job = (
        " Initial command:\n"
        " command-two\n"
        " Entering Link 1 = /gaussian/l1.exe\n"
        " -----\n"
        " # b3lyp/6-31g(d) opt\n"
        " -----\n"
        " Charge = 0 Multiplicity = 1\n"
        + _orientation_block("Standard", 2.0)
        + " SCF Done: E(RB3LYP) = -11.0 A.U.\n"
        " Optimization completed.\n"
        " Normal termination of Gaussian 16\n"
    )

    result = parse_gaussian_text(first_job + second_job, "multiple.log")

    assert len(result.conformers) == 2
    assert [conformer.source_job_index for conformer in result.conformers] == [1, 2]
    assert result.warnings[0].code is (
        ParseWarningCode.MULTIPLE_GROUND_STATE_OPTIMIZATIONS
    )


def test_segments_preserve_source_order_text_and_line_ranges() -> None:
    log_text = """Entering Gaussian System, Link 0=g16
 Initial command:
 command-one
 Entering Link 1 = /gaussian/l1.exe
 -----
 # b3lyp/6-31g(d) opt
 -----
 Normal termination of Gaussian 16
 Initial command:
 command-two
 Entering Link 1 = /gaussian/l1.exe
 -----
 # b3lyp/6-31g(d)
 -----
 Normal termination of Gaussian 16
"""

    jobs = segment_gaussian_jobs(log_text)

    assert len(jobs) == 2
    assert jobs[0].index == 1
    assert jobs[0].start_line == 1
    assert jobs[0].end_line == 8
    assert "command-one" in jobs[0].text
    assert jobs[1].index == 2
    assert jobs[1].start_line == 9
    assert jobs[1].end_line == 15
    assert "command-two" in jobs[1].text
    assert jobs[1].job_types == frozenset({GaussianJobType.SINGLE_POINT})


@pytest.mark.parametrize(
    "excited_keyword",
    (
        "td",
        "td=(nstates=5,root=2)",
        "tda=(nstates=3)",
        "cis(nstates=4)",
    ),
)
def test_excited_state_routes_are_not_ground_state_candidates(
    excited_keyword: str,
) -> None:
    log_text = f"""Entering Link 1 = /gaussian/l1.exe
 -----
 #p cam-b3lyp/6-31g(d) opt {excited_keyword}
 -----
"""

    (job,) = segment_gaussian_jobs(log_text)

    assert job.is_excited_state
    assert GaussianJobType.OPTIMIZATION in job.job_types
    assert not job.is_ground_state_optimization_candidate


def test_empty_text_contains_no_jobs() -> None:
    assert segment_gaussian_jobs("") == ()


def test_simple_xyz_file_parses_atom_symbols_and_coordinates() -> None:
    result = parse_xyz_text(
        """3
water
O  0.000000  0.000000  0.117300
h  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
""",
        "water.xyz",
    )

    assert result.warnings == ()
    assert len(result.conformers) == 1
    conformer = result.conformers[0]
    assert conformer.name == "water"
    assert conformer.charge is None
    assert conformer.multiplicity is None
    assert conformer.energy_hartree is None
    assert conformer.geometry_source is GeometrySource.XYZ_FILE
    assert [atom.atomic_number for atom in conformer.geometry.atoms] == [8, 1, 1]
    assert [atom.symbol for atom in conformer.geometry.atoms] == ["O", "H", "H"]
    assert conformer.geometry.atoms[1].coordinates == pytest.approx(
        (0.0, 0.7572, -0.4692)
    )


def test_xyz_invalid_symbol_fails_clearly() -> None:
    with pytest.raises(ValueError, match="Invalid element symbol 'Xx'"):
        parse_xyz_text(
            """1
bad
Xx 0.0 0.0 0.0
""",
            "bad.xyz",
        )


def test_xyz_malformed_coordinate_row_fails_clearly() -> None:
    with pytest.raises(ValueError, match="x, y, and z coordinates must be numeric"):
        parse_xyz_text(
            """1
bad
C 0.0 nope 0.0
""",
            "bad.xyz",
        )
