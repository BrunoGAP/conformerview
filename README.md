# Conformer Analyzer

Conformer Analyzer is a beginner-friendly local Python application for
comparing molecular conformers from Gaussian16 output files and simple XYZ
geometries.

It loads Gaussian `.log` files or simple `.xyz` files, extracts geometries,
validates that conformers are comparable, calculates fitted RMSD values, and
generates interactive aligned overlays similar to the examples in `exemplos/`.

## Project Status

All planned milestones are complete. The importable Python package, module
boundaries, test setup, typed molecular data models, Gaussian job segmentation,
geometry extraction, comparability validation, optimal rigid-body fitting, and
pairwise RMSD matrices are in place for all-atom and heavy-atom modes. Static
PNG export renders aligned conformer overlays in the target visual style, and
a py3Dmol page provides interactive rotation before final image selection. A
command-line and browser workflows connect Gaussian or XYZ parsing, validation,
both fitted RMSD matrices, report export, and automatic viewer launch. Every
reported RMSD is calculated after alignment.

## Workflow at a Glance

1. Choose two or more conformer files in the order you want them compared.
2. Parse Gaussian logs or XYZ geometries into the shared molecular data model.
3. Validate that every conformer has the same atom count and atomic-number
   sequence as the first conformer.
4. Align every pair by optimal rigid-body translation and rotation.
5. Report all-atom and heavy-atom fitted RMSD matrices.
6. Open an interactive 3D overlay with per-conformer colors and hydrogen
   visibility controls.

The tool intentionally stops before RMSD if atom ordering or composition does
not match. It does not guess atom mappings.

## Developer Setup

The project uses Python 3.11 or newer. From the project directory, create and
activate a virtual environment, then install the package and development tools:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --editable ".[dev]"
```

Run the tests with:

```powershell
python -m pytest
```

The editable installation means changes under `src/` are available immediately
without reinstalling the package.

## Project Structure

```text
conformer-analyzer/
|-- conformers/
|   `-- COOH/                     Example Gaussian16 log files
|-- exemplos/                      Reference overlay images
|-- src/
|   `-- conformer_analyzer/
|       |-- parsing.py             File reading and Gaussian/XYZ parsing
|       |-- elements.py            Element symbol and atomic-number lookup
|       |-- models.py              Typed molecular data structures
|       |-- validation.py          Atom count and order checks
|       |-- alignment.py           Rigid-body translation and rotation
|       |-- rmsd.py                RMSD calculations
|       |-- visualization.py       Overlay rendering and image export
|       |-- interactive.py         Interactive py3Dmol viewer export
|       |-- interface.py           User-facing command-line workflow
|       `-- web_interface.py       Local browser GUI
|-- tests/                         Automated scientific and behavior tests
|-- pyproject.toml                 Package, dependency, and pytest configuration
`-- PROJECT_PLAN.md                Milestone roadmap
```

Core scientific calculations belong in their dedicated modules rather than in
`interface.py`. Tests should mirror these responsibilities as functionality is
added.

## Molecular Data Models

The models in `models.py` use immutable dataclasses so parsed molecular data
cannot be reordered or edited accidentally during comparison:

- `Atom` stores an atomic number, optional symbol, and Cartesian coordinates in
  angstroms;
- `MoleculeGeometry` stores atoms as an ordered tuple, where tuple position
  defines atom correspondence between conformers;
- `Conformer` combines a geometry with its source path and optional Gaussian
  charge, multiplicity, and electronic energy in Hartree;
- `ParseWarning` records a non-fatal issue and can identify its Gaussian job and
  source line;
- `ParseResult` groups the conformers and warnings extracted from one file.

Atomic number is the authoritative atom identity. The optional element symbol
is retained for display and input convenience, but it must never be used to
silently reorder a geometry.

## Scientific Assumptions

Conformer Analyzer makes these assumptions explicit:

- Coordinates are Cartesian coordinates in angstroms.
- RMSD is meaningful only for conformers with the same atom count and the same
  atomic-number sequence.
- Atom order defines correspondence. The software does not reorder atoms,
  detect symmetry-equivalent atoms, or infer atom mappings.
- Gaussian conformers are selected from successfully optimized ground-state
  jobs when available.
- Excited-state and TD-DFT jobs are excluded from ground-state conformer
  extraction. The only exception is the documented TD-input fallback, which uses
  the first TD optimization input geometry and warns the user.
- RMSD is always reported after optimal rigid-body fitting. Unfitted RMSD is not
  a user-facing analysis mode.
- Heavy-atom mode removes hydrogen atoms from both fitting and RMSD, preserving
  the original order of all remaining atoms.
- Viewer and image connectivity is inferred from covalent radii, not read from
  Gaussian or XYZ bond records.

## Current Capabilities

The application can:

- load Gaussian16 `.log` files;
- load simple single-geometry `.xyz` files;
- extract the final converged optimized ground-state geometry;
- avoid excited-state geometries, including TD-DFT jobs, unless the documented
  TD-input fallback is required;
- compare conformers only when they contain the same atoms in the same order;
- calculate RMSD using all atoms;
- calculate RMSD excluding hydrogen atoms;
- calculate RMSD only after optimal rigid-body alignment and fitting;
- generate pairwise RMSD matrices for multiple conformers;
- export images of overlaid conformers;
- provide a beginner-friendly command-line interface and local browser GUI.

## Why Validation Matters

RMSD is meaningful only when each coordinate in one structure corresponds to the same atom in the other structure.

For this project, the first implementation will require:

- same number of atoms;
- same atomic numbers;
- same atom order.

The program must not silently reorder atoms. Atom mapping may be considered later, but it is a separate scientific feature and should be implemented explicitly.

The validation module enforces this rule before RMSD calculations:

- atom counts must match;
- atomic numbers must match at every tuple position;
- a mismatch reports the first differing one-based atom index and stops;
- no atom mapping or reordering is attempted;
- hydrogen exclusion uses stable filtering, preserving the relative order and
  coordinates of all remaining atoms.

The included `conf1.log` fallback geometry and `conf2.log` optimized geometry
both contain 58 atoms with matching atomic-number sequences. Removing hydrogen
from each leaves the same ordered sequence of 37 heavy atoms.

## RMSD Metric

After fitting, RMSD is evaluated over corresponding Cartesian coordinates:

```text
RMSD = sqrt((1 / N) * sum((dx)^2 + (dy)^2 + (dz)^2))
```

The calculation always validates atom correspondence and performs optimal
rigid-body fitting first. It supports either all atoms or stable removal of
atoms with atomic number `1`. The raw coordinate formula remains a private
numerical kernel; unfitted RMSD is not exposed or reported as an analysis mode.
An empty coordinate set, unequal coordinate count, invalid atom order, or
hydrogen-only geometry in heavy-atom mode produces a clear error.

## Optimal Rigid-Body Alignment

Alignment uses the Kabsch algorithm implemented with NumPy linear algebra:

1. Subtract the centroid from each coordinate set.
2. Use singular value decomposition to find the optimal rotation.
3. Reject reflection by enforcing a positive rotation determinant.
4. Apply one translation and rotation to the entire candidate geometry.

The transformation does not reorder atoms or alter internal geometry. Tests
verify that pairwise distances are preserved to numerical precision and that a
translated and rotated copy has near-zero aligned RMSD. A mirror reflection is
not accepted as a rotation because it would invert molecular chirality.

Heavy-atom mode fits the transform using only non-hydrogen atoms, then applies
that same transform to every atom. The aligned RMSD is subsequently calculated
over the heavy atoms only.

For the included `conf1.log` and `conf2.log` geometries:

| Mode | Fitted RMSD (angstrom) |
| --- | ---: |
| All atoms | `2.1199556775` |
| Heavy atoms | `1.5507099839` |

For this example, the largest change in any candidate atom-pair distance caused
by alignment was approximately `1.1e-14` angstrom, consistent with floating-point
roundoff rather than molecular distortion.

## Pairwise RMSD Matrices

Pairwise matrix functions accept any sequence of parsed conformers or raw
geometries and validate every geometry against the first geometry's atom
sequence before calculating values. Matrices are square and symmetric, with
exact zeroes on the diagonal. Every pair is rigid-body fitted before RMSD is
calculated, with support for all atoms or hydrogen exclusion.

The current inputs are `conformers/COOH/conf1.log` and
`conformers/COOH/conf2.log`. Their aligned matrices, in angstroms, are:

```text
All atoms                 Heavy atoms
        conf1       conf2         conf1       conf2
conf1  0.00000000  2.11995568    0.00000000  1.55070998
conf2  2.11995568  0.00000000    1.55070998  0.00000000
```

## Example Outputs

The repository includes example output reports under `outputs/`.

For the COOH Gaussian example:

```text
outputs/COOH/rmsd_matrix.txt
outputs/COOH/interactive_overlay.html
outputs/COOH/cooh_aligned_overlay.png
```

The RMSD report contains this result excerpt:

```text
All atoms (fit and RMSD use all atoms)
                     C1             C2
C1         0.0000000000   2.1199556775
C2         2.1199556775   0.0000000000

Heavy atoms (fit and RMSD exclude hydrogen)
                     C1             C2
C1         0.0000000000   1.5507099839
C2         1.5507099839   0.0000000000
```

For the CN Gaussian example:

```text
outputs/CN/rmsd_matrix.txt
outputs/CN/interactive_overlay.html
```

The CN report currently gives:

| Example | All-atom RMSD | Heavy-atom RMSD |
| --- | ---: | ---: |
| `conformers/CN/conf1.log` vs `conformers/CN/conf2.log` | `2.2736624911` | `1.7682491206` |

The included `conformers/CN/conf1.xyz` and `conformers/CN/conf2.xyz` files are
simple XYZ geometries with matching atom order. They are useful for testing the
XYZ path or for mixing one `.log` and one `.xyz` input in the GUI.

## Aligned Conformer Overlay

The visualization module rigidly aligns every conformer onto the first one and
exports a PNG with a white background, thick rounded molecular sticks, one
solid color per conformer, and a bottom legend. Hydrogens are hidden by default
to match the skeletal PyMOL examples, but they can be displayed explicitly.
Labels and conformer percentages are optional.

Connectivity is inferred conservatively from interatomic distances and
single-bond covalent radii. Every export returns a warning that inferred bonds
should be verified when explicit connectivity becomes available. Unsupported
elements and visible atoms without inferred bonds are reported clearly rather
than silently omitted.

The current red/blue aligned overlay is exported to
`outputs/cooh_aligned_overlay.png`.

### Interactive 3D preview

Generate the current interactive viewer from the project directory with:

```powershell
python -m conformer_analyzer.interactive
```

This writes `outputs/cooh_interactive_overlay.html`. Open that file in a web
browser to rotate, zoom, and translate the aligned structures. The toolbar can
show or hide each conformer, choose its color from a continuous color palette,
toggle hydrogens, reset the camera, and save the current view as a PNG. Color
changes are applied immediately to the selected conformer and its legend
swatch. The molecular coordinates and explicit inferred bonds are embedded in
the page; py3Dmol loads its matching 3Dmol.js library from the internet when
the viewer opens.

Alternative `.log` or `.xyz` inputs can be passed as positional arguments, and
the output path can be changed with `--output`.

## Command-Line Workflow

After activating the project virtual environment, analyze two conformer files
with one command:

```powershell
conformer-analyzer conformers/COOH/conf1.log conformers/COOH/conf2.log
```

The command validates atom correspondence, prints and writes the fitted
all-atom and heavy-atom RMSD matrices, reports parser and connectivity warnings,
creates `outputs/interactive_overlay.html`, and opens it in the default browser.
The RMSD report is written to `outputs/rmsd_matrix.txt`.

Inputs may be optimized Gaussian `.log` files, simple `.xyz` files, or a mix of
both. XYZ files provide geometry only, so charge, multiplicity, and energy are
reported as unknown.

Viewer alignment uses all atoms by default. To fit the displayed conformers
using heavy atoms only:

```powershell
conformer-analyzer conformers/COOH/conf1.log conformers/COOH/conf2.log --viewer-fit-atoms heavy
```

Hydrogens are hidden initially and can always be toggled in the viewer. Use
`--show-hydrogens` to show them initially, `--output-dir` to choose another
destination, or `--no-open` to generate files without opening a browser tab.

## Graphical Interface

Start the local browser interface after activating the project environment:

```powershell
conformer-analyzer-gui
```

The command opens `http://127.0.0.1:5000/` and keeps all processing on the
local computer. Select two or more `.log` or `.xyz` files, arrange them in the
desired C1, C2, ... comparison order, choose the viewer fitting options, and
run the analysis. The page displays:

- every extracted conformer and its available metadata;
- fitted all-atom and heavy-atom RMSD matrices;
- Gaussian selection and inferred-connectivity warnings;
- links to download the text report and reopen the interactive overlay.

After a successful run, the aligned py3Dmol visualization opens in a separate
browser tab. If the browser blocks or cannot open that tab, use the **Open 3D
overlay** button on the results page. Uploaded files and generated results are
stored under `outputs/gui_runs/`, which is excluded from version control.
This location is anchored to the project directory even when the GUI launcher
is started from `.venv/Scripts` or another working directory. The launcher also
refuses to start a second server on the same port, preventing an older running
instance from serving stale viewer code.

Use `conformer-analyzer-gui --no-open` if you do not want the home page opened
automatically. The server binds to `127.0.0.1` by default and is intended for
local use, not deployment on a public network.

## Gaussian Log Parsing Notes

Gaussian `.log` files can contain more than one job. For example, a single file may contain:

- an optimization job;
- a frequency job;
- a single-point calculation;
- an excited-state calculation using `td` or TD-DFT settings.

The application should extract conformer geometries and properties only from optimized ground-state jobs. A safe parser should segment the file by Gaussian job or route section, inspect each route, and select successful non-TD optimization jobs.

Useful Gaussian markers include:

- route section lines beginning with `#`;
- `Charge = ... Multiplicity = ...`;
- `Standard orientation`;
- `Input orientation`;
- `SCF Done`;
- `Optimization completed.`;
- `-- Stationary point found.`;
- `Normal termination of Gaussian`;
- `Error termination`.

The current parser segments appended jobs using their `Entering Link 1`
boundaries, extracts and unwraps each printed route section, and classifies
optimization, frequency, and single-point jobs. A route containing `td`, `tda`,
`cis`, or another explicitly supported excited-state keyword is marked as
excited-state and cannot be a ground-state optimization candidate. Convergence
is intentionally not decided from the route; that belongs to geometry
extraction.

Inspection of the included files established that:

- `conf1.log` contains a TD optimization followed by a TD frequency job, so it
  has no ground-state optimization candidate;
- `conf2.log` contains a ground-state optimization, a TD optimization, and a TD
  frequency job. Only its first job is a ground-state optimization candidate.

These classifications depend on complete route continuations. For example,
Gaussian wraps `freq` as `fre` plus `q` in these logs, which the parser rejoins
before detecting the job type.

## Geometry Selection

Geometry extraction follows this order:

1. Use geometries from non-excited optimization jobs that report successful
   convergence and normal Gaussian termination.
2. If no such geometry exists, use the first printed orientation from the first
   TD optimization job. This is the submitted geometry in Gaussian's input or
   standard orientation frame, not the final excited-state optimized geometry.
3. If neither source exists, return no conformer and report a warning.

The TD-input fallback records `TD_OPTIMIZATION_INPUT` as its geometry source,
leaves electronic energy unset, and returns a warning stating that no explicit
optimized ground-state geometry was found. The warning also states the required
assumption: the TD input geometry was previously optimized at the ground state.
The graphical interface displays this parser warning prominently above the
RMSD matrices.

## XYZ Parsing Notes

XYZ support is intentionally simple. The parser reads one geometry from a file
whose first line is an atom count, second line is a comment, and following rows
contain an element symbol plus x, y, and z coordinates in angstroms. Element
symbols are converted to atomic numbers and stored in the same `Atom` and
`MoleculeGeometry` models used by Gaussian parsing.

Malformed atom counts, unknown element symbols, missing coordinate columns, and
non-numeric coordinates fail with explicit errors. XYZ inputs can enter RMSD
calculations only after the normal comparability validation confirms the same
atom count, atomic numbers, and atom order as the reference conformer.

## Known Limitations

- Atom mapping is not implemented. If two files describe the same molecule with
  a different atom order, the comparison stops.
- Symmetry-equivalent atom matching is not implemented.
- Connectivity is inferred from covalent radii for visualization. It should be
  checked when exact bond data matters.
- Static PNG rendering and interactive viewing support the elements configured
  in the covalent-radii and element-symbol tables. Unsupported visualization
  elements fail explicitly.
- XYZ parsing reads a single geometry and does not parse charge, multiplicity,
  energy, bonds, trajectories, or multi-frame XYZ files.
- Gaussian parsing is intentionally conservative. It supports the markers used
  by the inspected Gaussian16 outputs and tested synthetic cases; unusual output
  formatting may need new parser tests.
- The interactive viewer embeds molecular data locally, but loads 3Dmol.js from
  the internet when the HTML page opens.

## Reference Visual Style

The images in `exemplos/` show the visual target:

- white background;
- overlaid conformers;
- thick rounded stick rendering;
- each conformer drawn in one solid color;
- no atom labels;
- bottom legend with colored line swatches;
- labels such as `C1`, `C2`, etc.;
- percentages shown beside conformer labels;
- downloadable static image output.

The visible evaluation watermark in the references appears to come from external software and should not be reproduced intentionally.

## Architecture

The code should be modular:

- parsing code should read files and extract molecular data;
- data models should describe atoms, geometries, conformers, and parse results;
- validation code should decide whether structures are comparable;
- alignment code should perform rigid-body fitting;
- RMSD code should calculate scientific values;
- visualization code should render/export images;
- interface code should call these modules without containing scientific calculations itself.

Type hints and tests should be used throughout the project.

## Development Approach

This project should grow through small milestones. Each milestone should include validation criteria so that scientific correctness is checked before more features are added.

See `PROJECT_PLAN.md` for the proposed roadmap.
