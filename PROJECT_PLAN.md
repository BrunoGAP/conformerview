# Project Plan

This plan divides development into small milestones. Each milestone should leave the project in a working, understandable state.

## Milestone 1: Project Skeleton

Status: Complete.

Goal: Create a clean Python project structure without implementing the full application.

Tasks:

- Create source and test directories.
- Choose simple module names for parsing, models, validation, alignment, RMSD, visualization, and interface.
- Add basic development dependencies.
- Add a minimal test runner configuration.

Validation criteria:

- The project imports cleanly.
- The test command runs, even if only placeholder tests exist.
- Documentation explains where each kind of code should live.

## Milestone 2: Molecular Data Models

Status: Complete.

Goal: Define typed structures for molecular data.

Tasks:

- Create an `Atom` model with atomic number, optional element symbol, and coordinates.
- Create a `Geometry` or `MoleculeGeometry` model containing an ordered atom list.
- Create a `Conformer` model containing a name, source file, geometry, charge, multiplicity, and energy when available.
- Decide how to represent parse warnings.

Validation criteria:

- Models are type hinted.
- Tests can construct small molecules manually.
- Atom order is preserved exactly as provided.

## Milestone 3: Gaussian Job Segmentation

Status: Complete.

Goal: Read Gaussian16 `.log` files as multi-job files.

Tasks:

- Split or index logs by route sections and Gaussian job boundaries.
- Extract route text for each job.
- Detect job type from the route section.
- Mark jobs containing `td`, `td=...`, or related excited-state indicators as not eligible for ground-state conformer extraction.

Validation criteria:

- `conf1.log` is recognized as containing an optimization job and a later frequency job.
- `conf2.log` is recognized as containing multiple jobs.
- TD or excited-state route sections are excluded in tests using synthetic log snippets.

## Milestone 4: Optimized Ground-State Geometry Extraction

Status: Complete.

Goal: Extract only successful optimized ground-state geometries.

Tasks:

- Parse `Standard orientation` blocks.
- Parse `Input orientation` blocks if present.
- Parse charge and multiplicity.
- Parse SCF electronic energies.
- Detect successful optimization convergence.
- Select the correct optimized ground-state geometry from a job segment.
- If none exists, use the first TD optimization's input geometry with an explicit warning and no associated energy.

Validation criteria:

- `conf1.log` uses the TD optimization input geometry and reports the ground-state-assumption warning.
- The first optimization in `conf2.log` is eligible for ground-state geometry extraction.
- `conf2.log` does not blindly use the last geometry in the whole file.
- Extracted geometries have 58 atoms for the inspected files.
- Charge and multiplicity are extracted as `0` and `1` for the inspected files.
- Parse warnings are reported when multiple eligible optimization jobs exist.

## Milestone 5: Geometry Comparability Validation

Status: Complete.

Goal: Prevent invalid RMSD comparisons.

Tasks:

- Validate equal atom counts.
- Validate identical atomic-number sequences.
- Provide clear error messages when structures are not comparable.
- Do not reorder atoms.

Validation criteria:

- Matching structures pass validation.
- Different atom counts fail validation.
- Different atomic-number order fails validation.
- Heavy-atom filtering preserves relative order.

## Milestone 6: RMSD Numerical Kernel

Status: Complete.

Goal: Implement the coordinate RMSD metric used after fitting.

Tasks:

- Implement the low-level RMSD formula for two coordinate sets.
- Support all-atom evaluation.
- Support hydrogen-excluded evaluation.
- Add tests with simple known coordinate examples.

Validation criteria:

- Identical geometries have RMSD `0`.
- The private numerical kernel gives the expected value for known coordinates.
- Hydrogen-excluded RMSD ignores atoms with atomic number `1`.
- Unfitted RMSD is not exposed as a project analysis mode.

## Milestone 7: Optimal Rigid-Body Alignment

Status: Complete.

Goal: Align conformers by translation and rotation before RMSD.

Tasks:

- Implement centroid subtraction.
- Implement Kabsch alignment or another well-tested rigid-body method.
- Apply rotation without changing bond lengths or internal geometry.
- Calculate RMSD after alignment.

Validation criteria:

- A molecule compared with a translated and rotated copy has near-zero aligned RMSD.
- The public RMSD calculation always performs alignment first.
- Tests cover all-atom and hydrogen-excluded alignment.

## Milestone 8: Pairwise RMSD Matrices

Status: Complete.

Goal: Compare more than two conformers.

Tasks:

- Accept a list of conformers.
- Validate all conformers against a reference atom sequence.
- Produce pairwise RMSD matrices.
- Support all-atom and hydrogen-excluded modes.
- Fit every conformer pair before calculating RMSD.

Validation criteria:

- Matrix is square.
- Diagonal values are zero.
- Matrix is symmetric.
- Invalid conformer sets fail with a clear message.

## Milestone 9: Basic Visualization Export

Status: Complete.

Goal: Generate static overlay images.

Tasks:

- Convert molecular geometries into a renderable representation.
- Infer bonds conservatively or use an explicit connectivity method.
- Overlay multiple conformers after mandatory rigid-body alignment.
- Use a white background.
- Render each conformer in one solid color.
- Add a bottom legend with conformer labels and percentages when provided.
- Export downloadable PNG images.
- Provide an interactive py3Dmol preview for rotation and view selection before
  final PNG export.

Validation criteria:

- Two aligned conformers can be rendered as red and blue overlays.
- Multiple conformers can use distinct colors.
- Output image resembles the visual style of `exemplos/`.
- Missing or uncertain connectivity is reported instead of hidden.
- The aligned conformers can be rotated, zoomed, hidden, and shown in a local
  interactive viewer.

## Milestone 10: Simple Interface

Status: Complete.

Goal: Provide a beginner-friendly way to use the tool locally.

Tasks:

- Provide a command-line interface first; consider a graphical interface later.
- Accept two or more Gaussian log paths in comparison order.
- Show parsed conformers and scientific warnings.
- Print and save fitted all-atom and heavy-atom RMSD matrices.
- Generate and open the interactive py3Dmol overlay.
- Allow all-atom or heavy-atom fitting for the viewer and provide an interactive
  hydrogen visibility toggle.

Validation criteria:

- Interface code calls parsing, validation, RMSD, alignment, and visualization modules.
- Core scientific calculations are not implemented directly inside the interface module.
- The user can run the workflow locally on the example logs.
- Automated tests cover output files, current scientific values, warnings, and
  browser-launch behavior.

## Milestone 11: Local Graphical Interface

Status: Complete.

Goal: Provide a beginner-friendly browser interface for the validated Gaussian
workflow before expanding the set of input formats.

Tasks:

- Run a local-only web application from a dedicated command.
- Accept two or more Gaussian `.log` uploads in comparison order.
- Let the user choose all-atom or heavy-atom fitting for the viewer and whether
  hydrogens are initially visible.
- Delegate parsing, validation, alignment, RMSD, and visualization to the
  existing scientific modules.
- Display extracted conformers, scientific warnings, and both fitted RMSD
  matrices on the results page.
- Provide downloads for the text report and interactive overlay.
- Open the interactive overlay in a separate browser tab after a successful
  analysis.

Validation criteria:

- Invalid or fewer than two uploads produce a clear page-level error.
- The included Gaussian examples produce the same tested RMSD values as the
  command-line workflow.
- Parser and connectivity warnings remain visible in the page interface.
- The generated py3Dmol viewer opens separately and contains all loaded
  conformers overlaid.
- Automated tests cover upload validation, displayed results, generated files,
  and viewer-launch behavior.

## Milestone 11.1: Per-Conformer Viewer Colors

Status: Complete.

Goal: Let users distinguish conformers with colors chosen directly in the
interactive visualization.

Tasks:

- Add an individual continuous-spectrum color picker for every conformer.
- Apply color changes immediately without resetting the camera or visibility.
- Keep each conformer's swatch and optional hydrogen sticks synchronized with
  its selected color.

Validation criteria:

- Every loaded conformer has its own color picker.
- Changing one picker recolors only its corresponding conformer.
- Color selection continues to work whether hydrogens are shown or hidden.
- Automated tests verify that the generated viewer includes the controls and
  live recoloring behavior.

## Milestone 12: Optional XYZ Support

Status: Complete.

Goal: Add simple `.xyz` geometry input after Gaussian parsing is reliable.

Tasks:

- Parse atom symbols and coordinates from `.xyz` files.
- Convert symbols to atomic numbers.
- Represent XYZ inputs using the same geometry model.
- Validate comparability with Gaussian-derived geometries.

Validation criteria:

- Simple XYZ files parse correctly.
- Invalid symbols or malformed coordinate rows produce clear errors.
- XYZ geometries can enter RMSD calculations only after validation.

## Milestone 13: Scientific Review and Documentation

Status: Complete.

Goal: Make the tool easier to trust and easier to learn from.

Tasks:

- Document how Gaussian jobs are selected.
- Document RMSD formulas.
- Document alignment assumptions.
- Document limitations, especially atom ordering and connectivity.
- Add example outputs.

Validation criteria:

- A beginner can understand the workflow from the documentation.
- Scientific assumptions are explicit.
- Known limitations are listed honestly.
