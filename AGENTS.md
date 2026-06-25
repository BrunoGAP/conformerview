# AGENTS.md

This project is a local Python application for analyzing molecular conformers from Gaussian16 output files and simple XYZ geometries.

The user is learning while building the project. Prefer clear explanations, small changes, and explicit scientific assumptions over large automatic rewrites.

## Current Project Goal

Build a beginner-friendly local Python application that:

- loads Gaussian16 `.log` files;
- loads simple single-geometry `.xyz` files;
- extracts optimized ground-state molecular geometries;
- compares conformers with the same molecular connectivity;
- calculates RMSD using all atoms;
- calculates RMSD with hydrogen atoms excluded;
- calculates RMSD only after optimal rigid-body alignment and fitting;
- produces fitted pairwise RMSD matrices for multiple conformers;
- generates downloadable images of overlaid conformers similar to the references in `exemplos/`;
- provides a command-line interface and local browser GUI.

## Important Scientific Rules

- Gaussian `.log` files may contain multiple jobs in the same file.
- A file may contain optimization, frequency, single-point, or excited-state jobs.
- Prefer geometries and properties from successfully optimized ground-state jobs.
- Do not treat excited-state calculations as conformer geometries. In particular, route sections containing `td`, `td=...`, or related TD-DFT settings should be excluded unless the user explicitly asks otherwise.
- Explicit fallback: if a file has no successfully optimized ground-state geometry, the initial input geometry of its first TD optimization may be used under the assumption that it was previously optimized at the ground state. Never use the final TD-optimized geometry for this fallback, never attach its excited-state energy, and always emit a user-visible warning explaining the assumption.
- Prefer a geometry associated with successful convergence:
  - `Optimization completed.`
  - `-- Stationary point found.`
  - normal termination of the relevant Gaussian job segment.
- Do not blindly use the last geometry in the whole file.
- Do not assume two geometries are comparable until their atom counts and atomic-number sequences have been validated.
- Do not silently reorder atoms.
- If atom ordering differs, report the issue and stop comparison unless a future explicit atom-mapping feature exists.
- RMSD with hydrogens excluded must preserve the original heavy-atom order after filtering.
- Alignment must be rigid-body only: translation and rotation are allowed; molecular geometry must not be distorted.
- Do not report or expose unfitted RMSD as an analysis mode. Raw coordinate RMSD may exist only as the private numerical kernel applied after fitting.
- XYZ files provide geometry only. They do not provide charge, multiplicity, energy, or explicit bond data.

## Architecture Guidance

Keep the code modular. Do not place core scientific calculations directly in the interface module.

Suggested modules:

- `parsing`: Gaussian16 log parsing and simple XYZ parsing.
- `elements`: element symbol and atomic-number lookup tables.
- `models`: typed molecular data structures such as atoms, geometries, conformers, and parse results.
- `validation`: atom-count, atomic-composition, and atom-order checks.
- `alignment`: centroid calculation and optimal rigid-body alignment.
- `rmsd`: RMSD calculations after mandatory rigid-body alignment.
- `visualization`: conformer overlay rendering and export.
- `interactive`: py3Dmol-based interactive inspection of aligned conformers.
- `interface`: command-line workflow.
- `web_interface`: local browser GUI.
- `tests`: tests for parsing, validation, alignment, and RMSD.

Use type hints for public functions and data models.

The current user-facing interfaces are the `conformer-analyzer` command and
the `conformer-analyzer-gui` local browser application. They must delegate
parsing, validation, RMSD, alignment, and visualization to their scientific
modules rather than reimplementing those calculations.

## Development Style

- Prefer small milestones.
- Add tests for scientific calculations before relying on results.
- Use simple, readable Python.
- Explain important scientific and programming choices in comments or documentation.
- Keep generated images, caches, virtual environments, and temporary outputs out of version control.
- Preserve user data files in `conformers/` and reference images in `exemplos/` unless the user asks otherwise.

## Known Initial Data

The initial project contains:

- `conformers/COOH/conf1.log`: Gaussian16 output with a TD optimization job and a later TD frequency job; neither is eligible as a ground-state conformer.
- `conformers/COOH/conf2.log`: Gaussian16 output with a ground-state optimization, a TD optimization, and a later TD frequency job.
- `conformers/CN/`: Gaussian logs and matching simple XYZ files for another two-conformer example.
- `exemplos/`: reference images showing overlaid conformers as thick colored molecular sticks on a white background with a bottom legend.

The inspected jobs contain `Standard orientation` blocks with 58 atoms and matching atomic-number sequences. Future code should still validate this instead of assuming it, and should use only `conf2.log` job 1 as a ground-state extraction candidate.
