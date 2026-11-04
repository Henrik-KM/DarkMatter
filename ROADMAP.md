# DarkMatter revival roadmap

This is a proposed software-development backlog for reviving the original personal numerical project. Estimates total **28-36 hours** and describe future work, not completed implementation or elapsed time.

## 1. Portable numerical core and explicit inputs: 10-12 hours

- Extract a documented configuration for particle mass, detector target, energy/momentum grids, exposure and velocity parameters.
- Define and validate array shapes and units at the numerical boundaries.
- Implement the core rate, velocity and background calculations in a small Python/NumPy/SciPy package, retaining the MATLAB source as the legacy reference.
- Resolve the outdated `Eta` call in the germanium path: `Eta.m` now expects velocity parameters, while `dRdE_Ge.m` calls the earlier three-argument interface.
- Preserve a documented numerical baseline and explain any intentional behavior changes.

**Completion evidence:** a small rate calculation runs from explicit inputs in a fresh process, with no saved workspace or MATLAB license required. The input contract and comparisons with legacy numerical fixtures are documented.

## 2. Numerical boundaries and likelihood regression tests: 8-10 hours

- Cover the exact `v_min = v_esc - v_E` and `v_min = v_esc + v_E` boundaries in `Eta.m`, where strict comparisons currently leave an output unassigned.
- Exercise zero and near-zero expected counts, zero signal rates and finite log-likelihood behavior.
- Replace poorly scaled optimization bounds/tolerances with a documented, testable parameterization.
- Introduce explicit random seeds and small Poisson fixtures for both single-detector and combined-detector calculations.
- Separate regression agreement with legacy results from scientific validation of the underlying physical assumptions.

**Completion evidence:** focused tests reproduce seeded examples and cover invalid, limiting and boundary inputs with documented numerical tolerances.

## 3. Reproducible parameter sweeps and saved outputs: 6-8 hours

- Replace manual script edits and hidden workspace state with a batch runner accepting explicit configuration.
- Record grids, units, seed, version and run status beside numerical outputs.
- Add checkpoint/resume and small-run limits so expensive sweeps can be interrupted and inspected.
- Reconcile the 50-point velocity grids in the sweep driver with the 51-point grids constructed by `speedplots.m`.
- Export only required numerical results and metadata, excluding serialized function handles and local paths.

**Completion evidence:** a reduced sweep completes from a clean checkout, resumes from a checkpoint and recreates a plot from its own saved configuration.

## 4. Result review, data provenance and examples: 4-6 hours

- Provide plot reconstruction and comparisons without assuming unrelated variables exist in the MATLAB workspace.
- Document the origin, attribution and redistribution terms for numerical input tables.
- Add small synthetic fixtures and portable examples for development.
- Expand installation, supported-environment and interpretation notes, including uncertainty and limitations.

**Completion evidence:** a reviewer can install the portable tool, run its small example and inspect a labeled result without using private files or manually restoring a historical workspace.

## Preparation findings

- The current numerical scripts are retained as a legacy baseline.
- Source review identified missing function arguments, unhandled equality boundaries, inconsistent sweep/plot grids and references to absent saved files.
- The archived MAT workspaces contain local path metadata inside serialized function handles. Publication preparation must also consider the original Git history.
- MATLAB R2022b is installed locally, but execution during preparation was blocked by an expired license. No end-to-end simulation or numerical correctness result is claimed.
