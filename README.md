# DarkMatter

DarkMatter is a personal MATLAB project for exploring dark-matter/electron scattering calculations in silicon and germanium detector scenarios. It combines numerical rate calculations, background estimates, Poisson simulations and likelihood-ratio calculations, with scripts for particle-mass and velocity-parameter sweeps.

The intended audience is the author and other developers interested in reproducible numerical simulation software. The original implementation dates from 2019; the present goal is to revive it as a portable, testable tool. The code performs ordinary scientific calculations and contains no machine-learning training or inference pipeline.

## Current status

This is a legacy numerical prototype. Three original commits from November 2019 are retained in the repository history. There is no production release or automated regression suite yet, and the archived numerical results have not been independently validated.

The preparation update adds documentation, a license, a small example and a proposed development plan. The existing numerical routines and archived MAT files are unchanged. The legacy driver scripts require maintenance before they can be treated as reproducible end-to-end runs; see [the roadmap](ROADMAP.md).

## Requirements

- MATLAB. The source is legacy MATLAB code; a supported release range has not yet been established.
- Statistics and Machine Learning Toolbox for `poissrnd` in `q0Dist.m` and `q0DistComb.m`. This dependency supplies random sampling; the project does not train or run AI models.
- A valid license for MATLAB and any toolbox being used.

The small example below uses only base MATLAB functions. GNU Octave compatibility has not been verified. A portable Python numerical implementation is proposed in the roadmap and is not implemented in this version.

## Start with the small example

Clone or download the repository, start MATLAB and set its current folder to the repository root. Run:

```matlab
run('examples/minimal_example.m')
```

The example calls the velocity helper at one non-boundary point and constructs a 500-bin silicon background. It loads no archived workspaces, performs no long simulation, and writes no output files. It should print the velocity-integral value and a completion message.

From a terminal with a licensed MATLAB installation:

```text
matlab -batch "run('examples/minimal_example.m')"
```

Execution of this command was attempted during preparation on 4 November 2026, but the local MATLAB R2022b license had expired before the example could run. The example has been inspected in source; it is not presented as a passed runtime test.

## Code map

| File | Purpose |
| --- | --- |
| `Eta.m` | Piecewise velocity-integral helper. |
| `dRdE_Si.m`, `dRdE_Ge.m` | Numerical scattering-rate routines. |
| `background.m`, `backgroundGe_asInPaper.m`, `getBackgroundGe.m` | Background construction and aggregation. |
| `q0Dist.m`, `q0DistComb.m` | Poisson sampling, likelihood optimization and likelihood-ratio distributions. |
| `ScatteringRateIter.m`, `ScatteringCCD.m` | Legacy simulation drivers and exploratory plotting. |
| `changingspeed.m`, `speedplots.m` | Velocity-parameter sweeps and archived-result plots. |
| `plasterTogetherDM.m` | Combines selected archived result arrays. |
| `examples/minimal_example.m` | Small example independent of saved MAT workspaces. |

Some legacy routines document energies and momenta in eV, velocities in cm/s, detector target masses in kg and exposure in kg-years. Consistent units, array contracts and parameter validation are work items, rather than a verified package-wide contract.

## Archived inputs and results

`DM-FF.mat` contains tabulated numerical arrays used by the rate calculations. `Background.mat` contains numerical background arrays. The other MAT files contain historical outputs or complete saved simulation workspaces. Their contents and limitations are described in [DATA.md](DATA.md).

The old saved workspaces include serialized function handles with local computer-folder paths. These metadata fields are retained in the original files and commit history. Future runs should save explicit numerical outputs and configuration metadata into a local `results/` folder.

## Planned work

[ROADMAP.md](ROADMAP.md) proposes approximately **28-36 hours** of software development: a portable numerical core, numerical and statistical regression tests, reproducible batch runs and clearer result visualization. These are estimates for proposed future work. The implementation is not included in this preparation update.

## License

The original software and documentation are licensed under the [MIT License](LICENSE). MATLAB and its toolboxes are separately licensed products. Any externally sourced numerical tables retain the rights of their original source; the legacy repository does not identify an upstream source for `DM-FF.mat`, so that provenance should be recorded before redistribution.
