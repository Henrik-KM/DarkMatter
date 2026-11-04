# Numerical inputs and historical workspaces

The repository contains 14 MATLAB MAT files inherited from the original 2019 project. The preparation review matched all 31 downloaded files to the Git blobs at commit `dd2a98534646695c5eaa34350bdecaa95a3bbb91`. Every MAT-file blob reachable from the three original commits is also present in that current snapshot.

## Inventory

| Files | Contents observed |
| --- | --- |
| `Background.mat` | Two numerical arrays. |
| `DM-FF.mat` | Four numerical arrays, including the tabulated arrays used by the rate drivers. |
| `FDM_1_v0.mat`, `FDM_1_vE.mat`, `FDM_1_vesc.mat` | One numerical result array per file. |
| `FDM_q2_v0.mat`, `FDM_q2_vE.mat`, `FDM_q2_vesc.mat` | One numerical result array per file. |
| `DMSimulationResults.mat` | Saved workspace with 96 top-level variables, including numerical values, flags and function handles. |
| `DMreqn0exposure0.1and1.mat` | Saved workspace with 93 top-level variables. |
| `DMreqn2exposure0.1first99correct.mat`, `DMreqn2exposure0.1restcorrect.mat`, `DMreqn2exposure1first79correct.mat`, `DMreqn2exposure1restcorrect.mat` | Saved workspaces with 94 top-level variables each. |

The input tables and historical output arrays are not independently validated physical reference data. Do not treat their presence as proof that the current scripts reproduce a published experimental result.

## Privacy review

The review examined all 20 unique text blobs across the three original commits and decoded the 14 MAT files with SciPy. No matches for credential patterns or email addresses were found in the reviewed text and decoded fields.

Six full saved workspaces contain local computer-folder paths in serialized function-handle metadata. Raw paths are deliberately omitted here. These paths exist in the historical blobs as well as the current snapshot; changing only the latest files would leave the old metadata in Git history. The repository was made public on 4 November 2026 at the owner's request, with the original files and commit history retained.

This was a targeted review, including nested decoded fields. It does not certify the absence of every possible sensitive value, including opaque internal binary payloads.

## Provenance and attribution

The legacy source does not identify the upstream origin or redistribution terms of the arrays in `DM-FF.mat`. Record that provenance before a public release; the new software license does not claim ownership of third-party data.

The comment in `backgroundGe_asInPaper.m` cites [Search for Low-Mass WIMPs with SuperCDMS, arXiv:1402.7137](https://arxiv.org/abs/1402.7137). The precise relationship between that reference, the parameters in the legacy function and the saved background arrays still needs to be checked.

## Future output format

New runs should save an explicit selection of numerical results, the parameter configuration, units, seed and software version. Prefer small synthetic fixtures for automated tests. Keep raw workspace dumps, local paths, credentials and unrelated files out of version control.
