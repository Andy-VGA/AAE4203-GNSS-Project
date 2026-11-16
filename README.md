# AAE4203 GNSS Positioning Project

This repository contains coursework code, processed GNSS observations, Single Point Positioning (SPP) outputs, RTKLIB comparison results, and benchmark ablation experiments prepared for AAE4203.

## Project contents

- **SPP implementation:** GPS-only weighted least-squares SPP using Gauss–Newton iteration.
- **Independent benchmark solver:** experimental GPS L1 broadcast-ephemeris WLS-SPP implementation.
- **RTKLIB comparison:** Python SPP and RTKLIB solutions matched by epoch, with ENU differences and summary statistics.
- **Benchmark V1–V4:** four solution variants and their epoch-level diagnostics, availability and trajectory-continuity comparisons.
- **Processed data:** epoch, satellite, pseudorange, satellite-position, clock, ionospheric, tropospheric, group-delay and elevation CSVs.
- **Figures:** comparison and benchmark diagnostic plots.

## Repository layout

```text
.
├── spp_baseline.py
├── spp_baseline_corrected.py
├── urbannav_benchmark_solver.py
├── rinex2csv.m
├── data/
│   ├── processed_inputs/
│   ├── spp_outputs/
│   └── comparisons/
├── benchmark_V1-V4/
│   ├── G16_benchmark_V1.csv ... G16_benchmark_V4.csv
│   ├── G16_V1_diagnostics.csv ... G16_V4_diagnostics.csv
│   └── G16_ablation_summary.csv
└── figures/
```

## Running the SPP baseline

The corrected GPS-only script expects the processed input CSV files in `outputs/own_inputs/`:

```bash
python spp_baseline_corrected.py
```

The script writes epoch-level solutions and retains failed epochs as NaN rather than silently dropping them. Check the constants and correction sign conventions against the observation model and lab instructions before reusing the results.

The independent benchmark solver is experimental and requires NumPy. Run its help first to see the supported input arguments:

```bash
python urbannav_benchmark_solver.py --help
```

## Benchmark V1–V4 summary

The current ablation summary reports:

| Version | Method | Solved epochs | Availability |
|---|---|---:|---:|
| V1 | GPS elevation-weighted WLS | 1,437 / 1,534 | 93.677% |
| V2 | GPS + Galileo + BeiDou elevation-weighted WLS | 1,533 / 1,534 | 99.935% |
| V3 | V2 + C/N0 weighting | 1,533 / 1,534 | 99.935% |
| V4 | V3 + residual outlier rejection | 1,533 / 1,534 | 99.935% |

The V1–V4 summary currently has **no instructor ground-truth RMSE values**. Therefore, RMSE is marked unavailable in that ablation table; availability and adjacent-epoch continuity statistics are not substitutes for positioning accuracy.

## Python vs RTKLIB comparison

The saved matched-epoch comparison contains 687 epochs. The recorded ENU difference statistics are:

| Metric | Value |
|---|---:|
| East RMS difference | 33.0682 m |
| North RMS difference | 4.1204 m |
| Up RMS difference | 9.4669 m |
| Horizontal RMS difference | 33.3239 m |
| 3D RMS difference | 34.6425 m |

These are **differences between the Python and RTKLIB solutions**, not errors against surveyed ground truth. The comparison also contains sensitivity results; inspect the CSVs and figures before drawing conclusions.

## Important notes

- This repository is a project record, not a claim of survey-grade accuracy.
- RTKLIB-vs-Python agreement does not establish absolute accuracy.
- The benchmark V1–V4 RMSE requires the instructor reference trajectory; do not report a missing RMSE as zero.
- Verify units, satellite-system selection, correction signs, epoch matching and reference-coordinate provenance before final submission.
- Do not commit credentials, private keys, account tokens or any data you are not permitted to redistribute.

## Main outputs

- `G16_ablation_summary.csv` — V1–V4 availability and continuity summary.
- `AAE4203_Python_vs_RTKLIB_statistics.csv` — Python/RTKLIB ENU difference statistics.
- `AAE4203_Python_vs_RTKLIB_all_epochs.csv` — epoch-level comparison.
- `AAE4203_ThreeWay_Comparison_Statistics.csv` — three-way comparison summary.
