# G16 UrbanNav GNSS SPP V1–V4 Ablation (experimental)

This package contains four independently generated Python SPP trajectories from `urbannav_deep1.obs` and `urbannav_deep1.nav`. No instructor or public benchmark ground truth was used. The four variants are **not** official-scored results.

| Version | Incremental change |
|---|---|
| V1 | GPS L1 C1C pseudorange; broadcast ephemerides; ionosphere/troposphere/earth rotation corrections; 15° elevation mask; elevation-weighted least squares |
| V2 | V1 plus Galileo E1 C1C and BeiDou B1I C2I (excluding BDS GEO); separate receiver clock offsets per constellation |
| V3 | V2 plus signal-strength C/N0 weights, `min(1,max(0.3,CN0/40))²` multiplied by elevation weights |
| V4 | V3 plus post-fit residual screening: remove largest residual if above 60 m and redundant measurements remain; at most two removals per epoch |

## Run

Requires Python 3.9+ and NumPy (pandas/matplotlib only used for postprocessing figures).

Place `G16_benchmark_ablation.py`, `urbannav_deep1.obs` and `urbannav_deep1.nav` in the same folder and run:

```bash
python G16_benchmark_ablation.py --obs urbannav_deep1.obs --nav urbannav_deep1.nav --outdir G16_ablation
```

The outputs are `G16_benchmark_V1.csv` ... `G16_benchmark_V4.csv`, per-version diagnostics, and `G16_ablation_summary.csv`. Each trajectory contains only solved epochs and has the five required CSV columns.

## Interpretation and limitations

- **Availability** is calculated locally as solved epochs / 1534 observed epochs. The official system matches GPS timestamps to its reference epochs; its reported availability can differ.
- `adjacent_jump_p95_m` is the 95th percentile of consecutive ECEF displacements for successive solved epochs within 1.1 seconds. This is an **internal trajectory continuity statistic, not 2D position error P95 or accuracy**. UrbanNav is a moving receiver, so true motion contributes to this statistic.
- **Official 2D RMSE, 3D RMSE, 2D error P95 and rank are unknown** without the instructor's private reference and evaluation feedback. Do not put the internal jump statistic in the official P95 error column.
- Galileo and BeiDou use broadcast ephemerides, separate inter-system clock offsets, GPS broadcast Klobuchar ionosphere as an approximation, and a simple standard-atmosphere troposphere model. These corrections and group-delay choices have not been independently verified against RTKLIB for this multi-GNSS dataset.
- BeiDou GEO satellites are excluded because their ECEF orbit model requires additional handling. GLONASS is excluded. A 15° elevation mask is fixed across all variants.
- The provided **previously submitted** `G16_benchmark.csv` used GPS-only **elevation + C/N0** weighting, whereas experimental V1 uses elevation-only weighting to make V3 a clean incremental comparison. The two GPS-only trajectories are not identical.
- These experiments were executed on the benchmark observations; the tutorial also asks for an example-dataset regression check before selecting a final submission. That independent example-data validation has **not** been completed here.
- V4 rejects outliers based on post-fit residuals only; no ground-truth-driven tuning was performed. No variant should be claimed to be more accurate until official feedback is available.

## Submission caution

Do **not** upload all four files as `G16_benchmark.csv`. Only a deliberate selected variant, renamed to the exact filename and checked with the course `tools/check_submission.py`, should replace the current submitted file. The teaching platform uses the latest recognized submission, even if it is worse.
