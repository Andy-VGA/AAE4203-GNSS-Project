# AAE4203 GNSS Positioning Project

This archive consolidates the available project source code, raw UrbanNav RINEX files, processed observation inputs, SPP/RTKLIB outputs, benchmark artifacts and figures.

## Raw GNSS data included
- `data/raw/urbannav_deep1.obs`
- `data/raw/urbannav_deep1.nav`

The earlier filenames `mydata.obs` and `mydata.nav` were not present in the accessible working files during this audit. The included UrbanNav files are the raw RINEX files that were actually available.

## Layout
- `src/`: Python SPP, benchmark solver and MATLAB RINEX conversion source.
- `data/raw/`: original OBS/NAV files available in the working set.
- `data/processed/`: epoch, satellite, pseudorange, satellite position, clock, ionosphere, troposphere, group-delay, elevation and reference CSV inputs.
- `results/`: SPP outputs, RTKLIB POS, comparison tables and V1-V4 benchmark outputs/diagnostics.
- `figures/`: organized diagnostic and comparison plots.
- `FILE_MANIFEST.csv`: file sizes and SHA-256 checksums.

