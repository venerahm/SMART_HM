# Requirements

The scripts were developed for R 4.x.

Required CRAN packages:

- `dplyr`
- `ggplot2`
- `MASS`
- `mvtnorm`
- `patchwork`
- `purrr`
- `readr`
- `scales`

The scripts also use the base R `stats` and `utils` packages.

## Validation

Parse all R files:

```bash
Rscript -e 'invisible(lapply(list.files(pattern = "\\\\.R$", recursive = TRUE), parse))'
```

Run a small primary-simulation smoke test:

```bash
SMART_HM_NITER=1 Rscript scripts/run_primary_simulation.R
```
