# Phase 0 audit results

Generated 2026-09-21 on R 4.5.2, scDiagnostics 1.7.3.

Timings are single runs on one machine. They are a baseline to beat, not a benchmark.

## Legacy startup cost

- Bundled data: **36.1 MB** on disk across 3 files
- Read at process start: **1.98 s**
- Resident in memory: **193.3 MB**

| object | dimensions | assays |
|---|---|---|
| query | 23341 x 651 | counts,logcounts |
| reference | 23341 x 1518 | counts,logcounts |
| reference subset | 23341 x 1337 | counts,logcounts |

## Parity: forked `Shiny/support/` vs package

| diagnostic | fork (s) | package (s) | agreement | detail |
|---|---:|---:|---|---|
| `projectPCA` | 4.07 | 0.55 | identical | 10845 shared fields, max abs diff 2.84e-14 |
| `detectAnomaly` | 4.72 | error | ERROR | one side failed to run |
| `calculateWassersteinDistance` | 22.87 | 21.77 | DIFFERENT | 1005 shared fields, max abs diff 12.5 |
| `calculateDiscriminantSpace` | 23.23 | 24.49 | identical | 11690 shared fields, max abs diff 9.95e-14 |
| `calculateGraphIntegration` | 3.78 | 2.16 | DIFFERENT | 127568 shared fields, max abs diff 206 |
| `calculateVarImpOverlap` | 14.94 | 9.23 | identical | 700245 shared fields, max abs diff 0.00e+00 |
| `plotCellTypePCA` | 1.99 | 0.94 | no numeric leaves | nothing comparable (plot-only result) |
| `plotGeneExpressionDimred` | 0.25 | 0.31 | no numeric leaves | nothing comparable (plot-only result) |

Errors and caveats:

- `detectAnomaly` package error: 'newdata' must be a data.frame, matrix, or sparse matrix.
- `detectAnomaly`: fork exposes a fixed anomaly_threshold only; the package adds threshold_method (MAD/absolute), mad_multiplier, n_hvgs and max_cells_* downsampling
- `calculateVarImpOverlap`: fork uses randomForest; the package Imports ranger
- `plotCellTypePCA`: ggplot results; the package downsamples to max_cells_* = 2000 by default, the fork does not
- `plotGeneExpressionDimred`: fork argument is `se_object`, package argument is `sce_object`

## Package functions on package preset data

Reference 392 x 1500, query 392 x 503, 4.2 MB in memory (vs 193.3 MB for the legacy bundle).

| diagnostic | seconds |
|---|---:|
| `detectAnomaly` | 0.35 |
| `calculateWassersteinDistance` | 8.85 |
| `plotCellTypePCA` | 0.36 |
| `calculateGraphIntegration` | 1.08 |
| `calculateGeneShifts` | 1.34 |

