# Phase 0 findings

What the baseline audit of the previous app turned up, and what was done about
it. Raw numbers are in [phase0-results.md](phase0-results.md), produced by
[phase0_baseline.R](phase0_baseline.R).

## 1. The old app did not compute what the package computes

`Shiny/support/` held 3,798 lines of hand-forked package functions. Run
side by side on identical inputs under identical seeds:

| diagnostic | agreement with the package |
|---|---|
| `projectPCA` | identical (max abs diff 2.8e-14) |
| `calculateDiscriminantSpace` | identical (max abs diff 1.0e-13) |
| `calculateVarImpOverlap` | identical |
| `calculateWassersteinDistance` | **different** — max abs diff 12.5 over 1,005 shared values |
| `calculateGraphIntegration` | **different** — max abs diff 206 over 127,568 shared values |
| `detectAnomaly` | not comparable: the package call failed (see §2) |
| `plotCellTypePCA`, `plotGeneExpressionDimred` | figures only |

So two of the six comparable diagnostics produced materially different numbers
from the published package the app was presenting as its engine. That is the
finding that justified deleting the fork rather than maintaining it.

## 2. A real bug in scDiagnostics, found by the audit

While running the comparison, the package's own `detectAnomaly()` failed on the
old app's demo data:

```
'newdata' must be a data.frame, matrix, or sparse matrix.
```

Isolated to `R/detectAnomaly.R`: the per-cell-type loop subsets the query
projection with `query_mat[which(query_cell_types %in% cell_type), ]` and hands
the result to `isotree`'s `predict`. When a reference cell type has **no** cells
in the query, that subset has zero rows and `predict` rejects it.

A reference cell type absent from the query is an ordinary situation — it is
precisely the scenario the old app's demo dataset was built to illustrate, with
promonocytes present in the reference and missing from the query.

**Fixed** in `scDiagnostics/R/detectAnomaly.R`:

- guard the empty query subset and report `numeric(0)` query scores for that
  cell type, keeping its reference results — the same idiom
  `calculateReconstructionError()` already used;
- add `drop = FALSE` to both subsets, so a single-cell subset stays a matrix;
- skip, with a warning, any cell type with fewer than two reference cells,
  which an isolation forest cannot be fitted to.

Regression tests are in
`scDiagnostics/tests/testthat/test-detectAnomaly-missing-celltype.R`. The same
pattern was checked across the rest of the package; `detectAnomaly()` was the
only function missing the guard.

## 3. The old app's data cost

| | value |
|---|---|
| bundled `.rds` files | 36.1 MB on disk |
| read at process start | ~2 s, before the first user sees anything |
| resident in memory | **193.3 MB per R process** |
| dimensions | 23,341 genes × up to 1,518 cells, `counts` **and** `logcounts` |

The package's own curated datasets are 392–943 genes and **4.2 MB** in memory
for reference plus query — the same biology at 2% of the footprint, because
they are already reduced to the genes the diagnostics use.

The new app bundles no data at all: presets come from `data()` in the installed
package, so they cannot drift from the package version they demonstrate.

## 4. Two further traps, found while building

Neither is visible from reading the old code; both came out of testing.

**Subsetting an SCE silently discards the PCA.** `sce[, idx]` drops the
`rotation`, `percentVar` and `varExplained` attributes carried on each
`reducedDim`. Any diagnostic run on the result then fails inside
`argumentCheck()` with a message about rotation matrices that does not mention
subsetting. The app routes every subset through `subset_cells()`, which
restores them, and the generated R code shows the same fix rather than hiding
it. This affects users too, not just the app.

**`comparePCA()` and `comparePCASubspace()` need matched gene sets.** Both
compare rotation matrices gene by gene, so they cannot run when reference and
query PCAs were fitted on different genes — which is the case for the MERFISH
preset. The app checks up front and explains, instead of surfacing "The genes
in the rotation matrices differ."

## 5. Corrections to the plan

- **Exporting `downsampleSCE()` turned out to be unnecessary.** Nearly every
  exported function already takes `max_cells_query` / `max_cells_ref` and
  downsamples internally, which is the package's supported mechanism and what
  the app uses. The app does its own stratified subsample only once, at upload,
  to bound memory for the whole session.
- **Precomputing was scoped down.** Warming every candidate diagnostic for
  every preset produced 34.7 MB to save 158 s — 4.6 s per MB, the same kind of
  bargain the old app's bundled data was. Applying an explicit rule (cache only
  results that take ≥ 1 s and serialise to ≤ 3 MB) gives 3.5 MB for the same
  158 s, or 45 s per MB. A `ggplot` object carries the data it was built from,
  so `boxplotPCA` weighs 17 MB and takes 0.2 s; the Wasserstein null is 30 KB
  and takes half a minute.

## 6. Baseline timings to beat

Single runs, one machine, from `phase0-results.md`. These were the numbers the
rebuild was measured against.

| diagnostic | old app (23,341 genes) | package data (392 genes) |
|---|---:|---:|
| `detectAnomaly` | 4.7 s | 0.35 s |
| `calculateWassersteinDistance` | 22.9 s | 8.9 s |
| `plotCellTypePCA` | 2.0 s | 0.36 s |
| `calculateGraphIntegration` | 3.8 s | 1.1 s |
| `projectPCA` | 4.1 s | 0.12 s |

Most of that difference is the data, not the code: the old app ran every
diagnostic against a 23,341-gene matrix. On top of it, the new app caches
preset results across sessions (a warm Wasserstein is 0.02 s instead of 18.6 s)
and runs the slow diagnostics in a background worker so one user's job does not
block another's.
