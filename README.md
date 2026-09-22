# scDiagnostics app

Interactive companion to the
[scDiagnostics](https://github.com/ccb-hms/scDiagnostics) Bioconductor package:
audit cell type annotations by comparing a query dataset against a reference
you trust.

**Live app:** <https://ccb.connect.hms.harvard.edu/scDiagnosticsApp/>

## Citation

Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026).
*scDiagnostics: systematic assessment of cell type annotation in single-cell
transcriptomics data.* **Briefings in Bioinformatics**, 27(5), bbag496.
doi:[10.1093/bib/bbag496](https://doi.org/10.1093/bib/bbag496)

## What it does

Automated annotation is only as trustworthy as the alignment between reference
and query. The app runs the package's diagnostics over both and shows where the
labels hold up and where they do not.

- **Guided audit** — the three-step workflow the paper's case studies follow:
  project the query onto the reference's PCA space, detect cells that do not
  look like their assigned type, and characterize the genes that explain the
  difference. Runs automatically on the data you have selected and ends in a
  plain-language verdict.
- **All 33 diagnostics** in the package, grouped the way the package
  documentation groups them, each with headline figures, the numbers as
  sortable tables, guidance on how to read the output, and the R code that
  reproduces it.
- **Four prepared datasets**, including the paper's COVID-19 and MERFISH case
  studies and a mouse brain dataset with ground-truth labels.
- **Your own data** — a `SingleCellExperiment` or `SpatialExperiment` as
  `.rds`, or an AnnData `.h5ad`. The app reports what it found before running
  anything, and computes a PCA if the object lacks one.
- **Report export** — a self-contained HTML file with every figure embedded,
  plus a standalone R script that regenerates them outside the app.

The app holds no copy of any statistical method: every number and every figure
comes from `scDiagnostics` itself, so results here match results in an R
session. A test suite enforces that (see below).

## Layout

```
app/                  the application
  app.R               entry point
  R/                  modules, registry, theme, cache (auto-sourced by Shiny)
  www/styles.css
  inst/precomputed/   warmed cache for the prepared datasets
  data-raw/           script that builds that cache
  tests/              seven test suites, from a parse check to a browser run
  manifest.json       generated for Posit Connect; never hand-edited
audit/                Phase 0 baseline of the previous app, and its results
run.R                 local launcher (checks the environment first)
PLAN.md               the redesign plan this was built from
```

## Running it locally

Requires R ≥ 4.4 and `scDiagnostics` ≥ 1.7.6 (see **Package version** below).
From the repository root:

```bash
Rscript run.R
```

`run.R` checks the R it is running under before starting, and if a package is
missing or `scDiagnostics` is too old it says so and points at any other R
installation on the machine. That check earns its keep: with more than one R
installed, the shell often resolves `Rscript` to the wrong one, and the
failure otherwise surfaces as a missing package deep in startup rather than as
the version problem it is. Check which one you are using with `which Rscript`
(`Get-Command Rscript` in PowerShell), and call the right one directly if
needed:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe" run.R
```

Avoid `Rscript -e '...'` for this. Shells disagree about quoting: PowerShell
strips the inner double quotes of `-e 'shiny::runApp("app")'`, leaving R with
a bare symbol and the unhelpful error `object 'app' not found`. If you do want
a one-liner, put the double quotes on the outside:

```powershell
Rscript -e "shiny::runApp('app', launch.browser = TRUE)"
```

Behaviour can be tuned with environment variables:

| variable | default | meaning |
|---|---|---|
| `SCDIAG_MAX_UPLOAD_MB` | 500 | largest file accepted |
| `SCDIAG_MAX_CELLS` | 20000 | uploads above this are subsampled, keeping cell type proportions |
| `SCDIAG_ASYNC` | `true` | run slow diagnostics in background workers |
| `SCDIAG_WORKERS` | 2 | number of workers |
| `SCDIAG_CACHE_DIR` | temp dir | shared cache for results from the prepared datasets |
| `SCDIAG_CACHE_MAX_MB` | 1024 | size cap for that cache |

## Tests

```bash
cd app && Rscript tests/run-all.R
```

Seven suites, fastest first:

| suite | what it checks |
|---|---|
| `boot` | everything parses and sources; the UI builds; every registry entry's declared arguments exist on its function |
| `test-modules` | reactive logic via `shiny::testServer`: data selection, caching tiers, context helpers |
| `test-async` | a background worker returns the same answer as the foreground, and is sent a descriptor rather than the data |
| `test-registry` | all 33 diagnostics run, **and the R code the app shows reproduces the result exactly** |
| `test-plots` | every result draws, on all four prepared datasets; anything inapplicable is gated with an explanation instead of an error |
| `test-app` | the whole app in a headless browser, including upload and report download |
| `test-rootrun` | launching the way this README says, from the repository root |

`test-registry` is the important one. It evaluates the code from each panel's
**R code** tab in a clean environment and compares every numeric value against
what the app computed. A reproducible-code panel that does not reproduce the
panel is worse than none.

## Package version

The app requires `scDiagnostics` ≥ 1.7.6 and checks at startup, showing a
banner if the installed version is older.

This matters for deployment: Bioconductor **release** currently ships an older
version than the app needs, so the host must provide Bioconductor **devel**, or
install the package from GitHub:

```r
BiocManager::install("ccb-hms/scDiagnostics")
```

The generated `manifest.json` records whichever version was installed when it
was built, so regenerate it on a machine with the intended version:

```r
options(repos = BiocManager::repositories())
rsconnect::writeManifest(appDir = "app", appPrimaryDoc = "app.R")
```

## Adding a diagnostic

Diagnostics are declared, not written. To expose a new package function, add an
entry to `app/R/21-registry-entries.R`:

```r
add_entry(diag_entry(
    id       = "calculateSomething",
    fn       = "calculateSomething",
    label    = "Something",
    category = "alignment",
    question = "The question this answers, in one line.",
    blurb    = "What it computes.",
    reading  = "How to interpret the output.",
    data_map = map_pair(),
    params   = list(param_int("n_iter", "Iterations", 100, 10, 1000)),
    heavy    = TRUE
))
```

The generic module derives the controls, the cache key, the plot, the tables,
the help text and the reproducible code from that. `boot` will fail if the
declared parameters do not match the function's actual arguments.
