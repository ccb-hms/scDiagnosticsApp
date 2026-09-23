# scDiagnostics App

Interactive companion to the
[scDiagnostics](https://github.com/ccb-hms/scDiagnostics) Bioconductor package.
It checks whether the cell type labels in your data hold up when compared
against a reference you trust.

**Use it here:** <https://ccb.connect.hms.harvard.edu/scDiagnosticsApp/> —
nothing to install.

## What it does

Automated annotation is only as good as the match between your data and the
reference it was labelled from. The scDiagnostics App compares the two and
shows you where the labels look right and where they do not.

- **Guided audit** — a three-step walkthrough: project your data onto the
  reference, find the cells that do not look like the type they were assigned,
  and see which genes explain the difference. It ends in a plain-language
  verdict.
- **All 33 diagnostics** from the package, each with figures, the numbers as
  sortable tables, a note on how to read the result, and the R code that
  produces it.
- **Four ready-made datasets** to explore, including the COVID-19 and MERFISH
  case studies from the paper and a mouse brain dataset with known-correct
  labels.
- **Your own data** — a `SingleCellExperiment` or `SpatialExperiment` saved as
  `.rds`, or an AnnData `.h5ad`. The app tells you what it found before it runs
  anything, and computes a PCA if your object does not have one.
- **A report you can keep** — one self-contained HTML file with every figure,
  plus an R script that recreates them outside the app.

Every number and every figure comes from `scDiagnostics` itself. The app
implements no statistics of its own, so what you see here is what you would get
in an R session.

## Citation

Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026).
*scDiagnostics: systematic assessment of cell type annotation in single-cell
transcriptomics data.* **Briefings in Bioinformatics**, 27(5), bbag496.
doi:[10.1093/bib/bbag496](https://doi.org/10.1093/bib/bbag496)

## Running your own copy

Most people should just use the live app. To run it locally you need R 4.4 or
newer and `scDiagnostics` 1.7.6 or newer. From the repository root:

```bash
Rscript run.R
```

`run.R` checks your setup before starting and says plainly what is missing.
That check matters when several versions of R are installed, since the shell
often picks a different one than you expect. If it reports the wrong R, call
the one you want directly:

```powershell
& "C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe" run.R
```

These environment variables change how it behaves:

| variable | default | meaning |
|---|---|---|
| `SCDIAG_MAX_UPLOAD_MB` | 500 | largest file accepted |
| `SCDIAG_MAX_CELLS` | 20000 | larger uploads are subsampled, keeping cell type proportions |
| `SCDIAG_ASYNC` | `true` | run slow diagnostics in background workers |
| `SCDIAG_WORKERS` | 2 | number of workers |
| `SCDIAG_CACHE_DIR` | temp dir | shared cache for the ready-made datasets |
| `SCDIAG_CACHE_MAX_MB` | 1024 | size cap for that cache |

---

The rest of this file is for people working on the app itself.

## Layout

```
app/                  the application
  app.R               entry point
  R/                  modules, registry, theme, cache (auto-sourced by Shiny)
  www/styles.css
  inst/precomputed/   warmed cache for the ready-made datasets
  data-raw/           script that builds that cache
  tests/              seven test suites, from a parse check to a browser run
  manifest.json       generated for Posit Connect; never hand-edited
audit/                baseline of the previous app, and its results
run.R                 local launcher
manifest.R            regenerates app/manifest.json, and checks it is deployable
PLAN.md               the redesign plan this was built from
```

## Tests

```bash
cd app && Rscript tests/run-all.R
```

| suite | what it checks |
|---|---|
| `boot` | everything parses and sources; the UI builds; every registry entry's declared arguments exist on its function |
| `test-modules` | reactive logic via `shiny::testServer`: data selection, caching tiers, context helpers |
| `test-async` | a background worker returns the same answer as the foreground, and is sent a descriptor rather than the data |
| `test-registry` | all 33 diagnostics run, **and the R code the app shows reproduces the result exactly** |
| `test-plots` | every result draws, on all four ready-made datasets; anything inapplicable is gated with an explanation instead of an error |
| `test-app` | the whole app in a headless browser, including upload and report download |
| `test-rootrun` | launching from the repository root, the way this README says |

`test-registry` is the one that matters most: it runs the code from each
panel's **R code** tab in a clean session and compares every number against
what the app computed.

## Package version and deployment

The app needs `scDiagnostics` 1.7.6 or newer and shows a banner if the
installed version is older. Bioconductor release currently ships an older
version, so the deployment installs the package from GitHub:

```r
BiocManager::install("ccb-hms/scDiagnostics")
```

`app/manifest.json` records whatever was installed on the machine that
generated it, so regenerate it from a library holding the version you intend to
deploy:

```bash
Rscript manifest.R
```

That script also checks every package it recorded can actually be downloaded at
the version written down, which catches a manifest that would fail on the
server. Regenerate it whenever you add or remove a dependency or a file under
`app/`.

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
the help text and the reproducible code from that. `boot` fails if the declared
parameters do not match the function's actual arguments.
