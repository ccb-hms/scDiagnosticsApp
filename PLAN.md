# scDiagnostics App — Redesign Plan

Companion app for [scDiagnostics](https://github.com/ccb-hms/scDiagnostics)
(Briefings in Bioinformatics 27(5), bbag496).

Goal: an app that is fast, looks like a modern research tool, lets users bring
their own data, exposes the full package, and provably agrees with the package
and the paper.

---

## Status: built

The plan below was carried out. The app is in [`app/`](app/); see the
[README](README.md) for how to run and test it, and
[`audit/findings.md`](audit/findings.md) for what Phase 0 turned up.

| | planned | delivered |
|---|---|---|
| Functions exposed | 33 of 33 | 33 of 33 |
| Forked package code | 0 lines | 0 lines (`Shiny/support/` deleted) |
| Bundled data | 0 MB | 0 MB; presets come from `data()` |
| Reproducible code per panel | yes | yes, and **verified**: `test-registry` evaluates the code each panel shows and compares every numeric value against what the app computed — 33 of 33 match exactly |
| Async | `ExtendedTask` | yes, verified to run in a separate process and agree with the foreground |
| Caching | two-tier | yes, plus a 3.5 MB shipped cache saving 158 s of cold-start compute |
| Tests | smoke + parity | six suites, including a 64-check headless-browser run |

Three things the plan did not anticipate:

1. **A real bug in `scDiagnostics` itself.** `detectAnomaly()` failed whenever a
   reference cell type had no cells in the query — the exact scenario the old
   app's demo data was built to show. Fixed upstream with regression tests; see
   [`audit/findings.md`](audit/findings.md) §2.
2. **Two of the six comparable forked diagnostics gave materially different
   numbers** from the package (Wasserstein distance and graph integration), not
   merely stale ones. That is the evidence for deleting rather than maintaining
   the fork.
3. **Subsetting a `SingleCellExperiment` silently drops the PCA rotation**,
   which breaks any diagnostic downstream. Handled throughout, and surfaced in
   the generated R code rather than hidden.

Two plan items were revised on evidence, both recorded in `audit/findings.md` §5:
exporting `downsampleSCE()` proved unnecessary, and precomputing was cut from
34.7 MB to 3.5 MB for the same benefit.

Still open, because they need answers from the deployment host: Connect's
RAM/CPU and upload limits (the caps are currently env-configurable defaults),
and whether `zellkonverter` installs there (h5ad support degrades cleanly to
rds-only if not). One new constraint: the app needs `scDiagnostics` ≥ 1.7.6,
which is newer than the current Bioconductor release — see **Package version**
in the README.

---

## 1. Diagnosis

Measured against the current `Shiny/` tree at commit `a52a226`.

### 1.1 The app does not run the published package

`Shiny/support/` is **3,648 lines — 55% of the app's 6,625 lines** — of
hand-forked copies of package functions. They have drifted:

| Function | App copy | Package |
|---|---:|---:|
| `detectAnomaly` | 157 lines | 367 lines |
| `calculateDiscriminantSpace` | 234 | (package version) |
| `calculateGraphIntegration` | 507 | (package version) |
| plotting helpers (`plot*_app.R`) | 2,148 | S3 `plot()` methods |

Consequences: results shown in the app can differ from the package and from the
figures in the paper; every package bug fix has to be re-applied by hand; the
app silently freezes whatever the package looked like in mid-2025.

This is the single most important thing to fix, and it deletes more code than
anything else we add.

### 1.2 Data is ~60× heavier than it needs to be

`Shiny/data/` is **37.9 MB**, three SCEs at **23,341 genes** carrying both
`counts` and `logcounts`, all `readRDS()`'d eagerly in `global.R` at process
start.

The package already ships curated datasets — **5.6 MB total, 392–943 genes,
PCA precomputed with rotation and percentVar**:

| Dataset | Genes × cells | Why it matters |
|---|---|---|
| `reference_data` / `query_data` / `qc_data` | 392×1500 / 392×503 / 500×750 | Intro vignette; `qc_data` has QC + annotation scores |
| `zeisel_reference_data` / `zeisel_query_data` | 230×2103 / 230×902 | **Ground truth labels** — lets the app show the benchmark |
| `covid_reference_data` / `covid_query_data` | 414×900 / 414×1330 | Paper case study 3 |
| `merfish_reference_data` / `merfish_query_data` | 943×1100 / 943×1300 | Paper case study 4, spatial |

The app currently offers none of these. It ships a bespoke marrow-myeloid subset
instead, so a visitor who read the paper cannot reproduce a single paper figure.

### 1.3 No caching, no async

No `bindCache`, `memoise`, `future`, or `promises` anywhere. Every button press
recomputes from scratch, and because Shiny is single-threaded per R process, one
user running graph integration blocks every other user sharing that process on
Posit Connect.

### 1.4 Coverage gaps

8 of 33 exported functions are exposed. Missing, among others,
**`calculateGeneShifts()` — step 3 of the package's own three-step workflow**,
plus `calculateReconstructionError`, `comparePCA`, `comparePCASubspace`,
`calculateSIRSpace`, `boxplotPCA`, `plotCellTypeMDS`, `calculateCellDistances`,
all four statistical tests, and the entire QC/annotation-score family.

### 1.5 Structure and UX

- No Shiny modules. `server_main.R` carries **seven near-identical `observe`
  blocks (~200 lines)** updating cell-type-column choices, one per tab, because
  each tab hand-namespaces its inputs (`disc_`, `wass_`, `marker_`, …).
- The user re-picks reference dataset, query dataset, both cell-type columns,
  assay, and PC subset **on every single tab**.
- `shinydashboard` pins the app to Bootstrap 3 / AdminLTE 2 — effectively frozen,
  no dark mode, poor mobile behaviour.
- Results are plots plus `verbatimTextOutput` console dumps. No tables, no
  export beyond a PNG, no reproducible code, no report.
- No tests, no CI. `manifest.json` is hand-maintained and has needed emergency
  dependency fixes twice in the last three commits.

---

## 2. Decisions taken

| Decision | Choice |
|---|---|
| Function coverage | **Tiered**: ~12 curated core panels + a declarative registry that auto-generates UI for the remaining ~21 |
| Deployment | **Posit Connect only** (`ccb.connect.hms.harvard.edu`), driven by `renv` + generated manifest |
| Upstream changes | **Allowed** — fix gaps in `scDiagnostics` itself, pin the app to that version |
| Migration | **Fresh `app/` directory**, port panel by panel, keep `Shiny/` deployed until parity, then delete |

---

## 3. Target architecture

```
scDiagnosticsApp/
├── app.R                     # thin: shinyApp(ui, server)
├── renv.lock                 # pinned deps; manifest.json generated from this
├── manifest.json             # GENERATED — never hand-edited
├── R/                        # auto-sourced by Shiny
│   ├── registry.R            # declarative diagnostic registry (see §4)
│   ├── params.R              # param_int / param_choice / param_genes / ...
│   ├── data_context.R        # the single shared reactive data contract
│   ├── cache.R               # cache key construction + cachem backends
│   ├── theme.R               # bslib theme, thematic, ggplot defaults
│   ├── presets.R             # preset catalogue built on package data()
│   ├── repro.R               # deparse a registry call into runnable R
│   ├── mod_data.R            # data selection + upload + pre-flight
│   ├── mod_diagnostic.R      # GENERIC module driven by a registry entry
│   ├── mod_audit.R           # guided Project → Detect → Characterize
│   ├── mod_report.R          # Quarto report assembly + download
│   └── mod_<core>.R          # ~12 bespoke overrides for Tier 1
├── inst/
│   ├── quarto/report.qmd     # report template
│   └── precomputed/          # warmed cache for presets (build artifact)
├── data-raw/precompute.R     # builds inst/precomputed/
├── tests/testthat/           # shinytest2 + package-parity tests
└── www/                      # css, logo, favicon
```

**`support/` is deleted.** The app calls `scDiagnostics::` directly. If a
package function's output is not app-friendly, the fix goes upstream (§8), not
into a fork. A parity test (§9) enforces this permanently.

### 3.1 The data context — fixes §1.5 in one object

One reactive replaces the seven duplicated observers and the per-tab re-picking:

```r
# R/data_context.R
data_context <- function(ref, query, ref_col, query_col,
                         cell_types, assay, pc_subset,
                         ref_id, query_id) { ... }
```

`ref_id` / `query_id` are stable identifiers (`"covid@scDiagnostics-1.7.4"` for
a preset, a content hash for an upload) used for cache keys and for the
reproducible-code output. Chosen once in a persistent sidebar; every panel
inherits it and may override locally without touching the global state.

### 3.2 Caching

Every expensive reactive gets `bindCache()` keyed on
`list(fn_id, ref_id, query_id, ref_col, query_col, cell_types, pc_subset, assay, params)`.

Two backends, deliberately separated:

- **Presets** → app-level `cachem::cache_disk()`, shared across all sessions and
  pre-warmed from `inst/precomputed/`. First visit to a preset panel is instant.
- **Uploads** → session-scoped `cachem::cache_mem()`, never written to shared
  disk. A shared disk cache keyed on user data would leak one user's results to
  another; this split prevents that by construction.

### 3.3 Async

Long jobs (graph integration, isolation forests, Wasserstein resampling, gene
shifts) run under `shiny::ExtendedTask` with `bslib::input_task_button`, on
`future::plan(multisession)`. The UI stays live, the button shows its own
progress state, and one user no longer blocks the process. `future.globals.maxSize`
raised to accommodate SCE objects crossing to workers.

---

## 4. The registry — how we get all 33 functions without 33 tabs

Nearly every `scDiagnostics` function shares a signature shape
(`query_data`, `reference_data`, `*_cell_type_col`, `cell_types`, `pc_subset`,
`assay_name`) and returns an object with an S3 `plot()` method. That regularity
is the lever.

```r
# R/registry.R
diag_entry(
  id       = "wasserstein",
  fn       = scDiagnostics::calculateWassersteinDistance,
  label    = "Wasserstein distance",
  category = "alignment",
  tier     = 1,
  needs    = c("ref", "query"),
  params   = list(
    param_int("n_resamples", "Resamples", value = 300, min = 50, max = 1000)
  ),
  render   = c("plot", "table"),
  heavy    = TRUE,          # -> ExtendedTask + task button
  doc      = "calculateWassersteinDistance"
)
```

`mod_diagnostic` consumes an entry and derives, with no per-function code:

1. **The control panel** — each `param_*` knows how to render itself
   (`as_input(p, ns)`) and how to read itself back (`read(p, input)`).
2. **The call** — arguments assembled from the data context plus params.
3. **The cache key** — from the same assembled arguments.
4. **The output** — `plot()` method, `DT` table of the result's tidy slots, or both.
5. **The help text** — pulled at build time from the package's own Rd via
   `tools::Rd_db("scDiagnostics")`, so app documentation cannot drift from
   package documentation, plus a link to the pkgdown reference page.
6. **The reproducible R code** — `deparse()` of the exact call, shown in a
   copyable block under every result.

Points 5 and 6 are free consequences of the design, and 6 is a genuine
credibility feature for a paper companion.

**Tier 1 (~12)** entries additionally supply a bespoke `ui`/`server` override
for narrative layout, interactive plotly/DT, and cross-panel linking. **Tier 2
(~21)** run on the generic module alone. Adding a future package function to the
app becomes a ~10-line registry entry.

---

## 5. Information architecture

`bslib::page_navbar()` with a persistent sidebar holding the data context.

- **Data** — preset picker, upload, pre-flight report (§6). Sets the context once.
- **Audit** ⭐ — the flagship. The paper's three-step workflow on one scrolling
  page: **Project** (`plotCellTypePCA`) → **Detect** (`detectAnomaly`) →
  **Characterize** (`calculateGeneShifts`), with `value_box` headline numbers
  (cells flagged, worst cell type, top shifted genes) and a plain-language
  verdict. A visitor who read the paper reproduces its central figure in three
  clicks.
- **Diagnostics** — navbar menu grouped by the package's own five categories, so
  the app, the vignette, and the pkgdown index all agree:
  - Visualization · Alignment & statistics · Anomaly & distances ·
    Marker genes · QC & annotation scores
- **Report** — tick the panels you want, download a Quarto HTML/PDF plus the
  standalone `.R` script that regenerates every figure outside the app.
- **About** — citation, DOI, package/app version, links, `sessionInfo()`.

Visual direction: `bslib` Bootstrap 5, a single accent palette shared with the
pkgdown site (`flatly`), `thematic::thematic_shiny()` so ggplot output matches
the theme automatically, light/dark toggle, `bslib::card` + `value_box` instead
of `box()` + `verbatimTextOutput`, skeleton loaders instead of spinners.

---

## 6. Data: presets and uploads

### 6.1 Presets

Drop all 37.9 MB from `Shiny/data/`. Presets are built from `data()` in the
installed package, so they can never fall out of sync with it:

| Preset | Source | Hook |
|---|---|---|
| Bone marrow (intro) | `reference_data` / `query_data` | The intro vignette, expert vs. SingleR labels |
| Bone marrow + QC | `qc_data` | Unlocks the QC/annotation-score panels |
| Brain (Zeisel) | `zeisel_*` | **Ground truth** — show diagnostics scored against truth |
| COVID-19 monocytes | `covid_*` | Paper case study 3 |
| MERFISH colitis | `merfish_*` | Paper case study 4, spatial (SpatialExperiment) |

Optionally keep one large marrow-myeloid preset, pruned to HVGs and `logcounts`
only, as the "realistic size" demo — but loaded lazily on selection, never at
startup.

### 6.2 Uploads

Current upload accepts `.rds`, checks `is(SingleCellExperiment)`, and stops
there. New pipeline:

1. **Formats**: `.rds` (SCE / SpatialExperiment) and **`.h5ad` via
   `zellkonverter`** — most users arriving from Scanpy have h5ad, and that is a
   real adoption barrier today.
2. **Pre-flight report card** before anything runs: dimensions, assays,
   reducedDims, and **gene overlap between reference and query** — the single
   most common silent failure, currently invisible until a downstream error.
3. **Cell-type column auto-detection**: propose `colData` columns that are
   character/factor with 2–100 levels, ranked, instead of dumping every column.
4. **PCA**: if no valid PCA, run `processPCA(n_hvgs = 2000, max_cells = ...)`
   once, in the data context, and reuse it everywhere.
5. **Guardrails**: file-size cap tuned to the Connect instance; above ~20,000
   cells offer `downsampleSCE()` with a slider and state plainly what is being
   dropped; above ~5,000 genes subset to HVGs for PCA.
6. **Errors as guidance**: every validation failure names the problem and the
   exact R line that fixes it.

---

## 7. Phases

Each phase ends deployable. Nothing half-migrated reaches production, because
`Shiny/` stays live until Phase 5.

### Phase 0 — Baseline and foundation
- Profile the current app with `profvis` on a fixed script (preset load, PCA
  panel, anomaly panel, graph integration) and **record the numbers**, so later
  performance claims are measured rather than asserted.
- Create `app/`, `renv` init, `bslib` skeleton that boots with the package
  loaded and a single preset.
- Verify each of the 8 existing diagnostics against the real package function on
  a preset, and write down every discrepancy the fork introduced.

**Done when:** baseline timings are committed, and the discrepancy list exists.

### Phase 1 — Core engine
- `data_context.R`, `cache.R`, `params.R`, `registry.R`, `mod_diagnostic.R`.
- `mod_data.R` with presets + upload + pre-flight.
- Async harness and one heavy function running under `ExtendedTask`.

**Done when:** one Tier-2 function works end to end through the generic module
with caching, async, help text, and reproducible code — proving the registry.

### Phase 2 — Shell and Tier 1
- Full navbar IA, theme, dark mode.
- ~12 core panels with bespoke UI, replacing the 8 forked ones and adding
  `calculateGeneShifts`, `calculateReconstructionError`, `comparePCASubspace`, QC.
- **Delete `Shiny/support/`.**

**Done when:** every Tier-1 panel matches a direct package call (parity test
green) and the 3,648 forked lines are gone.

### Phase 3 — Full coverage
- Registry entries for the remaining ~21 exported functions.
- Category landing pages with a "which diagnostic do I need?" chooser mirroring
  the vignette's function-finder section.

**Done when:** all 33 exported functions are reachable.

### Phase 4 — Audit and report
- Guided three-step Audit page with value boxes and verdict.
- Quarto report export and standalone `.R` script export.
- `data-raw/precompute.R` warming the preset cache; ship `inst/precomputed/`.

**Done when:** a first-time visitor reproduces a paper figure in three clicks,
and preset panels render from warm cache.

### Phase 5 — Harden and ship
- `shinytest2` smoke tests per preset per panel; package-parity tests; `lintr`.
- GitHub Actions: lint, test, `renv::snapshot()` →
  `rsconnect::writeManifest()` → deploy to Connect on `main`. **`manifest.json`
  is never hand-edited again.**
- Rewrite `README.md`; add an in-app tour; update the pkgdown site and the paper
  repo to link the new app.
- **Delete `Shiny/`.**

---

## 8. Upstream `scDiagnostics` changes

Small, and they benefit package users too:

1. **Export `downsampleSCE()`** — it exists in `R/downsampleSCE.R` but is absent
   from `NAMESPACE`, so the app cannot use it for upload guardrails.
2. **Audit the `plot()` methods** against what the 2,148 lines of `plot*_app.R`
   added. Anything genuinely useful there (facet options, per-cell-type
   faceting, annotation) becomes a documented argument on the S3 method.
3. **Return tidy data alongside plots** so the app can render interactive
   plotly/DT versions without reimplementing any statistics.
4. Optional: a progress callback on `processPCA()` and `calculateGraphIntegration()`
   so the app can show real progress instead of an indeterminate spinner.

The app pins to the resulting version; Connect installs that version.

---

## 9. Testing

- **Parity tests** — for every registry entry, the app's computation path and a
  direct `scDiagnostics::` call on the same preset produce identical results.
  This is the permanent guard against the fork returning.
- **`shinytest2` smoke tests** — app boots; each preset loads; each panel
  renders without error; upload of a fixture `.rds` and `.h5ad` succeeds.
- **Registry validation test** — every entry's declared params match the actual
  formals of its function, so a package signature change fails CI instead of
  failing a user.
- **Performance regression** — assert Phase 0 baselines are beaten on the
  preset path.

---

## 10. Expected outcome

| | Now | After |
|---|---|---|
| App code | 6,625 lines, 55% forked package code | ~1,800 lines, 0 forked |
| Bundled data | 37.9 MB eager at startup | 0 MB; presets from `data()`, lazy |
| Functions exposed | 8 of 33 | 33 of 33 |
| Paper case studies | none | COVID, MERFISH, Zeisel ground truth |
| Caching / async | none / none | disk + memory / `ExtendedTask` |
| Correctness vs. package | unverified, known to differ | enforced by parity tests in CI |
| Deploy | hand-edited `manifest.json` | `renv` → generated manifest → CI |

---

## 11. Risks

- **Performance gains are projected, not measured.** Phase 0 exists to fix that
  before any claim is made publicly.
- **Connect resource limits** (RAM per process, worker count) will bound upload
  size and async worker count. Needs a real number from the Connect admin before
  Phase 1 sets the caps.
- **`zellkonverter` on Connect** pulls a basilisk Python environment; confirm it
  installs on that host before committing to h5ad support, and degrade to
  rds-only if not.
- **Bioconductor release cadence** — upstream changes land in devel and reach
  release on Bioc's schedule. The app pins an explicit version so it is never
  caught mid-cycle.
- **Scope**: Phase 3 is the easiest place to overrun. If time is short, ship
  Phases 0–2 plus 4; Tier 2 coverage is additive and can land later without
  rework, which is precisely why the registry exists.

---

## 12. Open items

- Connect instance RAM/CPU limits and max upload size.
- Whether to keep a large "realistic size" preset, and if so which tissue.
- Whether the Report export should target HTML only or HTML + PDF (PDF needs
  TinyTeX on the Connect host).
