# Module-level tests with shiny::testServer.
#
# These exercise the reactive logic directly - no browser - so a broken
# context, a stale selection or a diagnostic that never runs is caught in
# seconds rather than through the end-to-end test.
#
# Run from the app directory:  Rscript tests/test-modules.R

suppressPackageStartupMessages({
    library(shiny); library(bslib); library(bsicons)
    library(scDiagnostics); library(SingleCellExperiment); library(ggplot2)
    library(promises); library(future)
})

for (f in sort(list.files("R", pattern = "[.][Rr]$", full.names = TRUE))) {
    source(f, local = globalenv())
}
options(sc.async.ready = FALSE)   # keep these tests in-process

failures <- character(0)
check <- function(label, cond, detail = "") {
    if (isTRUE(cond)) {
        cat(sprintf("  ok    %s\n", label))
    } else {
        failures <<- c(failures, label)
        cat(sprintf("  FAIL  %s  %s\n", label, detail))
    }
}

# --- data module -----------------------------------------------------------

cat("== data module ==\n")
shiny::testServer(data_server, args = list(id = "data"), {
    check("context starts empty", !isTRUE(session$getReturned()$ctx()$ready))

    # Load a preset the way the Home page's demo button does.
    session$getReturned()$load_preset("covid", notify = FALSE)
    session$flushReact()

    ctx <- session$getReturned()$ctx()
    check("preset loads a reference", !is.null(ctx$ref),
          "ref is NULL after load_preset")
    check("preset loads a query", !is.null(ctx$query))
    check("context is ready", isTRUE(ctx$ready),
          sprintf("ref_col = %s", ctx$ref_col %||% "NULL"))
    check("reference column defaults correctly",
          identical(ctx$ref_col, PRESETS$covid$ref$col), ctx$ref_col %||% "NULL")
    check("query column defaults correctly",
          identical(ctx$query_col, PRESETS$covid$query$col), ctx$query_col %||% "NULL")
    check("shared cell types found", length(ctx$shared_types) > 0,
          paste(ctx$shared_types, collapse = ", "))
    check("caption mentions both datasets",
          grepl("Reference", context_caption(ctx)) && grepl("Query", context_caption(ctx)))

    # Switching preset must replace, not merge.
    session$getReturned()$load_preset("marrow", notify = FALSE)
    session$flushReact()
    ctx2 <- session$getReturned()$ctx()
    check("switching preset swaps the data",
          identical(ctx2$ref_desc$key, PRESETS$marrow$ref$data), ctx2$ref_desc$key)
    check("switching preset swaps the columns",
          identical(ctx2$ref_col, PRESETS$marrow$ref$col), ctx2$ref_col %||% "NULL")

    # An explicit column choice must win over the preset default.
    session$setInputs(ref_col = "expert_annotation", assay = "logcounts",
                      pcs = c(1, 3), cell_types = character(0))
    ctx3 <- session$getReturned()$ctx()
    check("PC range follows the slider", identical(ctx3$pc_subset, 1:3),
          paste(ctx3$pc_subset, collapse = ","))
})

# --- uploads ----------------------------------------------------------------
#
# The upload path has to cope with objects that are not what the app wants:
# no PCA, too many cells, the wrong class, genes that do not match the other
# dataset. Each of those should produce an explanation, not a stack trace.

cat("\n== upload validation ==\n")

ref0 <- materialize(preset_descriptor("reference_data"))
qry0 <- materialize(preset_descriptor("query_data"))

status_of <- function(checks, label) {
    hit <- Filter(function(c) c$label == label, checks)
    if (!length(hit)) NA_character_ else hit[[1]]$status
}

ok_obj <- inspect_upload(ref0, qry0)
check("a good object passes every check",
      !upload_blocked(ok_obj) && all(vapply(ok_obj, `[[`, character(1), "status") == "ok"),
      paste(vapply(ok_obj, function(c) paste0(c$label, "=", c$status), character(1)),
            collapse = "; "))

not_sce <- inspect_upload(data.frame(a = 1))
check("a non-SCE upload is rejected with an explanation",
      upload_blocked(not_sce) && grepl("SingleCellExperiment", not_sce[[1]]$detail))

no_pca <- ref0
SingleCellExperiment::reducedDims(no_pca) <- list()
np <- inspect_upload(no_pca)
check("a missing PCA warns rather than blocks",
      !upload_blocked(np) && identical(status_of(np, "PCA"), "warn"))

raw_only <- ref0
SummarizedExperiment::assayNames(raw_only) <- "counts"
ro <- inspect_upload(raw_only)
check("a counts-only object warns about normalisation",
      identical(status_of(ro, "Assays"), "warn") &&
          grepl("logNormCounts", ro[[which(vapply(ro, `[[`, character(1), "label") == "Assays")]]$detail))

disjoint <- qry0
rownames(disjoint) <- paste0("X", seq_len(nrow(disjoint)))
dj <- inspect_upload(disjoint, ref0)
check("disjoint gene names are caught before anything runs",
      upload_blocked(dj) &&
          grepl("gene identifiers", dj[[length(dj)]]$detail))

# Rename a quarter of the genes: enough overlap to proceed, little enough to
# be worth saying out loud.
few_shared <- qry0
n_rename <- round(nrow(few_shared) * 0.25)
rownames(few_shared)[seq_len(n_rename)] <- paste0("X", seq_len(n_rename))
fs <- inspect_upload(few_shared, ref0)
check("partial gene overlap warns rather than blocks",
      identical(status_of(fs, "Gene overlap with the other dataset"), "warn"),
      status_of(fs, "Gene overlap with the other dataset"))

cat("\n== upload processing ==\n")

# An object with no PCA must come out of the pipeline usable.
processed <- scDiagnostics::processPCA(no_pca, n_hvgs = 200)
check("processPCA gives the object a usable PCA", has_valid_pca(processed))

# Subsampling must keep the object usable too, which is the step that would
# silently discard the PCA if subset_cells() were not used.
big <- subsample_sce(ref0, 300, "expert_annotation")
check("a subsampled upload still has a usable PCA", has_valid_pca(big))
check("a subsampled upload still projects",
      !is_failure(capture_run(scDiagnostics::projectPCA(
          query_data = qry0, reference_data = big,
          query_cell_type_col = "SingleR_annotation",
          ref_cell_type_col = "expert_annotation", pc_subset = 1:5))))

# A file descriptor must identify by content, so two uploads of the same
# object share cache entries and an edited file does not return a stale one.
f1 <- tempfile(fileext = ".rds"); saveRDS(ref0, f1)
f2 <- tempfile(fileext = ".rds"); saveRDS(ref0, f2)
f3 <- tempfile(fileext = ".rds"); saveRDS(qry0, f3)
check("identical uploads get the same identity",
      identical(descriptor_id(file_descriptor(f1, "a")), descriptor_id(file_descriptor(f2, "b"))))
check("different uploads get different identities",
      !identical(descriptor_id(file_descriptor(f1, "a")), descriptor_id(file_descriptor(f3, "c"))))
check("an uploaded object round-trips through its descriptor",
      identical(dim(materialize(file_descriptor(f1, "a"))), dim(ref0)))
unlink(c(f1, f2, f3))

# --- diagnostic module -----------------------------------------------------

cat("\n== diagnostic module ==\n")

fixture_ctx <- local({
    p <- PRESETS$marrow
    rd <- preset_descriptor(p$ref$data); qd <- preset_descriptor(p$query$data)
    build_context(ref = materialize(rd), query = materialize(qd),
                  ref_desc = rd, query_desc = qd,
                  ref_col = p$ref$col, query_col = p$query$col,
                  pc_subset = 1:5, preset_id = "marrow")
})

run_panel <- function(entry_id, inputs = list()) {
    entry <- REGISTRY[[entry_id]]
    out <- list()
    shiny::testServer(
        diagnostic_server,
        args = list(id = "d", entry = entry, ctx_r = shiny::reactive(fixture_ctx),
                    session_cache = new_session_cache(),
                    active = shiny::reactive(TRUE)),
        {
            if (length(inputs)) do.call(session$setInputs, inputs)
            session$flushReact()
            if (entry$heavy) {
                session$setInputs(run = 1)
                session$flushReact()
            }
            out <<- result_store()
        }
    )
    out
}

# A light panel updates on its own; a heavy one waits for the button.
light <- run_panel("plotCellTypePCA")
check("light panel computes without a click",
      !is.null(light) && !is_failure(light$value),
      if (is.null(light)) "no result" else light$value$message %||% "")

heavy <- run_panel("detectAnomaly", list(n_tree = 50))
check("heavy panel computes after the run button",
      !is.null(heavy) && !is_failure(heavy$value),
      if (is.null(heavy)) "no result" else heavy$value$message %||% "")
if (!is.null(heavy) && !is_failure(heavy$value)) {
    check("heavy panel honours its parameters",
          identical(heavy$params$n_tree, 50L), paste(heavy$params$n_tree))
    check("summary produces headline figures",
          length(REGISTRY$detectAnomaly$summarise(heavy$value)) > 0)
}

# A gated panel must not run at all.
gated_ctx <- fixture_ctx
gated <- local({
    entry <- REGISTRY$plotQCvsAnnotation
    zeisel <- PRESETS$zeisel
    rd <- preset_descriptor(zeisel$ref$data); qd <- preset_descriptor(zeisel$query$data)
    ctx <- build_context(ref = materialize(rd), query = materialize(qd),
                         ref_desc = rd, query_desc = qd,
                         ref_col = zeisel$ref$col, query_col = zeisel$query$col,
                         pc_subset = 1:5, preset_id = "zeisel")
    entry_ready(entry, ctx)
})
check("QC panel is gated on data without numeric metadata",
      !is.null(gated) && grepl("colData", gated), gated %||% "not gated")

# --- caching ---------------------------------------------------------------

cat("\n== caching ==\n")
key1 <- cache_key("detectAnomaly", fixture_ctx, list(n_tree = 500))
key2 <- cache_key("detectAnomaly", fixture_ctx, list(n_tree = 500))
key3 <- cache_key("detectAnomaly", fixture_ctx, list(n_tree = 100))
check("cache key is stable for identical inputs", identical(key1, key2))
check("cache key changes with a parameter", !identical(key1, key3))

ctx_other <- fixture_ctx
ctx_other$pc_subset <- 1:10
check("cache key changes with the component range",
      !identical(key1, cache_key("detectAnomaly", ctx_other, list(n_tree = 500))))

tier <- cache_for(fixture_ctx, new_session_cache())
check("preset results are shareable across sessions", isTRUE(tier$shareable))

upload_ctx <- fixture_ctx
tmp <- tempfile(fileext = ".rds"); saveRDS(materialize(fixture_ctx$query_desc), tmp)
upload_ctx$query_desc <- file_descriptor(tmp, "uploaded.rds")
tier2 <- cache_for(upload_ctx, new_session_cache())
check("uploaded results stay private to the session", !isTRUE(tier2$shareable))
unlink(tmp)

# --- context helpers -------------------------------------------------------

cat("\n== context helpers ==\n")
ref <- materialize(preset_descriptor("reference_data"))
check("candidate columns rank the annotation first",
      identical(candidate_celltype_cols(ref)[1], "expert_annotation"),
      paste(candidate_celltype_cols(ref), collapse = ", "))
check("valid PCA is recognised", has_valid_pca(ref))

sub <- subset_cells(ref, 1:100)
check("subsetting keeps the PCA rotation", has_valid_pca(sub))
check("subsetting keeps the cells asked for", ncol(sub) == 100)

small <- subsample_sce(ref, 200, "expert_annotation")
check("subsampling respects the cap", ncol(small) <= 200, paste(ncol(small)))
check("subsampling keeps every cell type",
      setequal(unique(small$expert_annotation), unique(ref$expert_annotation)))
check("subsampling keeps the PCA rotation", has_valid_pca(small))

cat(sprintf("\n%d failures\n", length(failures)))
if (length(failures)) {
    for (f in failures) cat("  -", f, "\n")
    quit(status = 1)
}
cat("MODULES OK\n")
