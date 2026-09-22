# Preset datasets.
#
# Every preset is built from data() in the installed scDiagnostics package
# rather than from files bundled with the app, so presets can never drift from
# the package they are meant to demonstrate. Objects are loaded lazily on
# selection and memoised for the life of the process.

PKG <- "scDiagnostics"

.preset_cache <- new.env(parent = emptyenv())

#' Load a dataset shipped with scDiagnostics, memoised per process
load_pkg_data <- function(name) {
    if (!is.null(.preset_cache[[name]])) return(.preset_cache[[name]])
    e <- new.env(parent = emptyenv())
    utils::data(list = name, package = PKG, envir = e)
    obj <- get(name, envir = e)
    .preset_cache[[name]] <- obj
    obj
}

PRESETS <- list(

    marrow = list(
        id = "marrow",
        label = "Bone marrow (introductory)",
        tagline = "Expert labels vs. SingleR, with QC and confidence scores",
        description = paste(
            "Peripheral blood and bone marrow cells carrying both an expert",
            "annotation and a SingleR transfer. The only preset with QC metrics",
            "and annotation confidence scores in colData, so it is the one to",
            "use for the QC panels."
        ),
        look_for = paste(
            "Compare the expert annotation against the SingleR call for the same",
            "cells. Where the two disagree, the diagnostics should light up."
        ),
        ref = list(data = "reference_data", col = "expert_annotation"),
        query = list(data = "query_data", col = "SingleR_annotation"),
        qc = list(
            source = "query",
            cell_type_col = "SingleR_annotation",
            qc_col = "percent_mito",
            score_col = "annotation_scores",
            gene_set_col = "gene_set_scores"
        ),
        pc_subset = 1:5,
        vignette = "https://ccb-hms.github.io/scDiagnostics/articles/Introduction.html"
    ),

    zeisel = list(
        id = "zeisel",
        label = "Mouse brain, Zeisel (ground truth)",
        tagline = "Known true labels, so diagnostics can be scored",
        description = paste(
            "Mouse cortex and hippocampus with curated true cell types in both",
            "reference and query. Because the truth is known, anomaly detection",
            "can be scored rather than merely inspected."
        ),
        look_for = paste(
            "Both objects carry true_cell_type. Corrupt the query labels",
            "deliberately, or compare flagged cells against the truth, to see",
            "how detection behaves when you already know the answer."
        ),
        ref = list(data = "zeisel_reference_data", col = "true_cell_type"),
        query = list(data = "zeisel_query_data", col = "true_cell_type"),
        pc_subset = 1:5,
        vignette = "https://ccb-hms.github.io/scDiagnostics/articles/ZeiselBenchmarking.html"
    ),

    covid = list(
        id = "covid",
        label = "COVID-19 monocytes",
        tagline = "Paper case study: an interferon-driven query state",
        description = paste(
            "A healthy reference annotated by the original authors, and a severe",
            "COVID-19 query annotated by Azimuth. CD14 monocytes in the query",
            "occupy a disease state that the reference does not contain."
        ),
        look_for = paste(
            "Project the query onto the reference PCA space and look at CD14",
            "monocytes: they shift along PC1/PC2, a subset is flagged as",
            "anomalous, and the shifted genes are interferon-response genes."
        ),
        ref = list(data = "covid_reference_data", col = "author_cell_type_merged"),
        query = list(data = "covid_query_data", col = "azimuth_celltype_l1_merged"),
        pc_subset = 1:5,
        focus_cell_type = "CD14 mono",
        vignette = "https://ccb-hms.github.io/scDiagnostics/articles/COVIDCaseStudy.html"
    ),

    merfish = list(
        id = "merfish",
        label = "MERFISH colitis (spatial)",
        tagline = "Paper case study: spatial data, same diagnostics",
        description = paste(
            "MERFISH spatial transcriptomics of healthy and inflamed colon,",
            "stored as a SpatialExperiment. The diagnostics apply unchanged",
            "because SpatialExperiment extends SingleCellExperiment."
        ),
        look_for = paste(
            "The query's tier2_merged column separates inflamed fibroblasts and",
            "inflamed smooth muscle, states absent from the reference. Annotated",
            "at the coarser tier, those cells should look anomalous."
        ),
        ref = list(data = "merfish_reference_data", col = "cell_type_merged"),
        query = list(data = "merfish_query_data", col = "cell_type_merged"),
        pc_subset = 1:5,
        vignette = "https://ccb-hms.github.io/scDiagnostics/articles/MERFISHCaseStudy.html"
    )
)

#' Choices for a preset picker
preset_choices <- function() {
    stats::setNames(names(PRESETS), vapply(PRESETS, `[[`, character(1), "label"))
}

# ---------------------------------------------------------------------------
# Data source descriptors
#
# A descriptor says *where* data comes from without holding it. The main
# process and any future worker can both materialise one, so a background job
# never has to receive a serialised SingleCellExperiment: a preset is loaded
# from the package inside the worker, and an upload is read from the file the
# session already wrote to disk. Serialising a large SCE to every worker is the
# thing that makes naive async slower than running synchronously.
# ---------------------------------------------------------------------------

#' @param kind "preset" or "file"
#' @param key dataset name (preset) or file path (file)
descriptor <- function(kind, key, label = key) {
    structure(list(kind = kind, key = key, label = label), class = "sc_descriptor")
}

preset_descriptor <- function(dataset_name) {
    descriptor("preset", dataset_name, dataset_name)
}

file_descriptor <- function(path, label) {
    descriptor("file", normalizePath(path, winslash = "/", mustWork = TRUE), label)
}

#' Stable identity of a descriptor, for cache keys.
#'
#' Presets are identified by name plus package version, so a package upgrade
#' correctly invalidates cached results. Files are identified by content hash,
#' so two uploads of the same object share cache entries and a modified file
#' never returns a stale result.
descriptor_id <- function(d) {
    stopifnot(inherits(d, "sc_descriptor"))
    if (d$kind == "preset") {
        paste0("preset:", d$key, "@", as.character(utils::packageVersion(PKG)))
    } else {
        paste0("file:", substr(digest::digest(file = d$key, algo = "xxhash64"), 1, 16))
    }
}

#' Materialise a descriptor into a SingleCellExperiment.
#'
#' Safe to call inside a future worker; memoised per process.
materialize <- function(d) {
    stopifnot(inherits(d, "sc_descriptor"))
    if (d$kind == "preset") return(load_pkg_data(d$key))
    key <- paste0("file:", d$key)
    if (!is.null(.preset_cache[[key]])) return(.preset_cache[[key]])
    obj <- readRDS(d$key)
    .preset_cache[[key]] <- obj
    obj
}

#' Is this descriptor shareable across users?
#'
#' Only preset results may go in the cross-session disk cache. Results derived
#' from an uploaded file stay in that session's memory, because the disk cache
#' is shared by every visitor to the app.
is_shareable <- function(d) inherits(d, "sc_descriptor") && d$kind == "preset"
