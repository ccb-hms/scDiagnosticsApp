# Registry entries: every exported scDiagnostics diagnostic.
#
# Defaults match the package defaults, so a user who changes nothing sees what
# the vignettes show. Parameters that exist mainly for tuning are marked
# advanced and collapse out of the way.

# --- shared summary helpers ------------------------------------------------

#' Headline figure for a value box
figure <- function(label, value, caption = NULL, tone = "default") {
    list(label = label, value = value, caption = caption, tone = tone)
}

#' Per-cell-type anomaly counts, shared by detectAnomaly and reconstruction error
anomaly_table <- function(x, flag_field, score_field) {
    types <- setdiff(names(x), "Combined")
    if (!length(types)) return(list())
    rows <- lapply(types, function(ct) {
        el <- x[[ct]]
        q <- el[[flag_field]]
        s <- el[[score_field]]
        data.frame(
            cell_type = ct,
            query_cells = length(q),
            flagged = sum(q, na.rm = TRUE),
            flagged_pct = if (length(q)) 100 * mean(q, na.rm = TRUE) else NA_real_,
            median_score = if (length(s)) stats::median(s, na.rm = TRUE) else NA_real_,
            threshold = el[["applied_threshold"]] %||% NA_real_,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(-out$flagged_pct), , drop = FALSE]
}

#' Worst offender across cell types, used for the plain-language verdict
worst_type <- function(df, col = "flagged_pct") {
    df <- df[!is.na(df[[col]]) & df$query_cells > 0, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    df[which.max(df[[col]]), , drop = FALSE]
}

#' Reproducible-code prelude for diagnostics restricted to one cell type
#'
#' Emits the attribute-restoring helper as well as the subsetting, because
#' subsetting an SCE drops the PCA rotation these functions require. Users hit
#' this too, so the generated script shows the fix rather than hiding it.
single_type_prelude <- function(ctx, args, param_values, ref_var, query_var) {
    ct <- (param_values$cell_types %||% ctx$cell_types)[1]
    if (is.null(ct)) return(NULL)
    c(
        "# Subsetting a SingleCellExperiment drops the PCA rotation and",
        "# percentVar attributes that the projection needs, so restore them.",
        "keep_pca <- function(sub, orig) {",
        '    m <- reducedDim(sub, "PCA")',
        '    a <- attributes(reducedDim(orig, "PCA"))',
        '    attr(m, "rotation") <- a[["rotation"]]',
        '    attr(m, "percentVar") <- a[["percentVar"]]',
        '    reducedDim(sub, "PCA") <- m',
        "    sub",
        "}",
        "",
        sprintf("# %s compares a single cell type at a time.",
                ctx$entry_label %||% "This diagnostic"),
        sprintf('query <- keep_pca(%s[, %s[["%s"]] == "%s"], %s)',
                query_var, query_var, ctx$query_col, ct, query_var),
        sprintf('reference <- keep_pca(%s[, %s[["%s"]] == "%s"], %s)',
                ref_var, ref_var, ctx$ref_col, ct, ref_var)
    )
}

REGISTRY <- list()

# This file is sourced into one environment and fills REGISTRY as it goes.
# add_entry() names that environment, rather than letting `<<-` search for
# REGISTRY up the scope chain from wherever it happens to be called.
registry_env <- environment()

add_entry <- function(e) {
    entries <- get("REGISTRY", envir = registry_env)
    entries[[e$id]] <- e
    assign("REGISTRY", entries, envir = registry_env)
    invisible(e)
}

# ===========================================================================
# Visualization
# ===========================================================================

add_entry(diag_entry(
    id = "plotCellTypePCA",
    fn = "plotCellTypePCA",
    label = "PCA projection",
    category = "visualization",
    tier = 1L,
    question = "Does the query sit where the reference says it should?",
    blurb = paste(
        "Projects the query onto the reference's PCA space and draws every pair",
        "of principal components, coloured by cell type."
    ),
    reading = paste(
        "Each panel is one pair of principal components; the diagonal shows the",
        "distribution along a single component. A cell type whose query cloud",
        "(dashed) sits on top of its reference cloud (solid) is well aligned. A",
        "shifted or broadened cloud is the first sign that the annotation for",
        "that cell type deserves a closer look. Shifts along a high-variance",
        "component matter more than shifts along a low-variance one, so read the",
        "percentage in each axis label."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20),
        param_celltypes(),
        param_choice("lower_facet", "Lower panels",
                     c("scatter", "contour", "ellipse", "blank"), "scatter"),
        param_choice("diagonal_facet", "Diagonal panels",
                     c("ridge", "density", "boxplot"), "ridge"),
        param_choice("upper_facet", "Upper panels",
                     c("blank", "scatter", "contour", "ellipse"), "blank"),
        param_int("max_cells_query", "Max query cells plotted", 2000, 100, 20000,
                  help = "Plotting every cell is slow and hides density. Sampling does not change the projection.",
                  advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells plotted", 2000, 100, 20000, advanced = TRUE)
    ),
    returns = "plot",
    plot_height = 640
))

add_entry(diag_entry(
    id = "boxplotPCA",
    fn = "boxplotPCA",
    label = "PC score distributions",
    category = "visualization",
    question = "How does each cell type spread along each component?",
    blurb = "Boxplots or violins of principal component scores, split by cell type and dataset.",
    reading = paste(
        "A reference and query box that do not overlap on a high-variance",
        "component indicates a systematic shift for that cell type, which is",
        "usually a batch effect or a genuine state difference."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20),
        param_celltypes(),
        param_choice("shape", "Shape", c("box", "violin"), "box"),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    returns = "plot"
))

add_entry(diag_entry(
    id = "plotCellTypeMDS",
    fn = "plotCellTypeMDS",
    label = "MDS projection",
    category = "visualization",
    question = "How do cell types relate in a distance-preserving 2-D view?",
    blurb = "Multidimensional scaling of reference and query cells, coloured by cell type.",
    reading = paste(
        "MDS lays cells out so that plotted distance approximates expression",
        "distance. Unlike PCA it has no interpretable axes, so read only the",
        "relative positions: which cell types sit near each other, and whether",
        "query cells land inside their reference neighbourhood."
    ),
    data_map = map_pair(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_int("max_cells_query", "Max query cells", 5000, 100, 20000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    returns = "plot",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateDiscriminantSpace",
    fn = "calculateDiscriminantSpace",
    label = "Discriminant space",
    category = "visualization",
    tier = 1L,
    question = "In the space that best separates reference cell types, where does the query fall?",
    blurb = paste(
        "Fits a discriminant projection on the reference, chosen to maximise",
        "separation between cell types, then projects the query into it."
    ),
    reading = paste(
        "PCA finds directions of maximum variance, which need not separate cell",
        "types. This projection is fitted specifically to separate them, so it",
        "is the more sensitive view for misannotation: a query cell landing in",
        "the wrong cluster here is strong evidence, because the axes were built",
        "to tell those clusters apart. The eigenvalues say how much separating",
        "power each discriminant carries."
    ),
    data_map = map_ref_first(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_int("n_tree", "Trees", 500, 50, 2000, advanced = TRUE),
        param_int("n_top", "Top genes", 20, 5, 200, advanced = TRUE),
        param_num("eigen_threshold", "Eigenvalue threshold", 0.1, 0, 1, 0.01, advanced = TRUE),
        param_bool("calculate_metrics", "Compute Mahalanobis / cosine metrics", FALSE, advanced = TRUE),
        param_num("alpha", "Alpha", 0.01, 0.001, 0.2, 0.001, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("lower_facet", "Lower panels",
                     c("scatter", "contour", "ellipse", "blank"), "scatter"),
        param_choice("diagonal_facet", "Diagonal panels",
                     c("ridge", "density", "boxplot", "blank"), "ridge"),
        param_choice("upper_facet", "Upper panels",
                     c("blank", "scatter", "contour", "ellipse"), "blank")
    ),
    summarise = function(x) {
        ev <- x[["discriminant_eigenvalues"]]
        if (is.null(ev)) return(list())
        list(
            figure("Discriminants", length(ev), "above the eigenvalue threshold"),
            figure("Leading eigenvalue", sprintf("%.1f", ev[1]),
                   sprintf("%s of total separating power", pct(ev[1] / sum(ev))))
        )
    },
    tables = function(x) {
        ev <- x[["discriminant_eigenvalues"]]
        list(
            "Eigenvalues" = data.frame(
                discriminant = paste0("DV", seq_along(ev)),
                eigenvalue = ev,
                proportion = ev / sum(ev)
            ),
            "Reference projection" = x[["ref_proj"]],
            "Query projection" = x[["query_proj"]]
        )
    },
    heavy = TRUE,
    plot_height = 620
))

add_entry(diag_entry(
    id = "calculateSIRSpace",
    fn = "calculateSIRSpace",
    label = "Sliced inverse regression space",
    category = "visualization",
    question = "Which directions of expression best predict cell type?",
    blurb = paste(
        "Sliced inverse regression finds a low-dimensional subspace that carries",
        "the relationship between expression and cell type."
    ),
    reading = paste(
        "SIR is a supervised alternative to PCA: it looks for directions along",
        "which the cell type label changes, rather than directions of largest",
        "variance. The loadings view names the genes driving each direction."
    ),
    data_map = map_pair(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_bool("multiple_cond_means", "Multiple conditional means", TRUE, advanced = TRUE),
        param_num("cumulative_variance_threshold", "Cumulative variance threshold",
                  0.7, 0.1, 0.99, 0.05, advanced = TRUE),
        param_int("n_neighbor", "Neighbours", 1, 1, 50, advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("plot_type", "View", c("scores", "loadings"), "scores"),
        param_choice("lower_facet", "Lower panels",
                     c("scatter", "contour", "ellipse", "blank"), "scatter"),
        param_choice("diagonal_facet", "Diagonal panels",
                     c("ridge", "density", "boxplot", "blank"), "ridge"),
        param_choice("upper_facet", "Upper panels",
                     c("blank", "scatter", "contour", "ellipse"), "blank"),
        param_int("n_top", "Top loadings", 10, 3, 50, advanced = TRUE)
    ),
    summarise = function(x) {
        pv <- x[["percent_var"]]
        if (is.null(pv)) return(list())
        list(figure("First two directions", pct(sum(pv[1:2]) / 100),
                    "of the cell-type-relevant variance"))
    },
    heavy = TRUE,
    plot_height = 620
))

add_entry(diag_entry(
    id = "plotGeneExpressionDimred",
    fn = "plotGeneExpressionDimred",
    label = "Gene expression on a reduction",
    category = "visualization",
    question = "Where is this gene expressed?",
    blurb = "Expression of one gene painted onto PCA, t-SNE or UMAP coordinates.",
    reading = paste(
        "Useful as a sanity check after a diagnostic flags a gene: if a marker",
        "is expressed in the wrong region of the embedding, the annotation of",
        "the cells in that region is suspect."
    ),
    data_map = map_single(cell_types = TRUE, assay = TRUE),
    params = list(
        param_gene("feature", "Gene"),
        param_reduction("method", "Reduction"),
        param_pcs(value = 1:3, max = 10,
                  help = "Used only when the reduction is PCA."),
        param_celltypes(),
        param_int("max_cells", "Max cells plotted", 2000, 100, 20000, advanced = TRUE)
    ),
    needs = "focus",
    returns = "plot"
))

# ===========================================================================
# Alignment and statistical comparison
# ===========================================================================

add_entry(diag_entry(
    id = "calculateWassersteinDistance",
    fn = "calculateWassersteinDistance",
    label = "Wasserstein distance",
    category = "alignment",
    tier = 1L,
    question = "Is the query further from the reference than the reference is from itself?",
    blurb = paste(
        "Compares the reference-to-query Wasserstein distance against a null",
        "built by repeatedly splitting the reference against itself."
    ),
    reading = paste(
        "The null distribution (reference vs. reference) is what distance you",
        "would see from sampling alone. The query distribution sitting to the",
        "right of it means the query is further away than sampling explains.",
        "The probability of superiority is the chance that a random",
        "reference-query distance exceeds a random reference-reference one:",
        "0.5 means indistinguishable, values near 1 mean clearly separated."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20),
        param_celltypes(),
        param_int("n_resamples", "Resamples", 300, 50, 2000,
                  help = "More resamples give a smoother null and a slower run."),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(param_celltypes("plot_cell_types", "Cell types to draw")),
    summarise = function(x) {
        ps <- x[["probability_superiority"]]
        if (is.null(ps) || !length(ps)) return(list())
        worst <- names(ps)[which.max(ps)]
        list(
            figure("Cell types compared", length(ps)),
            figure("Least aligned", worst,
                   sprintf("probability of superiority %.2f", max(ps)),
                   tone = if (max(ps) > 0.9) "warning" else "default")
        )
    },
    tables = function(x) {
        ps <- x[["probability_superiority"]]
        rr <- x[["ref_ref_dist"]]
        rq <- x[["ref_query_dist"]]
        list("Per cell type" = data.frame(
            cell_type = names(ps),
            probability_superiority = as.numeric(ps),
            median_ref_ref = vapply(names(ps), function(k) stats::median(rr[[k]]), numeric(1)),
            median_ref_query = vapply(names(ps), function(k) stats::median(rq[[k]]), numeric(1)),
            row.names = NULL
        ))
    },
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateGraphIntegration",
    fn = "calculateGraphIntegration",
    label = "Graph integration",
    category = "alignment",
    tier = 1L,
    question = "Do query cells share neighbourhoods with reference cells of the same type?",
    blurb = paste(
        "Builds a shared nearest-neighbour graph over both datasets, finds",
        "communities, and reports where annotations disagree with local structure."
    ),
    reading = paste(
        "Three distinct problems are reported. A community that is almost",
        "entirely query cells means a population the reference does not contain.",
        "Cross-type mixing means query cells of one type sitting among reference",
        "cells of another. Local inconsistency is per-cell: a cell whose",
        "reference neighbours mostly carry a different label, with a confidence",
        "score for the suggested relabelling. High modularity means the graph",
        "split cleanly, so the communities are worth trusting."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:10, max = 30),
        param_celltypes(),
        param_int("k_neighbors", "Neighbours (k)", 30, 5, 200,
                  help = "Larger k gives smoother, coarser communities."),
        param_num("resolution", "Community resolution", 0.1, 0.01, 2, 0.01,
                  help = "Higher values split the graph into more, smaller communities."),
        param_int("min_cells_per_community", "Min cells per community", 10, 2, 500, advanced = TRUE),
        param_int("min_cells_per_celltype", "Min cells per cell type", 20, 2, 500, advanced = TRUE),
        param_num("high_query_prop_threshold", "Query-dominated threshold", 0.9, 0.5, 1, 0.01, advanced = TRUE),
        param_num("cross_type_threshold", "Cross-type threshold", 0.15, 0.01, 1, 0.01, advanced = TRUE),
        param_num("local_consistency_threshold", "Local consistency threshold", 0.6, 0.1, 1, 0.05, advanced = TRUE),
        param_num("local_confidence_threshold", "Local confidence threshold", 0.2, 0, 1, 0.05, advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("plot_type", "View",
                     c("community_network", "cell_network", "community_data",
                       "summary", "local_issues", "annotation_issues"),
                     "summary"),
        param_choice("color_by", "Colour by", c("cell_type", "community_type"), "cell_type"),
        param_int("max_nodes", "Max nodes drawn", 2000, 100, 20000, advanced = TRUE),
        param_bool("exclude_reference_only", "Hide reference-only communities", FALSE, advanced = TRUE)
    ),
    summarise = function(x) {
        m <- x[["overall_metrics"]]
        if (is.null(m)) return(list())
        list(
            figure("Communities", m[["total_communities"]],
                   sprintf("modularity %.2f", m[["modularity"]])),
            figure("Query-only communities", m[["high_query_prop_communities"]],
                   "populations missing from the reference",
                   tone = if ((m[["high_query_prop_communities"]] %||% 0) > 0) "warning" else "default"),
            figure("Locally inconsistent cells", m[["total_locally_inconsistent_cells"]],
                   sprintf("%s of query cells", pct(m[["mean_local_inconsistency_rate"]])),
                   tone = if ((m[["mean_local_inconsistency_rate"]] %||% 0) > 0.05) "warning" else "default")
        )
    },
    tables = function(x) {
        list(
            "Per cell type" = x[["local_inconsistency_summary"]],
            "Annotation consistency" = x[["annotation_consistency"]],
            "Communities" = x[["community_composition"]],
            "Cells to re-check" = x[["local_annotation_inconsistencies"]],
            "Cross-type mixing" = x[["cross_type_mixing"]]
        )
    },
    heavy = TRUE,
    plot_height = 600
))

add_entry(diag_entry(
    id = "comparePCASubspace",
    fn = "comparePCASubspace",
    label = "PCA subspace comparison",
    category = "alignment",
    tier = 1L,
    question = "Do reference and query PCA span the same subspace for this cell type?",
    blurb = paste(
        "Compares the principal subspaces fitted separately on reference and",
        "query cells of a single cell type."
    ),
    reading = paste(
        "This asks a stricter question than whether cells overlap: whether the",
        "axes of variation themselves agree. Disagreement means the two datasets",
        "vary along different biological or technical directions, which usually",
        "precedes a failed label transfer."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        pc_subset = "pc_subset"
    ),
    params = list(
        param_pcs(value = 1:5, max = 20),
        param_celltypes("cell_types", "Cell type", multiple = FALSE, allow_empty = FALSE,
                        help = "This diagnostic compares one cell type at a time."),
        param_int("n_top_vars", "Top genes per component", 50, 5, 500, advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        ct <- args$cell_types
        if (is.null(ct) || !length(ct)) {
            stop("Choose exactly one cell type: this diagnostic compares a single population.")
        }
        args$cell_types <- NULL
        args$query_data <- subset_to_types(args$query_data, ctx$query_col, ct[1])
        args$reference_data <- subset_to_types(args$reference_data, ctx$ref_col, ct[1])
        args
    },
    repro_prelude = single_type_prelude,
    available = needs_matched_rotation,
    heavy = TRUE
))

add_entry(diag_entry(
    id = "comparePCA",
    fn = "comparePCA",
    label = "PCA loading comparison",
    category = "alignment",
    question = "Do the individual components match, one by one?",
    blurb = "Similarity between each reference component and each query component, for a single cell type.",
    reading = paste(
        "Read the heatmap along the diagonal: a strong diagonal means PC1",
        "matches PC1, PC2 matches PC2 and so on. Strong off-diagonal entries",
        "mean the same variation is present but ordered differently, which is",
        "far less alarming than no correspondence at all."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        pc_subset = "pc_subset"
    ),
    params = list(
        param_pcs(value = 1:5, max = 15),
        param_celltypes("cell_types", "Cell type", multiple = FALSE, allow_empty = FALSE,
                        help = "This diagnostic compares one cell type at a time."),
        param_int("n_top_vars", "Top genes per component", 50, 5, 500, advanced = TRUE),
        param_choice("metric", "Metric", c("cosine", "correlation"), "cosine"),
        param_choice("correlation_method", "Correlation method",
                     c("spearman", "pearson"), "spearman", advanced = TRUE),
        param_int("n_permutations", "Permutations for significance", 0, 0, 2000,
                  help = "Zero skips the permutation test.", advanced = TRUE)
    ),
    plot_params = list(
        param_bool("show_values", "Show values", TRUE),
        param_bool("show_significance", "Show significance", TRUE),
        param_num("significance_threshold", "Significance threshold", 0.05, 0.001, 0.2, 0.005, advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        ct <- args$cell_types
        if (is.null(ct) || !length(ct)) {
            stop("Choose exactly one cell type: this diagnostic compares a single population.")
        }
        args$cell_types <- NULL
        args$query_data <- subset_to_types(args$query_data, ctx$query_col, ct[1])
        args$reference_data <- subset_to_types(args$reference_data, ctx$ref_col, ct[1])
        args
    },
    repro_prelude = single_type_prelude,
    available = needs_matched_rotation,
    heavy = TRUE
))

add_entry(diag_entry(
    id = "regressPC",
    fn = "regressPC",
    label = "Variance attribution",
    category = "alignment",
    tier = 1L,
    question = "Is the variation in each component driven by biology or by batch?",
    blurb = paste(
        "Regresses each principal component on cell type, dataset of origin and",
        "their interaction, and reports how much variance each explains."
    ),
    reading = paste(
        "This is the cleanest single readout for a batch effect. If cell type",
        "explains most of the variance in the leading components, the two",
        "datasets are compatible. If the dataset term explains an appreciable",
        "share, the query differs from the reference in a way that has nothing",
        "to do with biology, and every downstream diagnostic will inherit that."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        cell_types = "cell_types", pc_subset = "pc_subset", assay_name = "assay"
    ),
    params = list(
        param_pcs(value = 1:10, max = 25),
        param_celltypes(),
        param_choice("adjust_method", "p-value adjustment",
                     c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                     "BH", advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("plot_type", "View",
                     c("r_squared", "variance_contribution", "coefficient_heatmap"),
                     "variance_contribution"),
        param_num("alpha", "Alpha", 0.05, 0.001, 0.2, 0.005, advanced = TRUE)
    ),
    summarise = function(x) {
        comp <- x[["var_contributions_components"]]
        if (is.null(comp)) return(list())
        tot <- x[["total_variance_explained"]] %||% NA_real_
        ds <- sum(comp[["dataset"]], na.rm = TRUE)
        ct <- sum(comp[["cell_type"]], na.rm = TRUE)
        list(
            figure("Variance explained", sprintf("%.1f%%", tot), "across the components shown"),
            figure("Attributed to cell type", sprintf("%.1f%%", ct), "biological signal"),
            figure("Attributed to dataset", sprintf("%.2f%%", ds), "batch signal",
                   tone = if (is.finite(ds) && is.finite(ct) && ct > 0 && ds > 0.1 * ct) "warning" else "success")
        )
    },
    tables = function(x) {
        comp <- x[["var_contributions_components"]]
        pcs <- names(x[["r_squared"]])
        list("Per component" = data.frame(
            component = pcs,
            r_squared = as.numeric(x[["r_squared"]]),
            var_cell_type = as.numeric(comp[["cell_type"]]),
            var_dataset = as.numeric(comp[["dataset"]]),
            var_interaction = as.numeric(comp[["interaction"]]),
            row.names = NULL
        ))
    }
))

add_entry(diag_entry(
    id = "calculateAveragePairwiseCorrelation",
    fn = "calculateAveragePairwiseCorrelation",
    label = "Average pairwise correlation",
    category = "alignment",
    question = "How correlated is each query cell type with each reference cell type?",
    blurb = "Mean correlation between cells of each query type and each reference type.",
    reading = paste(
        "Read the diagonal: a query cell type should correlate most strongly",
        "with the reference type of the same name. An off-diagonal maximum names",
        "the label the cells may actually deserve."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:10, max = 25),
        param_celltypes(),
        param_choice("correlation_method", "Method", c("spearman", "pearson"), "spearman"),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    heavy = TRUE
))

add_entry(diag_entry(
    id = "plotPairwiseDistancesDensity",
    fn = "plotPairwiseDistancesDensity",
    label = "Pairwise distance density",
    category = "alignment",
    question = "How do within-reference distances compare with reference-to-query distances?",
    blurb = "Density of pairwise distances or correlations for one cell type.",
    reading = paste(
        "Two overlapping densities mean the query is indistinguishable from the",
        "reference for this cell type. A right-shifted query density means query",
        "cells are systematically further from reference cells than reference",
        "cells are from each other."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        pc_subset = "pc_subset", assay_name = "assay"
    ),
    params = list(
        param_celltypes("cell_type", "Cell type", multiple = FALSE, allow_empty = FALSE,
                        help = "One cell type at a time."),
        param_pcs(value = 1:5, max = 20),
        param_choice("distance_metric", "Metric", c("correlation", "euclidean"), "correlation"),
        param_choice("correlation_method", "Correlation method",
                     c("spearman", "pearson"), "spearman", advanced = TRUE),
        param_num("bandwidth", "Bandwidth", 0.25, 0.01, 2, 0.01, advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        if (is.null(args$cell_type) || !length(args$cell_type)) {
            stop("Choose a cell type.")
        }
        args$cell_type <- args$cell_type[1]
        args
    },
    returns = "plot",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateCramerPValue",
    fn = "calculateCramerPValue",
    label = "Cramer test",
    category = "alignment",
    question = "Are the two multivariate distributions different, non-parametrically?",
    blurb = "Cramer two-sample test per cell type, on the PCA projections.",
    reading = paste(
        "A small p-value means the reference and query point clouds for that",
        "cell type are unlikely to be draws from the same distribution. With",
        "thousands of cells even trivial differences reach significance, so read",
        "the p-value alongside an effect size such as the Wasserstein distance."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20), param_celltypes(),
        param_int("max_cells_query", "Max query cells", 5000, 100, 20000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    tables = function(x) list("p-values" = data.frame(
        cell_type = names(x), p_value = as.numeric(x), row.names = NULL)),
    returns = "value",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateHotellingPValue",
    fn = "calculateHotellingPValue",
    label = "Hotelling T² test",
    category = "alignment",
    question = "Do the two groups have the same multivariate mean?",
    blurb = "Permutation-based Hotelling T-squared test per cell type.",
    reading = paste(
        "This tests the centre of each cloud only. A significant result with",
        "well-overlapping clouds means a small but consistent shift; a",
        "non-significant result does not rule out a difference in shape or",
        "spread, which the Cramer test is better placed to detect."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20), param_celltypes(),
        param_int("n_permutation", "Permutations", 500, 50, 5000),
        param_int("max_cells_query", "Max query cells", 5000, 100, 20000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    tables = function(x) list("p-values" = data.frame(
        cell_type = names(x), p_value = as.numeric(x), row.names = NULL)),
    returns = "value",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateMMDPValue",
    fn = "calculateMMDPValue",
    label = "Maximum mean discrepancy",
    category = "alignment",
    question = "Do the two distributions differ in any respect a kernel can see?",
    blurb = "Kernel two-sample test (MMD) per cell type.",
    reading = paste(
        "MMD is sensitive to differences in mean, spread and shape at once,",
        "which makes it a good catch-all. The permutation count sets the finest",
        "p-value obtainable: 100 permutations cannot produce a p below 0.01."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 20), param_celltypes(),
        param_int("n_permutation", "Permutations", 100, 20, 2000),
        param_choice("kernel_type", "Kernel", c("gaussian"), "gaussian", advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 20000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    tables = function(x) list("p-values" = data.frame(
        cell_type = names(x), p_value = as.numeric(x), row.names = NULL)),
    returns = "value",
    heavy = TRUE
))

# ===========================================================================
# Anomaly detection & cell distances
# ===========================================================================

add_entry(diag_entry(
    id = "detectAnomaly",
    fn = "detectAnomaly",
    label = "Anomaly detection",
    category = "anomaly",
    tier = 1L,
    question = "Which individual cells do not look like their assigned cell type?",
    blurb = paste(
        "Fits an isolation forest to each reference cell type and scores every",
        "query cell against the forest for the type it was assigned."
    ),
    reading = paste(
        "Each cell gets a score between 0 and 1; higher means easier to isolate,",
        "and therefore more unusual relative to the reference cells of that",
        "type. The threshold is derived from the reference scores themselves",
        "(median plus a multiple of the MAD), so it adapts to how tight each",
        "population is. A high flagged rate in one cell type and not others",
        "points at that annotation. A high rate everywhere usually means a batch",
        "effect rather than a labelling problem, so check variance attribution."
    ),
    data_map = map_ref_first(),
    params = list(
        param_pcs(value = 1:5, max = 25),
        param_celltypes(),
        param_int("n_tree", "Trees", 500, 50, 2000,
                  help = "More trees give a more stable score."),
        param_choice("threshold_method", "Threshold", c("MAD", "absolute"), "MAD",
                     help = "MAD adapts to each cell type; absolute uses one fixed cut."),
        param_num("mad_multiplier", "MAD multiplier", 2, 0.5, 6, 0.1,
                  help = "Higher is more conservative: fewer cells flagged."),
        param_num("anomaly_threshold", "Absolute threshold", 0.5, 0.05, 0.95, 0.05,
                  help = "Used only when the threshold method is 'absolute'."),
        param_int("n_hvgs", "Highly variable genes", 100, 20, 2000,
                  help = "Used only when no principal components are supplied.",
                  advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("cell_type", "Cell type to draw", c("Combined"), "Combined"),
        param_choice("data_type", "Dataset", c("query", "reference", "both"), "query"),
        param_choice("upper_facet", "Upper panels", c("blank", "contour", "ellipse"), "blank"),
        param_choice("diagonal_facet", "Diagonal panels",
                     c("density", "ridge", "boxplot", "blank"), "density")
    ),
    summarise = function(x) {
        tb <- anomaly_table(x, "query_anomaly", "query_anomaly_scores")
        cmb <- x[["Combined"]]
        total <- sum(tb$flagged, na.rm = TRUE)
        n <- sum(tb$query_cells, na.rm = TRUE)
        w <- worst_type(tb)
        out <- list(
            figure("Query cells flagged", n_fmt(total),
                   sprintf("%s of %s cells", pct(if (n) total / n else NA), n_fmt(n)),
                   tone = if (n && total / n > 0.1) "warning" else "default")
        )
        if (!is.null(w)) {
            out <- c(out, list(figure(
                "Most affected type", w$cell_type,
                sprintf("%s flagged", pct(w$flagged_pct / 100)),
                tone = if (w$flagged_pct > 20) "danger" else "default"
            )))
        }
        if (!is.null(cmb)) {
            out <- c(out, list(figure(
                "Threshold", sprintf("%.3f", cmb[["applied_threshold"]]),
                "isolation score cut for the combined fit"
            )))
        }
        out
    },
    tables = function(x) {
        tb <- anomaly_table(x, "query_anomaly", "query_anomaly_scores")
        cmb <- x[["Combined"]]
        cells <- NULL
        if (!is.null(cmb) && length(cmb[["query_anomaly_scores"]])) {
            cells <- data.frame(
                cell = names(cmb[["query_anomaly_scores"]]),
                anomaly_score = as.numeric(cmb[["query_anomaly_scores"]]),
                flagged = as.logical(cmb[["query_anomaly"]]),
                row.names = NULL
            )
            cells <- cells[order(-cells$anomaly_score), , drop = FALSE]
        }
        not_null(list("Per cell type" = tb, "Per cell (combined fit)" = cells))
    },
    heavy = TRUE,
    plot_height = 620
))

add_entry(diag_entry(
    id = "calculateReconstructionError",
    fn = "calculateReconstructionError",
    label = "Reconstruction error",
    category = "anomaly",
    tier = 1L,
    question = "Which cells cannot be rebuilt from the reference's own components?",
    blurb = paste(
        "Fits PCA on each reference cell type, projects query cells into it and",
        "back out, and measures what the round trip loses."
    ),
    reading = paste(
        "A cell that the reference subspace reconstructs poorly carries",
        "variation the reference does not contain. This is an independent second",
        "opinion on the isolation forest: the two methods fail in different",
        "ways, so a cell flagged by both is a much stronger candidate than a",
        "cell flagged by either alone."
    ),
    data_map = map_ref_first(),
    params = list(
        param_pcs(value = 1:5, max = 25),
        param_celltypes(),
        param_int("n_hvgs", "Highly variable genes", 100, 20, 2000),
        param_num("mad_multiplier", "MAD multiplier", 2, 0.5, 6, 0.1),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_choice("cell_type", "Cell type to draw", c("Combined"), "Combined"),
        param_choice("data_type", "Dataset", c("both", "query", "reference"), "both"),
        param_choice("plot_type", "Shape", c("violin", "boxplot", "ridge", "heatmap"), "violin")
    ),
    summarise = function(x) {
        tb <- anomaly_table(x, "query_anomaly", "query_reconstruction_errors")
        total <- sum(tb$flagged, na.rm = TRUE)
        n <- sum(tb$query_cells, na.rm = TRUE)
        w <- worst_type(tb)
        out <- list(figure("Query cells flagged", n_fmt(total),
                           sprintf("%s of %s cells", pct(if (n) total / n else NA), n_fmt(n)),
                           tone = if (n && total / n > 0.1) "warning" else "default"))
        if (!is.null(w)) {
            out <- c(out, list(figure("Most affected type", w$cell_type,
                                      sprintf("%s flagged", pct(w$flagged_pct / 100)))))
        }
        out
    },
    tables = function(x) list(
        "Per cell type" = anomaly_table(x, "query_anomaly", "query_reconstruction_errors")
    ),
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateCellDistances",
    fn = "calculateCellDistances",
    label = "Cell distances",
    category = "anomaly",
    question = "How far is each query cell from the reference population?",
    blurb = "Distances between query cells and reference cells of the same type, in PCA space.",
    reading = paste(
        "The reference-to-reference distances describe how spread out a healthy",
        "population is. A query cell whose distances sit well outside that range",
        "is far from every reference cell of its assigned type, not merely from",
        "the centroid."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 25), param_celltypes(),
        param_int("max_cells_query", "Max query cells", 5000, 100, 20000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    plot_fn = function(result, args, ctx) {
        ct <- args$ref_cell_type %|% names(result)[1]
        cells <- args$cell_names
        if (is.null(cells) || !length(cells)) {
            # Default to the query cells furthest from this reference population.
            d <- result[[ct]][["query_to_ref_distances"]]
            if (is.null(d)) stop("No distances for this cell type.")
            cells <- rownames(d)[utils::head(order(-rowMeans(d)), 4)]
        }
        graphics::plot(result, ref_cell_type = ct, cell_names = cells)
    },
    plot_params = list(
        param_choice("ref_cell_type", "Reference population", character(0)),
        param_cellnames("cell_names", "Query cells to draw", max_select = 8,
                        help = "Left empty, the four most distant query cells are shown.")
    ),
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateCellDistancesSimilarity",
    fn = "calculateCellDistancesSimilarity",
    label = "Cell distance similarity",
    category = "anomaly",
    question = "Do specific cells have the same distance profile as the reference?",
    blurb = paste(
        "Bhattacharyya coefficient and Hellinger distance between the distance",
        "distribution of chosen query cells and that of the reference."
    ),
    reading = paste(
        "A Bhattacharyya coefficient near 1 means the cell's distance profile is",
        "indistinguishable from a typical reference cell's. Values near 0 mean",
        "the cell sits in a different part of the space entirely."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        cell_types = "cell_types", pc_subset = "pc_subset", assay_name = "assay"
    ),
    params = list(
        param_cellnames("cell_names_query", "Query cells", max_select = 12,
                        help = "Pick the cells to profile. Anomaly detection is a good way to choose them."),
        param_pcs(value = 1:5, max = 25),
        param_celltypes(),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 20000, advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        if (is.null(args$cell_names_query) || !length(args$cell_names_query)) {
            stop("Choose at least one query cell.")
        }
        args
    },
    # Returns two data frames and has no plot method: the tables are the output.
    returns = "value",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "calculateCellSimilarityPCA",
    fn = "calculateCellSimilarityPCA",
    label = "Cell / component similarity",
    category = "anomaly",
    question = "Which components does an unusual cell load on?",
    blurb = "Cosine similarity between chosen cells and the top-loading genes of each component.",
    reading = paste(
        "Once a cell is flagged, this says which axes of variation make it",
        "unusual. Pair it with the gene shift diagnostic, which names the genes",
        "driving those same components."
    ),
    data_map = map_single(arg = "sce_object", col = NULL, cell_types = FALSE, assay = TRUE),
    params = list(
        param_cellnames("cell_names", "Cells", max_select = 12,
                        help = "Pick a handful of cells to profile."),
        param_pcs(value = 1:5, max = 25),
        param_int("n_top_vars", "Top genes per component", 50, 5, 500, advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        args$cell_type_col <- NULL
        if (is.null(args$cell_names) || !length(args$cell_names)) {
            stop("Choose at least one cell.")
        }
        args
    },
    needs = "focus"
))

# ===========================================================================
# Marker genes
# ===========================================================================

add_entry(diag_entry(
    id = "calculateGeneShifts",
    fn = "calculateGeneShifts",
    label = "Gene expression shifts",
    category = "markers",
    tier = 1L,
    question = "Which genes explain the difference between query and reference?",
    blurb = paste(
        "Tests the genes that load most heavily on each principal component for",
        "a shift in expression between reference and query."
    ),
    reading = paste(
        "This is the step that turns a geometric observation into biology. The",
        "components tell you the query has moved; the genes loading on those",
        "components tell you why. Read the adjusted p-value, not the raw one,",
        "and look for a coherent set of genes rather than a single hit: a whole",
        "pathway shifting is a real state difference, one gene moving is usually",
        "noise or a technical artefact. Turn on anomaly comparison to contrast",
        "the flagged cells against the rest of the query rather than against the",
        "reference."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:5, max = 15),
        param_celltypes(),
        param_int("n_top_loadings", "Top loading genes per component", 50, 5, 500),
        param_num("p_value_threshold", "p-value threshold", 0.05, 0.001, 0.2, 0.005),
        param_choice("adjust_method", "p-value adjustment",
                     c("fdr", "BH", "holm", "bonferroni", "BY", "none"), "fdr"),
        param_bool("detect_anomalies", "Run anomaly detection first", FALSE,
                   help = "Needed for the anomaly comparison below."),
        param_bool("anomaly_comparison", "Compare anomalous vs. other query cells", FALSE),
        param_choice("threshold_method", "Anomaly threshold", c("MAD", "absolute"), "MAD", advanced = TRUE),
        param_num("mad_multiplier", "MAD multiplier", 2, 0.5, 6, 0.1, advanced = TRUE),
        param_num("anomaly_threshold", "Absolute anomaly threshold", 0.5, 0.05, 0.95, 0.05, advanced = TRUE),
        param_int("n_tree", "Trees", 500, 50, 2000, advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(
        param_celltypes("cell_type", "Cell type", multiple = FALSE, allow_empty = FALSE),
        param_choice("plot_type", "View", c("heatmap", "barplot", "boxplot"), "heatmap"),
        param_choice("plot_by", "Rank genes by", c("p_adjusted", "top_loading"), "p_adjusted"),
        param_int("n_genes", "Genes shown", 10, 3, 60),
        param_bool("show_anomalies", "Mark anomalous cells", FALSE, advanced = TRUE),
        param_bool("pseudo_bulk", "Pseudo-bulk by cell type", FALSE, advanced = TRUE),
        param_bool("cluster_cols", "Cluster columns", FALSE, advanced = TRUE)
    ),
    summarise = function(x) {
        pcs <- grep("^PC[0-9]+$", names(x), value = TRUE)
        if (!length(pcs)) return(list())
        all_df <- do.call(rbind, lapply(pcs, function(p) x[[p]]))
        sig <- all_df[isTRUE_vec(all_df$significant), , drop = FALSE]
        top <- if (nrow(sig)) sig[order(sig$p_adjusted), ][seq_len(min(3, nrow(sig))), ] else NULL
        out <- list(
            figure("Genes tested", n_fmt(nrow(all_df)),
                   sprintf("across %s", paste(pcs, collapse = ", "))),
            figure("Significantly shifted", n_fmt(nrow(sig)),
                   pct(if (nrow(all_df)) nrow(sig) / nrow(all_df) else NA),
                   tone = if (nrow(all_df) && nrow(sig) / nrow(all_df) > 0.3) "warning" else "default")
        )
        if (!is.null(top)) {
            out <- c(out, list(figure("Strongest shifts", paste(unique(top$gene), collapse = ", "),
                                      "lowest adjusted p-values")))
        }
        out
    },
    tables = function(x) {
        pcs <- grep("^PC[0-9]+$", names(x), value = TRUE)
        tabs <- lapply(pcs, function(p) x[[p]])
        names(tabs) <- pcs
        c(tabs, list("Variance by cell type" = x[["cell_type_variance"]]))
    },
    heavy = TRUE,
    plot_height = 560
))

add_entry(diag_entry(
    id = "compareMarkers",
    fn = "compareMarkers",
    label = "Marker comparison",
    category = "markers",
    tier = 1L,
    question = "Do the same genes mark the same cell types in both datasets?",
    blurb = paste(
        "Derives marker genes for each cell type separately in reference and",
        "query, and measures how far the two sets agree."
    ),
    reading = paste(
        "Marker overlap is the fraction of top markers shared. A cell type with",
        "low overlap is defined by different genes in the two datasets, which",
        "means the label is being applied to different biology even if the name",
        "matches. The quality score summarises overlap and expression",
        "consistency together."
    ),
    data_map = map_pair(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_int("n_markers", "Markers per cell type", 50, 5, 500),
        param_int("min_cells", "Minimum cells per type", 10, 3, 200),
        param_choice("anomaly_filter", "Restrict to",
                     c("none", "anomalous_only", "non_anomalous_only"), "none",
                     help = "Recompute markers using only flagged, or only unflagged, query cells."),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    plot_params = list(param_celltypes("cell_types", "Cell types to draw")),
    summarise = function(x) {
        ov <- x[["marker_overlap"]]
        if (is.null(ov) || !length(ov)) return(list())
        worst <- names(ov)[which.min(ov)]
        list(
            figure("Median marker overlap", sprintf("%.2f", stats::median(ov)),
                   "1.0 means identical marker sets"),
            figure("Weakest cell type", worst, sprintf("overlap %.2f", min(ov)),
                   tone = if (min(ov) < 0.5) "warning" else "default")
        )
    },
    tables = function(x) {
        list("Per cell type" = data.frame(
            cell_type = x[["common_cell_types"]],
            marker_overlap = as.numeric(x[["marker_overlap"]]),
            expression_consistency = as.numeric(x[["expression_consistency"]]),
            quality = as.character(x[["quality_scores"]]),
            n_cells_query = as.numeric(x[["n_cells_query"]]),
            n_cells_reference = as.numeric(x[["n_cells_ref"]]),
            row.names = NULL
        ))
    },
    heavy = TRUE
))

add_entry(diag_entry(
    id = "plotMarkerExpression",
    fn = "plotMarkerExpression",
    label = "Marker expression",
    category = "markers",
    tier = 1L,
    question = "How does one gene's expression compare between the datasets?",
    blurb = "Density of a single gene's expression in reference and query, for one cell type.",
    reading = paste(
        "A shifted or differently shaped density means the gene behaves",
        "differently in the query. Normalisation matters here: z-score removes",
        "scale differences and shows distributional change, while 'none' keeps",
        "raw values and will show depth differences as well as biology."
    ),
    data_map = list(
        query_data = "query", reference_data = "ref",
        query_cell_type_col = "query_col", ref_cell_type_col = "ref_col",
        assay_name = "assay"
    ),
    params = list(
        param_celltypes("cell_type", "Cell type", multiple = FALSE, allow_empty = FALSE),
        param_gene("gene_name", "Gene"),
        param_choice("normalization", "Normalisation",
                     c("z_score", "min_max", "rank", "none"), "z_score")
    ),
    prepare = function(ctx, args) {
        if (is.null(args$cell_type) || !length(args$cell_type)) stop("Choose a cell type.")
        if (is.null(args$gene_name)) stop("Choose a gene.")
        args$cell_type <- args$cell_type[1]
        args
    },
    returns = "plot"
))

add_entry(diag_entry(
    id = "calculateHVGOverlap",
    fn = "calculateHVGOverlap",
    label = "Highly variable gene overlap",
    category = "markers",
    question = "Do the two datasets find the same genes variable?",
    blurb = "Overlap coefficient between the highly variable genes of each dataset.",
    reading = paste(
        "The app selects the top highly variable genes in each dataset with",
        "scran, then reports the overlap coefficient. A low value means the two",
        "datasets disagree about which genes carry signal at all, which",
        "undermines any shared embedding built from them."
    ),
    data_map = list(),
    params = list(
        param_int("n_hvgs", "Highly variable genes per dataset", 500, 50, 5000)
    ),
    prepare = function(ctx, args) {
        n <- args$n_hvgs %||% 500
        list(
            reference_genes = scran::getTopHVGs(ctx$ref, n = min(n, nrow(ctx$ref))),
            query_genes = scran::getTopHVGs(ctx$query, n = min(n, nrow(ctx$query)))
        )
    },
    repro_prelude = function(ctx, args, param_values, ref_var, query_var) {
        n <- param_values$n_hvgs %||% 500
        c(
            "# Highly variable genes are selected with scran before comparing.",
            sprintf("reference_genes <- scran::getTopHVGs(%s, n = %d)", ref_var, n),
            sprintf("query_genes <- scran::getTopHVGs(%s, n = %d)", query_var, n)
        )
    },
    repro_names = list(reference_genes = "reference_genes", query_genes = "query_genes"),
    summarise = function(x) list(
        figure("Overlap coefficient", sprintf("%.2f", as.numeric(x)),
               "1.0 means the same genes are variable in both",
               tone = if (as.numeric(x) < 0.5) "warning" else "success")
    ),
    returns = "value"
))

add_entry(diag_entry(
    id = "calculateVarImpOverlap",
    fn = "calculateVarImpOverlap",
    label = "Gene importance overlap",
    category = "markers",
    question = "Do the same genes distinguish each pair of cell types in both datasets?",
    blurb = "Random-forest gene importance for every pair of cell types, compared across datasets.",
    reading = paste(
        "For each pair of cell types, a forest is trained to tell them apart in",
        "each dataset, and the top genes are compared. A low overlap for a",
        "particular pair says the boundary between those two types is drawn",
        "using different genes in the query, which is exactly the situation in",
        "which a label transfer between them goes wrong."
    ),
    data_map = map_ref_first(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_int("n_tree", "Trees", 500, 50, 2000),
        param_int("n_top", "Top genes compared", 50, 5, 500),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    summarise = function(x) {
        cmp <- x[["var_imp_comparison"]]
        if (is.null(cmp) || !length(cmp)) return(list())
        list(
            figure("Median overlap", sprintf("%.2f", stats::median(cmp)), "across cell type pairs"),
            figure("Weakest pair", names(cmp)[which.min(cmp)],
                   sprintf("overlap %.2f", min(cmp)),
                   tone = if (min(cmp) < 0.4) "warning" else "default")
        )
    },
    tables = function(x) {
        cmp <- x[["var_imp_comparison"]]
        list("Cell type pairs" = data.frame(
            pair = names(cmp), overlap = as.numeric(cmp), row.names = NULL))
    },
    returns = "value",
    heavy = TRUE
))

# ===========================================================================
# QC & annotation scores
# ===========================================================================

add_entry(diag_entry(
    id = "plotQCvsAnnotation",
    fn = "plotQCvsAnnotation",
    label = "QC vs. annotation score",
    category = "qc",
    tier = 1L,
    question = "Are low-confidence annotations just low-quality cells?",
    blurb = "Scatter of a quality control metric against annotation confidence, by cell type.",
    reading = paste(
        "A downward trend means the annotator is least confident exactly where",
        "the data is worst, so those calls are a data quality problem rather",
        "than a biological one. No trend, with low-confidence cells spread",
        "across all quality levels, points instead at a genuine mismatch",
        "between reference and query."
    ),
    data_map = map_single(assay = FALSE),
    params = list(
        param_numeric_col("qc_col", "QC metric", prefer = c("percent_mito", "total", "sum", "detected")),
        param_numeric_col("score_col", "Annotation score", prefer = c("annotation_scores", "scores")),
        param_celltypes(),
        param_int("max_cells", "Max cells plotted", 5000, 100, 50000, advanced = TRUE)
    ),
    available = needs_numeric_cols(2, "a QC metric and an annotation score"),
    needs = "focus",
    returns = "plot"
))

add_entry(diag_entry(
    id = "histQCvsAnnotation",
    fn = "histQCvsAnnotation",
    label = "QC and score histograms",
    category = "qc",
    question = "How are quality and confidence distributed?",
    blurb = "Side-by-side histograms of a QC metric and an annotation score.",
    reading = paste(
        "A bimodal confidence histogram usually means two populations: cells the",
        "reference covers well and cells it does not. Compare where the second",
        "mode sits against the QC histogram to tell quality from biology."
    ),
    data_map = map_single(assay = FALSE),
    params = list(
        param_numeric_col("qc_col", "QC metric", prefer = c("percent_mito", "total", "sum", "detected")),
        param_numeric_col("score_col", "Annotation score", prefer = c("annotation_scores", "scores")),
        param_celltypes(),
        param_int("max_cells", "Max cells", 5000, 100, 50000, advanced = TRUE)
    ),
    available = needs_numeric_cols(2, "a QC metric and an annotation score"),
    needs = "focus",
    returns = "plot"
))

add_entry(diag_entry(
    id = "plotGeneSetScores",
    fn = "plotGeneSetScores",
    label = "Gene set scores",
    category = "qc",
    question = "Where do high gene-set scores sit in the embedding?",
    blurb = "A gene set or pathway score painted onto a dimensionality reduction.",
    reading = paste(
        "Score concentrated in one region that crosses annotation boundaries",
        "suggests a cell state cutting across the cell type labels, which is a",
        "common reason for confident but wrong annotations."
    ),
    data_map = map_single(assay = FALSE),
    params = list(
        param_numeric_col("score_col", "Score column", prefer = c("gene_set_scores", "annotation_scores")),
        param_reduction("method", "Reduction"),
        param_pcs(value = 1:5, max = 10, help = "Used only when the reduction is PCA."),
        param_celltypes(),
        param_int("max_cells", "Max cells plotted", 2000, 100, 20000, advanced = TRUE)
    ),
    available = needs_numeric_cols(1, "a numeric score column"),
    needs = "focus",
    returns = "plot"
))

add_entry(diag_entry(
    id = "calculateCategorizationEntropy",
    fn = "calculateCategorizationEntropy",
    label = "Categorization entropy",
    category = "qc",
    question = "How decisive is the annotator, cell by cell?",
    blurb = paste(
        "Entropy of the per-cell score vector across candidate cell types. Needs",
        "at least two numeric score columns in the chosen object."
    ),
    reading = paste(
        "Low entropy means the annotator strongly preferred one cell type; high",
        "entropy means it was close to indifferent. A bimodal histogram is",
        "healthy. A mass of high-entropy cells means a large set of calls that",
        "were effectively coin flips, whatever label they were finally given."
    ),
    data_map = list(),
    params = list(
        param_choice("score_columns", "Score columns", character(0), multiple = TRUE,
                     help = "Pick two or more numeric columns holding per-cell-type scores."),
        param_bool("inverseNormalTransformationform", "Inverse normal transform first", FALSE,
                   advanced = TRUE)
    ),
    prepare = function(ctx, args) {
        cols <- args$score_columns
        if (is.null(cols) || length(cols) < 2) {
            stop(paste(
                "This diagnostic needs a per-cell score for each candidate cell",
                "type, as two or more numeric columns in colData. The preset",
                "datasets carry a single confidence score rather than a full",
                "score matrix, so there is nothing to compute entropy over.",
                "Upload an object whose colData holds one score column per cell",
                "type to use it."
            ))
        }
        cd <- SummarizedExperiment::colData(ctx$focus)
        X <- t(as.matrix(as.data.frame(cd[, cols, drop = FALSE])))
        list(X = X, plot = TRUE, verbose = FALSE,
             inverseNormalTransformationform = args$inverseNormalTransformationform %||% FALSE)
    },
    repro_prelude = function(ctx, args, param_values, ref_var, query_var) {
        cols <- param_values$score_columns
        c(
            "# Scores are assembled into a cell-type-by-cell matrix.",
            sprintf("X <- t(as.matrix(as.data.frame(colData(sce)[, %s])))",
                    paste(deparse(cols), collapse = ""))
        )
    },
    repro_names = list(X = "X"),
    available = needs_numeric_cols(2, "two or more per-cell-type score columns"),
    needs = "focus",
    returns = "value"
))

# ===========================================================================
# Utilities
# ===========================================================================

add_entry(diag_entry(
    id = "projectPCA",
    fn = "projectPCA",
    label = "Project onto reference PCA",
    category = "utility",
    question = "What are the raw projected coordinates?",
    blurb = "The reference PCA scores and the query projected into the same space, as a table.",
    reading = paste(
        "This is the table every projection-based diagnostic is built on.",
        "Download it to reproduce a figure outside the app or to join the",
        "coordinates onto your own metadata."
    ),
    data_map = map_pair(),
    params = list(
        param_pcs(value = 1:10, max = 30), param_celltypes(),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    tables = function(x) list("Projected coordinates" = x),
    returns = "value"
))

add_entry(diag_entry(
    id = "projectSIR",
    fn = "projectSIR",
    label = "Project onto SIR space",
    category = "utility",
    question = "What are the raw sliced inverse regression coordinates?",
    blurb = "Reference and query cells projected into the sliced inverse regression space.",
    reading = "The coordinate table underlying the sliced inverse regression view.",
    data_map = map_pair(pcs = FALSE),
    params = list(
        param_celltypes(),
        param_bool("multiple_cond_means", "Multiple conditional means", TRUE, advanced = TRUE),
        param_num("cumulative_variance_threshold", "Cumulative variance threshold",
                  0.7, 0.1, 0.99, 0.05, advanced = TRUE),
        param_int("n_neighbor", "Neighbours", 1, 1, 50, advanced = TRUE),
        param_int("max_cells_query", "Max query cells", 5000, 100, 50000, advanced = TRUE),
        param_int("max_cells_ref", "Max reference cells", 5000, 100, 50000, advanced = TRUE)
    ),
    returns = "value",
    heavy = TRUE
))

add_entry(diag_entry(
    id = "processPCA",
    fn = "processPCA",
    label = "Compute PCA",
    category = "utility",
    question = "Does this object have a usable PCA, and what happens if it does not?",
    blurb = paste(
        "Checks the object's PCA and computes one from highly variable genes if",
        "it is missing or invalid."
    ),
    reading = paste(
        "The app runs this automatically on upload, so this panel is mainly for",
        "inspecting what it did: how many components exist, how much variance",
        "they carry, and whether the stored PCA was reused or recomputed."
    ),
    data_map = list(sce_object = "focus", assay_name = "assay"),
    params = list(
        param_int("n_hvgs", "Highly variable genes", 2000, 100, 10000),
        param_int("max_cells", "Max cells if recomputing", 5000, 200, 100000, advanced = TRUE)
    ),
    summarise = function(x) {
        pca <- SingleCellExperiment::reducedDim(x, "PCA")
        pv <- attributes(pca)[["percentVar"]]
        list(
            figure("Components", ncol(pca)),
            figure("Cells", n_fmt(nrow(pca))),
            figure("Variance in first 5", sprintf("%.1f%%", sum(utils::head(pv, 5))))
        )
    },
    tables = function(x) {
        pca <- SingleCellExperiment::reducedDim(x, "PCA")
        pv <- attributes(pca)[["percentVar"]]
        list("Variance explained" = data.frame(
            component = paste0("PC", seq_along(pv)),
            percent_variance = pv,
            cumulative = cumsum(pv)
        ))
    },
    needs = "focus",
    returns = "sce"
))

#' Logical coercion that tolerates NA, used by summaries
isTRUE_vec <- function(x) !is.na(x) & x

#' Subset an SCE to one or more cell types
subset_to_types <- function(sce, col, types) {
    keep <- which(as.character(SummarizedExperiment::colData(sce)[[col]]) %in% types)
    if (!length(keep)) {
        stop(sprintf("No cells of type '%s' in this dataset.", paste(types, collapse = ", ")))
    }
    subset_cells(sce, keep)
}
