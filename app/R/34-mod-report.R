# Report export.
#
# Assembles the chosen diagnostics into a single self-contained HTML file and
# a standalone R script. The HTML is built directly with htmltools and
# base64-embedded images rather than through pandoc, so it works on any
# deployment host without an external toolchain.

REPORT_DEFAULTS <- c("plotCellTypePCA", "detectAnomaly", "calculateGeneShifts",
                     "regressPC", "calculateWassersteinDistance")

report_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "sc-page",
        shiny::div(
            class = "sc-page-head",
            shiny::h2("Report", class = "sc-page-title"),
            shiny::p(class = "sc-page-lede",
                     paste("Pick the diagnostics to include. You get one self-contained HTML",
                           "file and an R script that regenerates every figure outside the app."))
        ),
        shiny::uiOutput(ns("gate")),
        bslib::card(
            bslib::card_header("Contents"),
            bslib::card_body(
                shiny::uiOutput(ns("picker")),
                shiny::div(
                    class = "row g-3 mt-2",
                    shiny::div(class = "col-md-8",
                               shiny::textInput(ns("title"), "Report title",
                                                value = "Annotation diagnostics",
                                                width = "100%")),
                    shiny::div(class = "col-md-4",
                               shiny::textInput(ns("author"), "Author", value = "", width = "100%"))
                ),
                shiny::div(
                    class = "sc-config-actions",
                    shiny::downloadButton(ns("dl_html"), "Download HTML report",
                                          class = "btn-primary"),
                    shiny::downloadButton(ns("dl_script"), "Download R script",
                                          class = "btn-outline-primary")
                ),
                shiny::p(class = "small text-muted mt-2",
                         "Building the report runs each selected diagnostic, reusing anything already computed.")
            )
        )
    )
}

report_server <- function(id, ctx_r, session_cache, active = shiny::reactive(TRUE)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$gate <- shiny::renderUI({
            ctx <- ctx_r()
            if (isTRUE(ctx$ready)) return(NULL)
            alert("Choose data on the Data page before building a report.", type = "info")
        })

        always_render(output, c("gate", "picker"))

        output$picker <- shiny::renderUI({
            ctx <- ctx_r()
            groups <- lapply(names(CATEGORIES), function(ck) {
                entries <- registry_by_category(ck)
                entries <- Filter(function(e) is.null(entry_ready(e, ctx)), entries)
                if (!length(entries)) return(NULL)
                choices <- stats::setNames(vapply(entries, `[[`, character(1), "id"),
                                           vapply(entries, `[[`, character(1), "label"))
                shiny::div(
                    class = "sc-report-group",
                    shiny::div(class = "sc-report-group-title", CATEGORIES[[ck]]$label),
                    shiny::checkboxGroupInput(
                        ns(paste0("pick_", ck)), NULL, choices = choices,
                        selected = intersect(choices, REPORT_DEFAULTS)
                    )
                )
            })
            bslib::layout_column_wrap(width = "260px", fixed_width = FALSE,
                                      !!!not_null(groups))
        })

        selected <- shiny::reactive({
            ids <- unlist(lapply(names(CATEGORIES), function(ck) input[[paste0("pick_", ck)]]))
            ids <- ids[!is.na(ids)]
            # Keep registry order so the report reads in a sensible sequence.
            intersect(names(REGISTRY), ids)
        })

        run_one <- function(entry_id, ctx) {
            entry <- registry_get(entry_id)
            pv <- default_param_values(entry$params, ctx)
            args <- tryCatch(build_args(entry, ctx, pv),
                             error = function(e) sc_failure(conditionMessage(e)))
            if (is_failure(args)) return(list(entry = entry, value = args, args = args, params = pv))
            tier <- cache_for(ctx, session_cache)
            key <- cache_key(entry$id, ctx, pv)
            value <- with_cache(key, tier$cache, tier$shareable, run_entry(entry, args))
            list(entry = entry, value = value, args = args, params = pv)
        }

        gather <- function() {
            ctx <- ctx_r()
            ids <- selected()
            shiny::validate(shiny::need(length(ids) > 0, "Select at least one diagnostic."))
            prog <- shiny::Progress$new(session, min = 0, max = length(ids))
            on.exit(prog$close(), add = TRUE)
            out <- list()
            for (i in seq_along(ids)) {
                prog$set(message = sprintf("Running %s", registry_get(ids[i])$label),
                         value = i - 1)
                out[[ids[i]]] <- run_one(ids[i], ctx)
            }
            list(ctx = ctx, results = out)
        }

        output$dl_html <- shiny::downloadHandler(
            filename = function() sprintf("scdiagnostics-report-%s.html", Sys.Date()),
            contentType = "text/html",
            content = function(file) {
                g <- gather()
                html <- build_report_html(g$ctx, g$results,
                                          title = input$title %|% "Annotation diagnostics",
                                          author = input$author)
                writeLines(html, file, useBytes = TRUE)
            }
        )

        output$dl_script <- shiny::downloadHandler(
            filename = function() sprintf("scdiagnostics-%s.R", Sys.Date()),
            content = function(file) {
                ctx <- ctx_r()
                ids <- selected()
                shiny::validate(shiny::need(length(ids) > 0, "Select at least one diagnostic."))
                writeLines(build_report_script(ctx, ids), file)
            }
        )
    })
}

# --- rendering ---------------------------------------------------------------

#' Render a plot to a base64 data URI
plot_data_uri <- function(expr, width = 1800, height = 1200, res = 170) {
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path), add = TRUE)
    grDevices::png(path, width = width, height = height, res = res)
    ok <- tryCatch({
        force(expr)
        TRUE
    }, error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("Plot failed:", conditionMessage(e)), cex = 0.9)
        TRUE
    })
    grDevices::dev.off()
    if (!ok || !file.exists(path)) return(NULL)
    raw <- readBin(path, "raw", file.info(path)$size)
    paste0("data:image/png;base64,", base64_encode(raw))
}

#' base64 without adding a dependency
#'
#' Prefers base64enc, and falls back to the encoder inside httpuv, which is a
#' hard dependency of shiny and therefore always present.
base64_encode <- function(raw) {
    if (requireNamespace("base64enc", quietly = TRUE)) {
        return(base64enc::base64encode(raw))
    }
    enc <- tryCatch(utils::getFromNamespace("base64encode", "httpuv"),
                    error = function(e) NULL)
    if (!is.null(enc)) return(enc(raw))
    stop("No base64 encoder available.", call. = FALSE)
}

report_css <- "
:root { --ink:#16242e; --muted:#5b6b78; --line:#dde3e7; --brand:#1a6985; }
* { box-sizing:border-box; }
body { font-family: system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;
       color:var(--ink); max-width:920px; margin:0 auto; padding:2.5rem 1.25rem 4rem;
       line-height:1.6; }
h1 { font-size:1.9rem; margin:0 0 .25rem; }
h2 { font-size:1.3rem; margin:2.5rem 0 .25rem; padding-top:1.5rem;
     border-top:1px solid var(--line); }
h3 { font-size:1rem; margin:1.5rem 0 .5rem; color:var(--muted);
     text-transform:uppercase; letter-spacing:.04em; }
.meta { color:var(--muted); font-size:.9rem; margin-bottom:2rem; }
.lede { color:var(--muted); margin:.25rem 0 1rem; }
.figs { display:flex; flex-wrap:wrap; gap:.75rem; margin:1rem 0; }
.fig { border:1px solid var(--line); border-radius:8px; padding:.7rem .9rem; min-width:170px; }
.fig .l { font-size:.72rem; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); }
.fig .v { font-size:1.4rem; font-weight:600; }
.fig .c { font-size:.78rem; color:var(--muted); }
.fig.warning { border-color:#e6c07a; background:#fdf8ee; }
.fig.danger  { border-color:#e0a3ab; background:#fdf1f2; }
.fig.success { border-color:#a3cbb6; background:#f1f9f4; }
img { max-width:100%; border:1px solid var(--line); border-radius:8px; }
table { border-collapse:collapse; width:100%; font-size:.85rem; margin:.75rem 0; }
th,td { border-bottom:1px solid var(--line); padding:.4rem .5rem; text-align:left; }
th { background:#f6f8f9; font-weight:600; }
td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; }
pre { background:#f6f8f9; border:1px solid var(--line); border-radius:8px;
      padding:.8rem; overflow-x:auto; font-size:.8rem; }
.err { background:#fdf1f2; border:1px solid #e0a3ab; border-radius:8px; padding:.8rem; }
.foot { margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line);
        font-size:.82rem; color:var(--muted); }
"

#' Small HTML table from a data frame
report_table <- function(df, max_rows = 25) {
    df <- utils::head(df, max_rows)
    num <- vapply(df, is.numeric, logical(1))
    fmt <- function(v, isnum) {
        if (!isnum) return(as.character(v))
        ifelse(is.na(v), "", ifelse(v == round(v) & abs(v) < 1e9,
                                    formatC(v, format = "d", big.mark = ","),
                                    formatC(signif(v, 4), format = "g")))
    }
    head_cells <- mapply(function(nm, isnum) {
        htmltools::tags$th(class = if (isnum) "num" else NULL, nm)
    }, names(df), num, SIMPLIFY = FALSE)
    rows <- lapply(seq_len(nrow(df)), function(i) {
        htmltools::tags$tr(mapply(function(j, isnum) {
            htmltools::tags$td(class = if (isnum) "num" else NULL, fmt(df[i, j], isnum))
        }, seq_along(df), num, SIMPLIFY = FALSE))
    })
    htmltools::tags$table(
        htmltools::tags$thead(htmltools::tags$tr(head_cells)),
        htmltools::tags$tbody(rows)
    )
}

#' Build the complete report document
build_report_html <- function(ctx, results, title, author = "") {
    sections <- lapply(names(results), function(id) {
        r <- results[[id]]
        e <- r$entry
        body <- if (is_failure(r$value)) {
            htmltools::div(class = "err",
                           htmltools::tags$b("Could not be computed. "), r$value$message)
        } else {
            figs <- if (!is.null(e$summarise)) {
                tryCatch(e$summarise(r$value), error = function(e) list())
            } else list()
            tabs <- if (!is.null(e$tables)) {
                tryCatch(e$tables(r$value), error = function(e) NULL)
            } else NULL
            if (is.null(tabs) || !length(tabs)) tabs <- auto_tables(r$value, max_tables = 2)
            tabs <- Filter(function(d) !is.null(d) && NROW(d) > 0, tabs)

            uri <- plot_data_uri(draw_result(build_plot(e, r$value, list(), ctx)))

            htmltools::tagList(
                if (length(figs)) {
                    htmltools::div(class = "figs", lapply(figs, function(f) {
                        htmltools::div(
                            class = paste("fig", f$tone %||% ""),
                            htmltools::div(class = "l", f$label),
                            htmltools::div(class = "v", f$value),
                            if (!is.null(f$caption)) htmltools::div(class = "c", f$caption)
                        )
                    }))
                },
                if (!is.null(uri)) htmltools::tags$img(src = uri, alt = e$label),
                if (length(tabs)) {
                    htmltools::tagList(
                        htmltools::tags$h3(names(tabs)[1]),
                        report_table(tabs[[1]])
                    )
                },
                htmltools::tags$h3("How to read this"),
                htmltools::tags$p(e$reading %||% e$blurb),
                htmltools::tags$h3("R code"),
                htmltools::tags$pre(repro_code(e, ctx, r$args, NULL, r$params))
            )
        }
        htmltools::tagList(
            htmltools::tags$h2(e$label),
            htmltools::tags$p(class = "lede", e$question %||% e$blurb),
            body
        )
    })

    doc <- htmltools::tags$html(
        lang = "en",
        htmltools::tags$head(
            htmltools::tags$meta(charset = "utf-8"),
            htmltools::tags$meta(name = "viewport",
                                 content = "width=device-width, initial-scale=1"),
            htmltools::tags$title(title),
            htmltools::tags$style(htmltools::HTML(report_css))
        ),
        htmltools::tags$body(
            htmltools::tags$h1(title),
            htmltools::div(
                class = "meta",
                paste0(
                    if (nzchar(author %|% "")) paste0(author, " · ") else "",
                    format(Sys.time(), "%d %B %Y, %H:%M"), " · ",
                    "scDiagnostics ", utils::packageVersion("scDiagnostics"), " · ",
                    R.version.string
                )
            ),
            htmltools::tags$h2("Data"),
            report_table(data.frame(
                field = c("Reference", "Reference annotation", "Query",
                          "Query annotation", "Assay", "Components", "Cell types"),
                value = c(
                    sprintf("%s (%s genes x %s cells)", ctx$ref_desc$label %||% "-",
                            n_fmt(nrow(ctx$ref)), n_fmt(ncol(ctx$ref))),
                    ctx$ref_col %||% "-",
                    if (is.null(ctx$query)) "-" else
                        sprintf("%s (%s genes x %s cells)", ctx$query_desc$label %||% "-",
                                n_fmt(nrow(ctx$query)), n_fmt(ncol(ctx$query))),
                    ctx$query_col %||% "-",
                    ctx$assay,
                    sprintf("PC%d-PC%d", min(ctx$pc_subset), max(ctx$pc_subset)),
                    if (length(ctx$cell_types)) paste(ctx$cell_types, collapse = ", ")
                    else paste(ctx$shared_types, collapse = ", ")
                ),
                stringsAsFactors = FALSE
            ), max_rows = 20),
            sections,
            htmltools::div(
                class = "foot",
                htmltools::HTML(paste(
                    "Generated by the scDiagnostics app. If you use these results, please cite:",
                    "Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026).",
                    "<em>scDiagnostics: systematic assessment of cell type annotation in",
                    "single-cell transcriptomics data.</em> Briefings in Bioinformatics, 27(5),",
                    "bbag496. <a href='https://doi.org/10.1093/bib/bbag496'>doi:10.1093/bib/bbag496</a>."
                ))
            )
        )
    )
    as.character(htmltools::doRenderTags(doc))
}

#' One script reproducing every selected diagnostic
build_report_script <- function(ctx, ids) {
    header <- c(
        "# Reproduces the diagnostics selected in the scDiagnostics app.",
        sprintf("# Generated %s with scDiagnostics %s.",
                format(Sys.Date()), utils::packageVersion("scDiagnostics")),
        "",
        "library(scDiagnostics)",
        "library(SingleCellExperiment)",
        ""
    )
    blocks <- lapply(ids, function(id) {
        e <- registry_get(id)
        pv <- default_param_values(e$params, ctx)
        args <- tryCatch(build_args(e, ctx, pv), error = function(err) NULL)
        if (is.null(args)) {
            return(c(sprintf("# --- %s: skipped, could not resolve arguments", e$label), ""))
        }
        code <- repro_code(e, ctx, args, NULL, pv)
        # The per-diagnostic code repeats the library and data lines; drop them
        # here because the script has its own header.
        lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
        lines <- lines[!grepl("^library\\(", lines)]
        c(sprintf("# --- %s -------------------------------------------------", e$label),
          sprintf("# %s", e$question %||% e$blurb),
          sub("^result <- ", sprintf("%s <- ", id), lines),
          "")
    })
    paste(c(header, unlist(blocks)), collapse = "\n")
}

# --- about -------------------------------------------------------------------

about_ui <- function() {
    shiny::div(
        class = "sc-page sc-about",
        shiny::div(
            class = "sc-page-head",
            shiny::h2("About", class = "sc-page-title")
        ),
        bslib::layout_column_wrap(
            width = "420px", fixed_width = FALSE,
            bslib::card(
                bslib::card_header("The package"),
                bslib::card_body(shiny::tagList(
                    shiny::p(paste(
                        "This app is a front end for the scDiagnostics Bioconductor package.",
                        "Every statistic and every figure comes from the package itself:",
                        "the app holds no copy of any method, so results here match results",
                        "in an R session."
                    )),
                    shiny::p(
                        shiny::a("Package documentation",
                                 href = "https://ccb-hms.github.io/scDiagnostics/",
                                 target = "_blank", rel = "noopener"), " · ",
                        shiny::a("Source",
                                 href = "https://github.com/ccb-hms/scDiagnostics",
                                 target = "_blank", rel = "noopener"), " · ",
                        shiny::a("Bioconductor",
                                 href = "https://bioconductor.org/packages/scDiagnostics",
                                 target = "_blank", rel = "noopener")
                    )
                ))
            ),
            bslib::card(
                bslib::card_header("Citation"),
                bslib::card_body(shiny::tagList(
                    shiny::p(
                        "Christidis A, Ghazi A, Chawla S, Turaga N, Gentleman R, Geistlinger L (2026).",
                        shiny::tags$em(paste("scDiagnostics: systematic assessment of cell type",
                                             "annotation in single-cell transcriptomics data.")),
                        "Briefings in Bioinformatics, 27(5), bbag496."
                    ),
                    shiny::p(shiny::a("doi:10.1093/bib/bbag496",
                                      href = "https://doi.org/10.1093/bib/bbag496",
                                      target = "_blank", rel = "noopener"))
                ))
            ),
            bslib::card(
                bslib::card_header("Privacy"),
                bslib::card_body(shiny::p(paste(
                    "Uploaded data is held in the memory of your own session and in a",
                    "temporary file on the server, both discarded when the session ends.",
                    "Results computed from uploads are cached only for your session.",
                    "Only results from the built-in datasets are cached across users."
                )))
            ),
            bslib::card(
                bslib::card_header("Session"),
                bslib::card_body(shiny::tags$pre(
                    class = "sc-snippet",
                    paste(
                        R.version.string,
                        paste("scDiagnostics", utils::packageVersion("scDiagnostics")),
                        paste("shiny", utils::packageVersion("shiny")),
                        paste("bslib", utils::packageVersion("bslib")),
                        sep = "\n"
                    )
                ))
            )
        )
    )
}
