# Visual theme.
#
# Bootstrap 5 via bslib, with a system font stack rather than a web font so the
# app has no network dependency at render time on the deployment host.

SC_COLORS <- list(
    primary   = "#1a6985",
    secondary = "#5b6b78",
    success   = "#2e7d5b",
    info      = "#3a7ca5",
    warning   = "#b5730c",
    danger    = "#a8323f",
    ink       = "#16242e",
    surface   = "#f6f8f9"
)

SANS <- paste(
    "system-ui", "-apple-system", "'Segoe UI'", "Roboto", "'Helvetica Neue'",
    "Arial", "sans-serif",
    sep = ", "
)
MONO <- paste("'SFMono-Regular'", "Menlo", "Consolas", "'Liberation Mono'", "monospace", sep = ", ")

sc_theme <- function() {
    bslib::bs_theme(
        version = 5,
        primary = SC_COLORS$primary,
        secondary = SC_COLORS$secondary,
        success = SC_COLORS$success,
        info = SC_COLORS$info,
        warning = SC_COLORS$warning,
        danger = SC_COLORS$danger,
        "body-color" = SC_COLORS$ink,
        "font-size-base" = "0.95rem",
        "headings-font-weight" = "600",
        "card-border-color" = "rgba(22, 36, 46, 0.12)",
        "card-cap-bg" = "rgba(26, 105, 133, 0.05)",
        "navbar-padding-y" = "0.5rem",
        base_font = SANS,
        heading_font = SANS,
        code_font = MONO
    )
}

#' A ggplot theme matched to the app chrome
#'
#' Applied as the global default so package plots inherit it without the app
#' having to modify any of them.
sc_ggtheme <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            text = ggplot2::element_text(colour = SC_COLORS$ink),
            plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.1),
            plot.subtitle = ggplot2::element_text(colour = SC_COLORS$secondary),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(colour = "#e3e8eb", linewidth = 0.35),
            strip.text = ggplot2::element_text(face = "bold", size = base_size * 0.9),
            axis.title = ggplot2::element_text(size = base_size * 0.9),
            legend.position = "bottom",
            legend.title = ggplot2::element_text(size = base_size * 0.85),
            plot.margin = ggplot2::margin(8, 10, 8, 8)
        )
}

#' Colour tone to a Bootstrap class suffix
tone_class <- function(tone) {
    switch(tone %||% "default",
        success = "sc-tone-success",
        warning = "sc-tone-warning",
        danger = "sc-tone-danger",
        "sc-tone-default")
}
