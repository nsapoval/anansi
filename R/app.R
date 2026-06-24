# app.R ---------------------------------------------------------------------
#
# The anansi Shiny web app: an in-package launcher plus the thin, testable
# helpers the bundled app (inst/shiny/app.R) relies on. shiny and the UI
# polish packages (bslib, colourpicker, svglite) are Suggests, so everything
# here is guarded and the helpers themselves carry no hard dependency on them.

# Treat an empty / all-blank selector as "unset" (NULL) so it falls through to
# the densinet() default. Internal.
.nz <- function(v) {
  if (is.null(v) || !length(v)) return(NULL)
  v <- v[!is.na(v)]
  v <- v[nzchar(as.character(v))]
  if (!length(v)) NULL else v
}

# Treat a blank / NA numeric input as "unset" (NULL). Internal. Used for the
# auto-scaled densinet() opacities (alpha, ret_alpha) whose default is NULL.
.nz_num <- function(v) {
  if (is.null(v) || length(v) != 1L || is.na(v)) NULL else v
}

#' Build a densinet() figure from a flat parameter list (web-app helper)
#'
#' Maps a flat list of UI parameters (e.g. Shiny inputs) onto a [densinet] call,
#' so the app server stays a one-liner and the input -> argument mapping is
#' unit-testable on its own. Empty selectors and blank text are treated as unset
#' and fall through to the [densinet] defaults.
#'
#' Any `outgroup` taxa are force-kept in the `keep` subset: an outgroup that is
#' pruned away would make the internal tip-ordering fail to pin it, so
#' `keep <- union(keep, outgroup)`.
#'
#' @param netset An [anansi_netset] (see [read_networks_csv] /
#'   [netset_from_enewick]).
#' @param params Named list of [densinet] arguments. Recognized names mirror the
#'   [densinet] formals; unknown names are ignored.
#' @return A `ggplot` object (as returned by [densinet]).
#' @seealso [run_anansi_app], [densinet]
#' @examples
#' ns <- read_networks_csv(anansi_example())
#' p <- build_densinet(ns, list(method = "first", reticulation_style = "hybrid"))
#' @export
build_densinet <- function(netset, params = list()) {
  if (!inherits(netset, "anansi_netset")) {
    rlang::abort("build_densinet() needs an anansi_netset (see read_networks_csv()).")
  }
  p <- params
  outgroup <- .nz(p$outgroup)
  keep     <- .nz(p$keep)
  # An outgroup taxon must survive pruning, or outgroup pinning cannot place it.
  if (!is.null(keep) && !is.null(outgroup)) keep <- union(keep, outgroup)
  top_n <- if (is.null(p$top_n) || is.na(p$top_n) || p$top_n <= 0) NULL else as.integer(p$top_n)

  densinet(
    netset,
    method             = p$method %||% "mode",
    reticulation_style = p$reticulation_style %||% "arrow",
    mode               = p$mode %||% "cladogram",
    layout             = p$layout %||% "slanted",
    consensus_p        = p$consensus_p %||% 0.5,
    top_n              = top_n,
    top_by             = .nz(p$top_by),
    outgroup           = outgroup,
    outgroup_position  = p$outgroup_position %||% "top",
    keep               = keep,
    title              = .nz(p$title),
    tip_labels         = p$tip_labels %||% TRUE,
    # --- advanced (NULL/unset falls through to densinet() defaults) ---
    snap_to_consensus  = p$snap_to_consensus %||% FALSE,
    color_by_support   = p$color_by_support %||% TRUE,
    consensus_ret      = p$consensus_ret %||% TRUE,
    consistent_only    = p$consistent_only %||% TRUE,
    consensus_ret_min  = p$consensus_ret_min %||% 0.1,
    ret_edge_frac      = p$ret_edge_frac %||% 0.1,
    linewidth          = p$linewidth %||% 0.3,
    consensus_linewidth = p$consensus_linewidth %||% 0.7,
    jitter             = p$jitter %||% 0,
    tip_size           = p$tip_size %||% 3,
    tip_offset         = p$tip_offset %||% 0.02,
    alpha              = .nz_num(p$alpha),
    ret_alpha          = .nz_num(p$ret_alpha),
    ret_linetype       = p$ret_linetype %||% "dashed",
    tree_color         = p$tree_color %||% "steelblue",
    ret_color          = p$ret_color %||% "firebrick",
    consensus_color    = p$consensus_color %||% "black",
    consensus_ret_color = p$consensus_ret_color %||% "darkred",
    # --- gradient cloud + explicit backbone (unset falls through to defaults) ---
    color_by           = .nz(p$color_by),
    color_palette      = p$color_palette %||% "viridis",
    color_direction    = p$color_direction %||% 1,
    color_low          = p$color_low %||% "grey80",
    color_high         = p$color_high %||% "firebrick",
    color_legend       = .nz(p$color_legend),
    color_trim         = .nz_num(p$color_trim),
    color_trim_action  = p$color_trim_action %||% "clamp",
    backbone           = .nz(p$backbone),
    backbone_color     = p$backbone_color %||% "black",
    backbone_ret_color = p$backbone_ret_color %||% "darkred")
}

#' Build a network set from extended-Newick strings
#'
#' Parses one or more extended-Newick strings (one per line) into an
#' [anansi_netset], the text-input complement to [read_networks_csv]. Blank
#' lines and lines beginning with `#` (comments) are skipped. Each network is
#' parsed with [parse_network] and gets an `index` metadata column.
#'
#' @param text A character vector of extended-Newick strings, or a single string
#'   with one network per line.
#' @param validate If TRUE (default), warn when the taxon set is not fixed (see
#'   [anansi_netset]).
#' @return An [anansi_netset].
#' @seealso [read_networks_csv], [parse_network], [run_anansi_app]
#' @examples
#' txt <- c("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);",
#'          "((A:1,(C:1)#H1:1::0.6):1,(B:1,#H1:1::0.4):1);")
#' netset_from_enewick(txt)
#' @export
netset_from_enewick <- function(text, validate = TRUE) {
  if (!is.character(text)) rlang::abort("`text` must be character.")
  if (length(text) == 1L) text <- strsplit(text, "\r?\n")[[1]]
  lines <- trimws(text)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) rlang::abort("No extended-Newick strings found in input.")
  nets <- lapply(seq_along(lines), function(i) {
    parse_network(lines[i], meta = list(index = i))
  })
  anansi_netset(nets, meta = data.frame(index = seq_along(lines)),
                validate = validate)
}

#' Launch the anansi Shiny web app
#'
#' Opens an interactive web interface for building a [densinet] figure: upload a
#' CSV of extended-Newick networks (see [read_networks_csv]), paste/upload raw
#' extended-Newick strings (see [netset_from_enewick]), or load the bundled
#' example ([anansi_example]); tune the key and advanced controls; preview the
#' figure; and download it as PNG/PDF (or SVG when \pkg{svglite} is installed).
#'
#' Requires the suggested package \pkg{shiny}. A richer UI also uses
#' \pkg{bslib} (layout) and \pkg{colourpicker} (colour inputs); when these are
#' absent the app degrades to a plain layout and text colour fields.
#'
#' @param ... Passed to [shiny::runApp] (e.g. `port`, `launch.browser`, `host`).
#' @return Invisibly `NULL`; called for the side effect of running the app.
#' @seealso [build_densinet], [densinet]
#' @examples
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   run_anansi_app()
#' }
#' @export
run_anansi_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    rlang::abort(paste0(
      "The anansi web app needs the 'shiny' package. Install it with:\n",
      "  install.packages(\"shiny\")\n",
      "For the full UI also install: bslib, colourpicker, svglite."))
  }
  app_dir <- system.file("shiny", package = "anansi")
  if (!nzchar(app_dir)) {
    rlang::abort("Could not find the bundled Shiny app (inst/shiny is missing).")
  }
  shiny::runApp(app_dir, ...)
  invisible(NULL)
}
