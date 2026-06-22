# netset.R ------------------------------------------------------------------
#
# The core data model: anansi_network (one network) and anansi_netset (a set).
# See docs/DESIGN.md (section 5) for the data model rationale.

#' The `anansi_network` class
#'
#' A single phylogenetic network, created by [parse_network()] (or
#' [read_networks_csv()] for a whole set). It is a list with elements:
#' \describe{
#'   \item{`evonet`}{the `ape::evonet` (backbone tree in `$edge`, hybrid edges
#'     in `$reticulation`).}
#'   \item{`enewick`}{the raw extended-Newick string.}
#'   \item{`gammas`}{a data.frame of inheritance probabilities (see
#'     [parse_gammas()]).}
#'   \item{`hybrid_nodes`}{named integer vector mapping `#Hk` labels to node ids.}
#'   \item{`meta`}{per-network metadata (rank, log-likelihood, ...).}
#'   \item{`issues`}{character vector of parsing issues, if any.}
#' }
#' @name anansi_network
#' @seealso [anansi_netset], [parse_network]
NULL

#' Construct a set of networks over a fixed taxon set
#'
#' Networks whose taxon set diverges from the modal set are kept but flagged in
#' `taxa_ok` (filter with `ns[ns$taxa_ok]`); with `validate = TRUE` a warning
#' summarizes how many diverge.
#'
#' @param networks A list of [anansi_network] objects (or a single one).
#' @param meta Optional data.frame of per-network metadata (one row per network).
#' @param validate If TRUE (default), warn when the taxon set is not fixed.
#' @return An `anansi_netset` with fields `networks`, `taxa` (the modal set),
#'   `taxa_ok` (logical per network), and `meta`.
#' @seealso [read_networks_csv], [taxa_consistency]
#' @export
anansi_netset <- function(networks, meta = NULL, validate = TRUE) {
  if (inherits(networks, "anansi_network")) networks <- list(networks)
  if (length(networks) == 0L) rlang::abort("No networks provided.")
  if (!all(vapply(networks, inherits, logical(1), "anansi_network"))) {
    rlang::abort("All elements of `networks` must be anansi_network objects.")
  }
  tc <- taxa_consistency(networks)
  out <- structure(list(networks = networks, taxa = tc$reference,
                        taxa_ok = tc$ok, meta = meta),
                   class = "anansi_netset")
  if (validate && tc$n_bad > 0L) {
    warning(sprintf(
      paste0("%d of %d networks do not match the modal taxon set (%d taxa) and ",
             "are flagged taxa_ok = FALSE (e.g. index %s). These are likely ",
             "read.evonet artifacts on nested/stacked reticulations; see ",
             "docs/DATA.md. Keep the consistent ones with ns[ns$taxa_ok]."),
      tc$n_bad, length(networks), length(tc$reference),
      paste(utils::head(tc$bad, 10L), collapse = ", ")), call. = FALSE)
  }
  out
}

#' @export
length.anansi_netset <- function(x) length(x$networks)

#' Subset a netset, preserving class and aligned metadata/flags
#' @param x An [anansi_netset].
#' @param i Indices (or logical) of networks to keep.
#' @export
`[.anansi_netset` <- function(x, i) {
  meta <- if (!is.null(x$meta)) x$meta[i, , drop = FALSE] else NULL
  structure(list(networks = x$networks[i], taxa = x$taxa,
                 taxa_ok = x$taxa_ok[i], meta = meta),
            class = "anansi_netset")
}

#' @export
print.anansi_network <- function(x, ...) {
  cat("<anansi_network>\n")
  cat("  taxa:         ", length(x$evonet$tip.label), "\n", sep = "")
  cat("  reticulations:", n_reticulations(x), "\n", sep = "")
  show <- intersect(c("dataset", "rank", "network_rank", "reticulations",
                      "num_reticulations", "log_probability",
                      "log_pseudo_likelihood"), names(x$meta))
  if (length(show)) {
    vals <- vapply(show, function(k) {
      v <- x$meta[[k]]
      if (is.numeric(v)) format(v, digits = 6) else as.character(v)
    }, character(1))
    cat("  meta:         ", paste(sprintf("%s=%s", show, vals), collapse = "  "),
        "\n", sep = "")
  }
  if (length(x$issues)) {
    cat("  issues:       ", paste(x$issues, collapse = "; "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.anansi_netset <- function(x, ...) {
  nr <- vapply(x$networks, n_reticulations, integer(1))
  tab <- table(nr)
  n_issues <- sum(vapply(x$networks, function(n) length(n$issues) > 0L, logical(1)))
  cat("<anansi_netset>\n")
  cat("  networks:", length(x$networks), "\n")
  cat("  taxa:    ", length(x$taxa), "(modal/reference set)\n")
  cat("  reticulations/network:",
      paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse = "  "),
      "\n")
  cat("  taxa_ok: ", sum(x$taxa_ok), "/", length(x$taxa_ok),
      " networks match the reference taxon set\n", sep = "")
  if (n_issues > 0L) {
    cat("  flagged: ", n_issues, " network(s) carry parsing issues",
        " (see $networks[[i]]$issues)\n", sep = "")
  }
  if (!is.null(x$meta)) {
    cat("  meta columns:", paste(names(x$meta), collapse = ", "), "\n")
  }
  invisible(x)
}
