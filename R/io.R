# io.R ----------------------------------------------------------------------
#
# Reading sets of phylogenetic networks and parsing extended Newick (eNewick).
# See docs/DATA.md for the formats / known parsing issues and docs/DESIGN.md
# (D8) for the reuse strategy.

#' Parse gamma inheritance probabilities from an extended-Newick string
#'
#' `ape::read.evonet()` parses the network *structure* but does not populate the
#' gamma (inheritance) probabilities, so we recover them from the raw string.
#' Each hybrid `#Hk` appears twice: the first occurrence is on the major (tree)
#' edge, the second on the minor (reticulation) edge. This generalizes the
#' reference `plot_network.R` logic to an arbitrary number of hybrids.
#'
#' @param enewick A single extended-Newick string.
#' @return A data.frame with columns `label`, `tree_gamma` (major edge) and
#'   `ret_gamma` (minor edge); zero rows if there are no hybrids.
#' @examples
#' parse_gammas("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);")
#' @export
parse_gammas <- function(enewick) {
  stopifnot(is.character(enewick), length(enewick) == 1L)
  empty <- data.frame(label = character(), tree_gamma = numeric(),
                      ret_gamma = numeric(), stringsAsFactors = FALSE)
  m <- regmatches(enewick,
                  gregexpr("#H\\d+:[^,()]+::[0-9.eE+-]+", enewick, perl = TRUE))[[1]]
  if (length(m) == 0L) return(empty)
  labels <- sub("^(#H\\d+):.*", "\\1", m)
  gammas <- as.numeric(sub("^#H\\d+:[^:]+::", "", m))
  do.call(rbind, lapply(unique(labels), function(lbl) {
    g <- gammas[labels == lbl]
    data.frame(label = lbl,
               tree_gamma = if (length(g) >= 1L) g[1] else NA_real_,
               ret_gamma  = if (length(g) >= 2L) g[2] else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

# Read one eNewick string to an evonet, robustly.
# - Networks with no hybrids make read.evonet() fail/degenerate; fall back to
#   read.tree() and attach an empty reticulation matrix.
# - Always guarantee a `$reticulation` matrix. Internal.
parse_one_evonet <- function(enewick) {
  net <- tryCatch(suppressWarnings(ape::read.evonet(text = enewick)),
                  error = function(e) NULL)
  degenerate <- is.null(net) || is.null(net$edge) ||
    nrow(net$edge) == 0L || length(net$tip.label) == 0L
  if (degenerate) {
    net <- ape::read.tree(text = enewick)
    if (!inherits(net, "evonet")) class(net) <- c("evonet", class(net))
  }
  if (is.null(net$reticulation)) net$reticulation <- matrix(integer(0), ncol = 2)
  net
}

# Legacy parse via ape::read.evonet + regex gammas (used only as a fallback if
# the native parser fails). Flags phantom tips it cannot resolve. Internal.
parse_network_legacy <- function(enewick, meta) {
  net <- parse_one_evonet(enewick)
  gammas <- parse_gammas(enewick)
  ntip <- length(net$tip.label)
  nl <- net$node.label
  hybrid_nodes <- integer(0)
  if (!is.null(nl)) {
    keep <- !is.na(nl) & nzchar(nl)
    if (any(keep)) hybrid_nodes <- stats::setNames(ntip + which(keep), nl[keep])
  }
  issues <- "parsed via legacy read.evonet fallback"
  n_empty <- sum(!nzchar(net$tip.label) | is.na(net$tip.label))
  if (n_empty > 0L) {
    issues <- c(issues, sprintf("%d unlabeled phantom tip(s)", n_empty))
  }
  structure(list(evonet = net, enewick = enewick, gammas = gammas,
                 hybrid_nodes = hybrid_nodes, meta = as.list(meta),
                 issues = issues),
            class = "anansi_network")
}

#' Parse one extended-Newick network into an `anansi_network`
#'
#' Uses anansi's native extended-Newick parser (which resolves nested/stacked
#' reticulations and hybrid-free trees), falling back to `ape::read.evonet` +
#' regex gamma parsing only if the native parser fails. The gamma table and
#' `#Hk`->node map come straight from the parse.
#'
#' @param enewick A single extended-Newick string.
#' @param meta Optional named list of per-network metadata (e.g. rank,
#'   log-likelihood) carried alongside the network.
#' @return An `anansi_network`: the `ape::evonet`, the raw string, the parsed
#'   gamma table, a `#Hk`->node map, parsing `issues` (character), and metadata.
#' @examples
#' net <- parse_network("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);")
#' net
#' @export
parse_network <- function(enewick, meta = list()) {
  stopifnot(is.character(enewick), length(enewick) == 1L)
  parsed <- tryCatch(enewick_to_anansi(enewick), error = function(e) NULL)
  ok <- !is.null(parsed) && length(parsed$evonet$tip.label) > 0L &&
    !any(is.na(parsed$evonet$tip.label) | !nzchar(parsed$evonet$tip.label))
  if (ok) {
    return(structure(list(evonet = parsed$evonet, enewick = enewick,
                          gammas = parsed$gammas,
                          hybrid_nodes = parsed$hybrid_nodes,
                          meta = as.list(meta), issues = parsed$issues),
                     class = "anansi_network"))
  }
  parse_network_legacy(enewick, meta)
}

# Resolve which CSV column holds the network strings. Internal.
detect_network_col <- function(df, network_col = NULL) {
  if (!is.null(network_col)) {
    if (!network_col %in% names(df)) {
      rlang::abort(sprintf("Column '%s' not found. Columns: %s",
                           network_col, paste(names(df), collapse = ", ")))
    }
    return(network_col)
  }
  candidates <- c("network", "newick", "enewick", "tree", "trees")
  for (cand in candidates) {
    hit <- names(df)[tolower(names(df)) == cand]
    if (length(hit)) return(hit[1])
  }
  rlang::abort(sprintf(
    "Could not auto-detect a network column among: %s. Pass `network_col=`.",
    paste(names(df), collapse = ", ")))
}

#' Path to a bundled example dataset
#'
#' Returns the path to a small example CSV shipped with anansi (40 inferred
#' networks over 15 lizard taxa, a subset of `Lacerta_agilis_inferred_networks`),
#' for use with [read_networks_csv].
#'
#' @param file Example file name (default `"lacerta_sample.csv"`).
#' @return The file path.
#' @examples
#' nets <- read_networks_csv(anansi_example())
#' nets
#' @export
anansi_example <- function(file = "lacerta_sample.csv") {
  system.file("extdata", file, package = "anansi", mustWork = TRUE)
}

#' Read a set of networks from a CSV of extended-Newick strings
#'
#' Reads a CSV in which one column holds extended-Newick network strings (named
#' `network`, `newick`, or `enewick`, auto-detected) and the remaining columns
#' are per-network metadata. Returns an [anansi_netset]. Networks whose taxon set
#' diverges from the modal set (e.g. parsing artifacts on nested reticulations)
#' are kept but flagged via `taxa_ok`; filter with `ns[ns$taxa_ok]`.
#'
#' @param path Path to the CSV file.
#' @param network_col Optional name of the network column (auto-detected if NULL).
#' @param validate If TRUE (default), warn when the taxon set is not fixed.
#' @param ... Passed to [utils::read.csv].
#' @return An [anansi_netset].
#' @export
read_networks_csv <- function(path, network_col = NULL, validate = TRUE, ...) {
  if (!file.exists(path)) rlang::abort(sprintf("File not found: %s", path))
  df <- utils::read.csv(path, stringsAsFactors = FALSE, ...)
  network_col <- detect_network_col(df, network_col)
  enwk <- df[[network_col]]
  meta_df <- df[, setdiff(names(df), network_col), drop = FALSE]
  networks <- lapply(seq_along(enwk), function(i) {
    parse_network(enwk[i], meta = as.list(meta_df[i, , drop = FALSE]))
  })
  anansi_netset(networks, meta = meta_df, validate = validate)
}
