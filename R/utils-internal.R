# utils-internal.R ----------------------------------------------------------
# Small internal helpers (not exported).

# Column names referenced inside ggplot2::aes() (non-standard evaluation), so
# R CMD check does not flag them as undefined globals.
utils::globalVariables(c("x", "y", "xend", "yend", ".net", "label", "support"))

# NULL-coalescing operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

# nrow() that treats NULL as 0.
nrow_or0 <- function(m) if (is.null(m)) 0L else nrow(m)

# Number of reticulations in an evonet/anansi_network.
n_reticulations <- function(x) {
  net <- if (inherits(x, "anansi_network")) x$evonet else x
  nrow_or0(net$reticulation)
}

# The backbone tree of an evonet: the tree edges only (reticulations dropped),
# keeping the original node numbering (so $reticulation indices stay valid).
# Built by stripping rather than as.phylo() to avoid S3-dispatch surprises.
# Internal.
backbone_tree <- function(ev) {
  tr <- ev
  tr$reticulation <- NULL
  class(tr) <- "phylo"
  tr
}

# Coerce to an evonet with a guaranteed (possibly empty) reticulation matrix.
# Accepts anansi_network, evonet, or phylo. Internal.
ensure_evonet <- function(x) {
  if (inherits(x, "anansi_network")) x <- x$evonet
  if (!inherits(x, "evonet")) {
    if (inherits(x, "phylo")) {
      class(x) <- c("evonet", class(x))
    } else {
      rlang::abort("Expected an anansi_network, evonet, or phylo.")
    }
  }
  if (is.null(x$reticulation)) x$reticulation <- matrix(integer(0), ncol = 2)
  x
}
