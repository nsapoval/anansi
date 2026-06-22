# taxa.R --------------------------------------------------------------------
#
# Taxon-set queries, consistency checks, and (Phase 5) restricting all networks
# to a subset of taxa. See docs/DESIGN.md (D7) for the pruning design and
# docs/DATA.md for the read.evonet taxon-divergence artifact.

#' Taxa of a network or network set
#' @param x An [anansi_network], [anansi_netset], or ape `evonet`/`phylo`.
#' @param ... Unused.
#' @return A character vector of taxon labels.
#' @export
network_taxa <- function(x, ...) UseMethod("network_taxa")

#' @rdname network_taxa
#' @export
network_taxa.anansi_network <- function(x, ...) x$evonet$tip.label

#' @rdname network_taxa
#' @export
network_taxa.anansi_netset <- function(x, ...) x$taxa

#' @rdname network_taxa
#' @export
network_taxa.phylo <- function(x, ...) x$tip.label

#' @rdname network_taxa
#' @export
network_taxa.default <- function(x, ...) {
  rlang::abort("No network_taxa() method for this object.")
}

#' Taxon-set consistency across a set of networks
#'
#' Compares every network's taxon set against the *modal* (most frequent) taxon
#' set in the collection and reports which networks diverge. Used to flag
#' parsing artifacts (e.g. phantom tips from nested reticulations) without
#' discarding data.
#'
#' @param x An [anansi_netset] or a list of [anansi_network].
#' @return A list with `ok` (logical, per network), `reference` (the modal taxon
#'   set), `n_bad`, and `bad` (indices of divergent networks).
#' @export
taxa_consistency <- function(x) {
  nets <- if (inherits(x, "anansi_netset")) x$networks else x
  if (!is.list(nets) || length(nets) == 0L) {
    rlang::abort("Expected an anansi_netset or a non-empty list of networks.")
  }
  keys <- vapply(nets, function(n)
    paste(sort(n$evonet$tip.label), collapse = "\t"), character(1))
  ref_key <- names(sort(table(keys), decreasing = TRUE))[1]
  ok <- unname(keys == ref_key)
  list(ok = ok,
       reference = strsplit(ref_key, "\t", fixed = TRUE)[[1]],
       n_bad = sum(!ok),
       bad = which(!ok))
}

# Reticulation (minor-edge) gamma for a recipient node, via the gamma table.
# Internal.
ret_gamma_of <- function(net, to) {
  if (length(net$hybrid_nodes) && !is.null(net$gammas)) {
    lbl <- names(net$hybrid_nodes)[match(to, net$hybrid_nodes)]
    if (length(lbl) && !is.na(lbl)) {
      gg <- net$gammas$ret_gamma[net$gammas$label == lbl]
      if (length(gg)) return(gg[1])
    }
  }
  NA_real_
}

# Reticulation-aware pruning of a single network. Internal; see restrict_taxa.
restrict_network <- function(net, keep) {
  ev <- ensure_evonet(net)
  tr <- backbone_tree(ev)
  keep <- intersect(keep, tr$tip.label)
  if (length(keep) < 2L) {
    rlang::abort("At least 2 of the requested taxa must be present in the network.")
  }

  # Record each reticulation by the surviving tips of its donor/recipient clades.
  dt <- descendant_tips(tr)
  R <- ev$reticulation
  specs <- list()
  if (nrow_or0(R) > 0L) {
    for (j in seq_len(nrow(R))) {
      specs[[j]] <- list(
        donor = intersect(dt[[R[j, 1]]], keep),
        recip = intersect(dt[[R[j, 2]]], keep),
        gamma = ret_gamma_of(net, R[j, 2]))
    }
  }

  # Prune the backbone (suppresses degree-2 nodes, renumbers).
  tr$node.label <- NULL
  tr2 <- ape::keep.tip(tr, keep)
  tr2$node.label <- NULL
  ntip2 <- length(tr2$tip.label)
  dt2 <- descendant_tips(tr2)

  mrca2 <- function(tips) {
    tips <- intersect(tips, tr2$tip.label)
    if (length(tips) == 0L) return(NA_integer_)
    if (length(tips) == 1L) return(match(tips, tr2$tip.label))
    as.integer(ape::getMRCA(tr2, tips))
  }

  # Re-locate surviving reticulations; drop degenerate ones.
  newR <- matrix(integer(0), ncol = 2)
  gam <- numeric(0)
  for (s in specs) {
    dn <- mrca2(s$donor); rc <- mrca2(s$recip)
    if (is.na(dn) || is.na(rc) || dn == rc) next
    dts <- dt2[[dn]]; rts <- dt2[[rc]]
    if (all(rts %in% dts) || all(dts %in% rts)) next   # one clade nests in other
    newR <- rbind(newR, c(dn, rc)); gam <- c(gam, s$gamma)
  }

  ev2 <- tr2
  if (nrow(newR) > 0L) {
    ev2$reticulation <- newR
    labels <- paste0("#H", seq_len(nrow(newR)))
    nl <- rep("", ev2$Nnode)
    for (k in seq_len(nrow(newR))) {
      rec <- newR[k, 2]
      if (rec > ntip2) nl[rec - ntip2] <- labels[k]
    }
    ev2$node.label <- nl
    hybrid_nodes <- stats::setNames(newR[, 2], labels)
    gammas <- data.frame(label = labels, tree_gamma = 1 - gam, ret_gamma = gam,
                         stringsAsFactors = FALSE)
  } else {
    ev2$reticulation <- matrix(integer(0), ncol = 2)
    hybrid_nodes <- integer(0)
    gammas <- data.frame(label = character(), tree_gamma = numeric(),
                         ret_gamma = numeric(), stringsAsFactors = FALSE)
  }
  class(ev2) <- c("evonet", "phylo")

  structure(list(evonet = ev2, enewick = NA_character_, gammas = gammas,
                 hybrid_nodes = hybrid_nodes, meta = net$meta,
                 issues = character(0)),
            class = "anansi_network")
}

#' Restrict networks to a subset of taxa (reticulation-aware pruning)
#'
#' Prunes every network to the taxa in `keep`. The backbone tree is pruned with
#' [ape::keep.tip] (suppressing degree-2 nodes); each reticulation is then
#' re-located by the MRCA of the surviving tips of its donor and recipient
#' clades, and dropped if either endpoint loses all its tips, the endpoints
#' collapse to the same node, or one clade nests within the other (see
#' docs/DESIGN.md D7).
#'
#' @param x An [anansi_netset] or a single [anansi_network].
#' @param keep Character vector of taxa to keep (>= 2 must be present).
#' @param validate Passed to [anansi_netset] when `x` is a set.
#' @return The same class as `x`, restricted to `keep`.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' ns <- anansi_netset(list(a, b))
#' restrict_taxa(ns, c("A", "B", "C"))
#' @export
restrict_taxa <- function(x, keep, validate = TRUE) {
  if (inherits(x, "anansi_network")) return(restrict_network(x, keep))
  if (!inherits(x, "anansi_netset")) {
    rlang::abort("restrict_taxa() needs an anansi_netset or anansi_network.")
  }
  pruned <- lapply(x$networks, function(n)
    tryCatch(restrict_network(n, keep), error = function(e) NULL))
  ok <- !vapply(pruned, is.null, logical(1))
  if (!any(ok)) rlang::abort("No networks could be restricted to the requested taxa.")
  if (!all(ok)) {
    warning(sprintf("Dropped %d network(s) lacking enough of the requested taxa.",
                    sum(!ok)), call. = FALSE)
  }
  meta <- if (!is.null(x$meta)) x$meta[ok, , drop = FALSE] else NULL
  anansi_netset(pruned[ok], meta = meta, validate = validate)
}

#' Validate that a set of networks shares one fixed taxon set
#'
#' @param x An [anansi_netset] or a list of [anansi_network].
#' @param error If TRUE abort on mismatch; if FALSE (default) warn and return
#'   FALSE.
#' @return TRUE invisibly if the taxon set is fixed; otherwise FALSE (or aborts).
#' @export
validate_taxon_set <- function(x, error = FALSE) {
  tc <- taxa_consistency(x)
  if (tc$n_bad > 0L) {
    msg <- sprintf(
      "%d of %d networks differ from the modal taxon set (e.g. index %s). See docs/DATA.md.",
      tc$n_bad, length(tc$ok), paste(utils::head(tc$bad, 10L), collapse = ", "))
    if (error) rlang::abort(msg)
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
