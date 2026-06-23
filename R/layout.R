# layout.R ------------------------------------------------------------------
#
# The shared-coordinate layout engine: the heart of the overlay (see
# docs/DESIGN.md D1-D3). Each network is decomposed into a backbone tree (which
# fixes tip order and node depths) plus reticulation edges drawn on the same
# frame. All networks are laid out on ONE tip order so they overlay.

# --- tip ordering ----------------------------------------------------------

# Tip labels of one network in its intrinsic (plotting) top-to-bottom order,
# i.e. the order tips appear in a cladewise traversal of the backbone. Internal.
network_tip_order <- function(net) {
  tr <- ape::reorder.phylo(backbone_tree(ensure_evonet(net)), "cladewise")
  ntip <- length(tr$tip.label)
  tip_nodes <- tr$edge[tr$edge[, 2] <= ntip, 2]
  tr$tip.label[tip_nodes]
}

# Mean topological (unit-branch) cophenetic distance over networks whose taxon
# set equals `ref`. Internal.
avg_tip_distance <- function(nets, ref) {
  acc <- matrix(0, length(ref), length(ref), dimnames = list(ref, ref))
  cnt <- 0L
  for (n in nets) {
    ev <- ensure_evonet(n)
    if (!setequal(ev$tip.label, ref)) next
    tr <- backbone_tree(ev)
    tr$edge.length <- rep(1, nrow(tr$edge))
    co <- ape::cophenetic.phylo(tr)[ref, ref]
    acc <- acc + co
    cnt <- cnt + 1L
  }
  if (cnt == 0L) rlang::abort("No networks with the reference taxon set.")
  acc / cnt
}

# 1-D MDS ordering of a distance matrix. Internal.
order_mds <- function(D) {
  mds <- stats::cmdscale(stats::as.dist(D), k = 1)
  rownames(D)[order(mds[, 1])]
}

# DensiTree's greedy closest-leaf ordering. Internal.
order_closest_leaf <- function(D) {
  labs <- rownames(D)
  n <- length(labs)
  if (n <= 2L) return(labs)
  diag(D) <- Inf
  ij <- which(D == min(D), arr.ind = TRUE)[1, ]
  ord <- c(ij[1], ij[2])
  used <- c(ij[1], ij[2])
  while (length(ord) < n) {
    cand <- setdiff(seq_len(n), used)
    dl <- D[ord[1], cand]
    dr <- D[ord[length(ord)], cand]
    if (min(dl) <= min(dr)) {
      add <- cand[which.min(dl)]; ord <- c(add, ord)
    } else {
      add <- cand[which.min(dr)]; ord <- c(ord, add)
    }
    used <- c(used, add)
  }
  labs[ord]
}

# Most frequent exact tip order across consistent networks. Internal.
order_mode <- function(nets, ref) {
  ords <- lapply(nets, function(n) tryCatch(network_tip_order(n), error = function(e) NULL))
  ords <- Filter(function(o) !is.null(o) && setequal(o, ref), ords)
  if (!length(ords)) return(ref)
  keys <- vapply(ords, paste, character(1), collapse = "\t")
  best <- names(sort(table(keys), decreasing = TRUE))[1]
  strsplit(best, "\t", fixed = TRUE)[[1]]
}

# Move outgroup taxa to one end of `order`, kept contiguous and in their
# current relative order. position "top" -> end of the vector (rendered at the
# top, since y = match(label, order)); "bottom" -> start. Internal.
pin_outgroup <- function(order, outgroup, position = c("top", "bottom")) {
  position <- match.arg(position)
  outgroup <- as.character(outgroup)
  miss <- setdiff(outgroup, order)
  if (length(miss)) {
    rlang::abort(sprintf("outgroup taxa not in the taxon set: %s",
                         paste(miss, collapse = ", ")))
  }
  og  <- order[order %in% outgroup]    # outgroup, current relative order
  ing <- order[!order %in% outgroup]   # ingroup, order preserved
  if (position == "top") c(ing, og) else c(og, ing)
}

# Rotate each internal node's children to follow `seed` (children ordered by the
# mean seed-rank of their descendant tips), then read off the resulting tip order
# by depth-first traversal. Every clade is contiguous by construction, so the
# order matches `seed` as closely as the topology allows with no crossings. Pure
# R (avoids ape::rotateConstr, whose C routine can crash on some consensus
# trees) and handles polytomies. Internal.
order_tree_to_seed <- function(tr, seed) {
  ntip <- length(tr$tip.label)
  rank <- match(tr$tip.label, seed)               # seed rank by tip index
  child_by_parent <- split(tr$edge[, 2], tr$edge[, 1])
  root <- setdiff(tr$edge[, 1], tr$edge[, 2])[1]   # node that is never a child
  meanrank <- function(node) {
    if (node <= ntip) return(rank[node])
    mean(vapply(child_by_parent[[as.character(node)]], meanrank, numeric(1)))
  }
  collect <- function(node) {
    if (node <= ntip) return(tr$tip.label[node])
    kids <- child_by_parent[[as.character(node)]]
    kids <- kids[order(vapply(kids, meanrank, numeric(1)))]
    unlist(lapply(kids, collect), use.names = FALSE)
  }
  collect(root)
}

# Snap a seed tip order to the majority-rule consensus topology: build the
# consensus of the (consistent) backbone trees and rotate it to follow `seed`
# (see order_tree_to_seed), so every consensus clade is contiguous and the
# consensus backbone has no crossings. Internal.
snap_order_to_consensus <- function(nets, ref, seed, p = 0.5) {
  ok <- vapply(nets, function(n) setequal(ensure_evonet(n)$tip.label, ref),
               logical(1))
  nets <- nets[ok]
  if (!length(nets)) return(seed)
  # collapse.singles: drop unary backbone nodes that crash ape's consensus().
  trees <- lapply(nets, function(n)
    ape::collapse.singles(backbone_tree(ensure_evonet(n))))
  class(trees) <- "multiPhylo"
  cons <- ape::consensus(trees, p = p, rooted = TRUE)
  order_tree_to_seed(cons, seed)
}

#' Choose a shared tip ordering for a set of networks
#'
#' Derives a single linear tip order (applied to every network so they overlay)
#' from the backbone trees. See docs/DESIGN.md (D3).
#'
#' @param x An [anansi_netset] or a list of [anansi_network].
#' @param method One of `"mode"` (most frequent tip order; default), `"mds"`
#'   (1-D MDS of mean leaf distance), `"closest_leaf"` (DensiTree greedy
#'   heuristic), `"first"` (order of the first network), or `"manual"`.
#' @param order For `method = "manual"`, the explicit tip order (character).
#' @param outgroup Optional character vector of taxa to pin to one end of the
#'   order (kept contiguous, in their existing relative order).
#' @param outgroup_position Where to pin the `outgroup`: `"top"` (default) or
#'   `"bottom"` of the figure.
#' @param snap_to_consensus If TRUE, rotate the majority-rule consensus tree to
#'   follow the chosen order as closely as its topology allows (each node's
#'   children ordered by the mean rank of their descendant tips) and return that
#'   order instead, so every consensus clade is contiguous (no crossings in the
#'   consensus backbone). Default FALSE.
#' @param consensus_p Majority threshold for the consensus tree used by
#'   `snap_to_consensus` (default 0.5).
#' @return A character vector of taxa in the chosen order.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' consensus_tip_order(anansi_netset(list(a, b)), method = "mds")
#' @export
consensus_tip_order <- function(x,
                                method = c("mode", "mds", "closest_leaf",
                                           "first", "manual"),
                                order = NULL, outgroup = NULL,
                                outgroup_position = c("top", "bottom"),
                                snap_to_consensus = FALSE, consensus_p = 0.5) {
  method <- match.arg(method)
  outgroup_position <- match.arg(outgroup_position)
  nets <- if (inherits(x, "anansi_netset")) x$networks else x

  if (method == "manual") {
    if (is.null(order)) rlang::abort("method = 'manual' requires `order`.")
    base <- as.character(order)
    ref  <- if (inherits(x, "anansi_netset")) x$taxa else sort(base)
  } else {
    if (!is.list(nets) || length(nets) == 0L) {
      rlang::abort("Expected an anansi_netset or a non-empty list of networks.")
    }
    ref <- if (inherits(x, "anansi_netset")) x$taxa
           else sort(ensure_evonet(nets[[1]])$tip.label)
    base <- switch(method,
      first = network_tip_order(nets[[1]]),
      mode  = order_mode(nets, ref),
      mds   = order_mds(avg_tip_distance(nets, ref)),
      closest_leaf = order_closest_leaf(avg_tip_distance(nets, ref)))
  }

  if (!is.null(outgroup)) base <- pin_outgroup(base, outgroup, outgroup_position)

  if (snap_to_consensus) {
    if (!is.list(nets) || length(nets) == 0L) {
      rlang::abort("snap_to_consensus requires an anansi_netset or non-empty list of networks.")
    }
    base <- tryCatch(snap_order_to_consensus(nets, ref, base, p = consensus_p),
                     error = function(e) base)
  }
  base
}

# --- per-network layout ----------------------------------------------------

# Per-node coordinates of one network's backbone on a shared tip order.
# Returns a data.frame(node, x, y, is_tip, label). Internal.
compute_node_coords <- function(ev, tip_order, mode = c("cladogram", "phylogram")) {
  mode <- match.arg(mode)
  tr <- backbone_tree(ensure_evonet(ev))
  ntip <- length(tr$tip.label)
  maxnode <- ntip + tr$Nnode
  edge <- tr$edge

  # Height above tips (x for cladogram; also orders the y pass children-first).
  h <- rep(0L, maxnode)
  epo <- ape::reorder.phylo(tr, "postorder")$edge
  for (k in seq_len(nrow(epo))) {
    par <- epo[k, 1]; ch <- epo[k, 2]
    if (h[par] < h[ch] + 1L) h[par] <- h[ch] + 1L
  }
  xvec <- if (mode == "cladogram") max(h) - h else ape::node.depth.edgelength(tr)

  yvec <- rep(NA_real_, maxnode)
  yvec[seq_len(ntip)] <- match(tr$tip.label, tip_order)
  child_by_parent <- split(edge[, 2], edge[, 1])
  internal <- (ntip + 1L):maxnode
  internal <- internal[order(h[internal])]
  for (nd in internal) {
    kids <- child_by_parent[[as.character(nd)]]
    if (length(kids)) yvec[nd] <- mean(yvec[kids])
  }

  lab_vec <- rep(NA_character_, maxnode)
  lab_vec[seq_len(ntip)] <- tr$tip.label
  if (!is.null(tr$node.label) && length(tr$node.label) == tr$Nnode) {
    lab_vec[(ntip + 1L):maxnode] <- tr$node.label
  }
  data.frame(node = seq_len(maxnode), x = xvec, y = yvec,
             is_tip = seq_len(maxnode) <= ntip, label = lab_vec,
             stringsAsFactors = FALSE)
}

# Rescale a pair of x-vectors jointly to [0, 1]. Internal.
scale_x01 <- function(x, xend) {
  vals <- c(x, xend)
  mn <- min(vals, na.rm = TRUE); mx <- max(vals, na.rm = TRUE)
  rg <- if (mx > mn) mx - mn else 1
  list(x = (x - mn) / rg, xend = (xend - mn) / rg)
}

#' Lay out one network on a shared tip order
#'
#' Computes node coordinates on a fixed tip order and returns a tidy table of
#' drawable segments (backbone tree edges and reticulation edges). Tips are
#' pinned to their rank in `tip_order` (shared across networks); internal-node y
#' is the mean of its children; x is the backbone node depth (tip-aligned for
#' `"cladogram"`, branch-length-scaled for `"phylogram"`). See docs/DESIGN.md
#' (D1, D2).
#'
#' @param net An [anansi_network] or ape `evonet`/`phylo`.
#' @param tip_order Character vector giving the shared tip order.
#' @param mode `"cladogram"` (default, uniform depths) or `"phylogram"`
#'   (branch-length scaled).
#' @return A data.frame of segments with columns `from`, `to`, `x`, `y`, `xend`,
#'   `yend`, `kind` (`"tree"`/`"reticulation"`), `is_tip`, `label`.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' seg <- layout_network(a, tip_order = c("A", "B", "C", "D"))
#' head(seg)
#' @export
layout_network <- function(net, tip_order, mode = c("cladogram", "phylogram")) {
  mode <- match.arg(mode)
  ev <- ensure_evonet(net)
  nc <- compute_node_coords(ev, tip_order, mode)

  nm <- as.character(nc$node)
  xpos <- stats::setNames(nc$x, nm)
  ypos <- stats::setNames(nc$y, nm)
  lab  <- stats::setNames(nc$label, nm)
  tipf <- stats::setNames(nc$is_tip, nm)

  seg_from_pairs <- function(parent, child, kind) {
    if (length(parent) == 0L) return(NULL)
    data.frame(
      from = parent, to = child,
      x    = xpos[as.character(parent)], y    = ypos[as.character(parent)],
      xend = xpos[as.character(child)],  yend = ypos[as.character(child)],
      kind = kind, stringsAsFactors = FALSE)
  }

  seg <- seg_from_pairs(ev$edge[, 1], ev$edge[, 2], "tree")
  R <- ev$reticulation
  if (nrow_or0(R) > 0L) {
    seg <- rbind(seg, seg_from_pairs(R[, 1], R[, 2], "reticulation"))
  }

  seg$is_tip <- unname(tipf[as.character(seg$to)])
  seg$label  <- ifelse(seg$is_tip, unname(lab[as.character(seg$to)]), NA_character_)
  rownames(seg) <- NULL
  seg
}

#' Lay out a whole set of networks on one shared tip order
#'
#' Applies [layout_network] to every network in the set using a common tip
#' order, and (by default) normalizes each network's x to `[0, 1]` so backbones
#' overlay. Networks whose taxon set diverges (`taxa_ok = FALSE`) are dropped
#' with a warning by default. See docs/DESIGN.md (D2).
#'
#' @param x An [anansi_netset].
#' @param tip_order Optional explicit tip order; if NULL it is computed via
#'   [consensus_tip_order] with `method`.
#' @param method Tip-ordering method passed to [consensus_tip_order].
#' @param mode `"cladogram"` (default) or `"phylogram"`.
#' @param scale_x If TRUE (default), rescale each network's x to `[0, 1]`
#'   (root = 0, tips = 1) so tips align across networks.
#' @param jitter Max per-network x-offset (evenly spread over `[-jitter, jitter]`)
#'   to separate exactly overlapping edges; 0 (default) disables it.
#' @param consistent_only If TRUE (default), drop networks with `taxa_ok = FALSE`.
#' @param outgroup,outgroup_position,snap_to_consensus,consensus_p Forwarded to
#'   [consensus_tip_order] when `tip_order` is NULL (outgroup pinning and
#'   consensus-topology crossing reduction).
#' @return A data.frame of segments (as [layout_network]) with an added integer
#'   `.net` column; the chosen order is attached as `attr(., "tip_order")`.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' seg <- layout_netset(anansi_netset(list(a, b)), method = "first")
#' attr(seg, "tip_order")
#' @export
layout_netset <- function(x, tip_order = NULL, method = "mode",
                          mode = "cladogram", scale_x = TRUE, jitter = 0,
                          consistent_only = TRUE, outgroup = NULL,
                          outgroup_position = c("top", "bottom"),
                          snap_to_consensus = FALSE, consensus_p = 0.5) {
  if (!inherits(x, "anansi_netset")) rlang::abort("layout_netset() needs an anansi_netset.")
  ns <- x
  if (consistent_only && !all(ns$taxa_ok)) {
    warning(sprintf("Dropping %d network(s) with divergent taxon sets (taxa_ok = FALSE).",
                    sum(!ns$taxa_ok)), call. = FALSE)
    ns <- ns[ns$taxa_ok]
  }
  if (length(ns$networks) == 0L) rlang::abort("No networks left to lay out.")
  if (is.null(tip_order)) {
    tip_order <- consensus_tip_order(ns, method = method, outgroup = outgroup,
                                     outgroup_position = outgroup_position,
                                     snap_to_consensus = snap_to_consensus,
                                     consensus_p = consensus_p)
  }

  segs <- lapply(seq_along(ns$networks), function(i) {
    s <- layout_network(ns$networks[[i]], tip_order = tip_order, mode = mode)
    s$.net <- i
    s
  })
  out <- do.call(rbind, segs)

  if (scale_x) {
    for (i in unique(out$.net)) {
      sel <- out$.net == i
      s <- scale_x01(out$x[sel], out$xend[sel])
      out$x[sel] <- s$x
      out$xend[sel] <- s$xend
    }
  }

  if (jitter > 0) {
    nets <- unique(out$.net)
    nn <- length(nets)
    off <- if (nn > 1) jitter * (2 * (seq_len(nn) - 1) / (nn - 1) - 1)
           else stats::setNames(0, nets)
    off <- stats::setNames(off, as.character(nets))
    o <- off[as.character(out$.net)]
    out$x <- out$x + o
    out$xend <- out$xend + o
  }

  attr(out, "tip_order") <- tip_order
  out
}
