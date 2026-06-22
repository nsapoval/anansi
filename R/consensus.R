# consensus.R ---------------------------------------------------------------
#
# Consensus / discrepancy summaries across a set of networks: backbone clade
# frequencies, reticulation-event frequencies, and a consensus ("root-canal")
# network. See docs/DESIGN.md (D4, D5).

# The consistent sub-set of a netset (taxa_ok networks only). Internal.
as_consistent_netset <- function(x) {
  if (!inherits(x, "anansi_netset")) rlang::abort("Expected an anansi_netset.")
  if (all(x$taxa_ok)) x else x[x$taxa_ok]
}

# Descendant tip labels of every node of a tree (list indexed by node). Internal.
descendant_tips <- function(tr) {
  ntip <- length(tr$tip.label)
  maxnode <- ntip + tr$Nnode
  res <- vector("list", maxnode)
  for (i in seq_len(ntip)) res[[i]] <- tr$tip.label[i]
  epo <- ape::reorder.phylo(tr, "postorder")$edge
  for (k in seq_len(nrow(epo))) {
    par <- epo[k, 1]; ch <- epo[k, 2]
    res[[par]] <- c(res[[par]], res[[ch]])
  }
  res
}

# Backbone trees of the consistent networks, as a multiPhylo. Internal.
backbone_multiphylo <- function(ns) {
  trees <- lapply(ns$networks, function(n) backbone_tree(ensure_evonet(n)))
  class(trees) <- "multiPhylo"
  trees
}

# Mean root-to-node distance for each node of a consensus tree, averaged over
# the source trees' branch lengths (for phylogram consensus). NA-safe. Internal.
consensus_node_heights <- function(cons_tree, trees) {
  dt <- descendant_tips(cons_tree)
  nn <- length(cons_tree$tip.label) + cons_tree$Nnode
  depthlist <- lapply(trees, function(tr)
    if (is.null(tr$edge.length)) NULL else ape::node.depth.edgelength(tr))
  vapply(seq_len(nn), function(nd) {
    tips <- dt[[nd]]
    vals <- vapply(seq_along(trees), function(i) {
      d <- depthlist[[i]]
      if (is.null(d)) return(NA_real_)
      tr <- trees[[i]]
      node <- if (length(tips) == 1L) match(tips, tr$tip.label) else {
        m <- ape::getMRCA(tr, tips); if (is.null(m)) NA_integer_ else m
      }
      if (is.na(node)) NA_real_ else d[node]
    }, numeric(1))
    mean(vals, na.rm = TRUE)
  }, numeric(1))
}

#' Backbone clade frequencies across a set of networks
#'
#' Tallies, across the backbone trees, how often each clade (the tip set below an
#' internal node) appears. This is the majority-rule consensus building block.
#'
#' @param x An [anansi_netset] (divergent networks are excluded).
#' @return A data.frame with `clade` (comma-separated taxa), `size`, `count` and
#'   `freq` (count / number of networks), ordered by decreasing frequency.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' clade_frequencies(anansi_netset(list(a, b)))
#' @export
clade_frequencies <- function(x) {
  ns <- as_consistent_netset(x)
  trees <- backbone_multiphylo(ns)
  N <- length(trees)
  pp <- ape::prop.part(trees)
  cnt <- attr(pp, "number")
  labs <- attr(pp, "labels")
  df <- data.frame(
    clade = vapply(pp, function(idx) paste(sort(labs[idx]), collapse = ","),
                   character(1)),
    size  = vapply(pp, length, integer(1)),
    count = as.integer(cnt),
    freq  = as.numeric(cnt) / N,
    stringsAsFactors = FALSE)
  df <- df[order(-df$freq, -df$size), ]
  rownames(df) <- NULL
  df
}

#' Reticulation-event frequencies across a set of networks
#'
#' Keys each reticulation by the pair of backbone clades it connects -- the
#' donor clade (tips below the source node) and the recipient clade (tips below
#' the hybrid node) -- so equivalent reticulations across networks are tallied
#' together (see docs/DESIGN.md D4). Reports frequency and mean inheritance
#' probability (gamma) per event.
#'
#' @param x An [anansi_netset] (divergent networks are excluded).
#' @return A data.frame with `key`, `donor`, `recipient`, `count`, `freq` and
#'   `mean_gamma`, ordered by decreasing frequency. Zero rows if there are no
#'   reticulations.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
#' reticulation_frequencies(anansi_netset(list(a, b)))
#' @export
reticulation_frequencies <- function(x) {
  ns <- as_consistent_netset(x)
  N <- length(ns$networks)
  empty <- data.frame(key = character(), donor = character(),
                      recipient = character(), count = integer(),
                      freq = numeric(), mean_gamma = numeric(),
                      stringsAsFactors = FALSE)
  rows <- list()
  for (net in ns$networks) {
    ev <- ensure_evonet(net)
    R <- ev$reticulation
    if (nrow_or0(R) == 0L) next
    dt <- descendant_tips(backbone_tree(ev))
    for (j in seq_len(nrow(R))) {
      from <- R[j, 1]; to <- R[j, 2]
      donor <- paste(sort(dt[[from]]), collapse = ",")
      recip <- paste(sort(dt[[to]]), collapse = ",")
      gamma <- NA_real_
      if (length(net$hybrid_nodes) && !is.null(net$gammas)) {
        lbl <- names(net$hybrid_nodes)[match(to, net$hybrid_nodes)]
        if (length(lbl) && !is.na(lbl)) {
          gg <- net$gammas$ret_gamma[net$gammas$label == lbl]
          if (length(gg)) gamma <- gg[1]
        }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        key = paste(donor, "->", recip), donor = donor, recipient = recip,
        gamma = gamma, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(empty)
  df <- do.call(rbind, rows)
  out <- do.call(rbind, lapply(split(df, df$key), function(d) {
    data.frame(key = d$key[1], donor = d$donor[1], recipient = d$recipient[1],
               count = nrow(d), freq = nrow(d) / N,
               mean_gamma = mean(d$gamma, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  out <- out[order(-out$freq), ]
  rownames(out) <- NULL
  out
}

#' Consensus ("root-canal") network for a set of networks
#'
#' Builds a majority-rule consensus backbone tree (via [ape::consensus]) with
#' per-clade support, plus the frequent reticulation events
#' ([reticulation_frequencies]). This is the representative topology drawn on top
#' of the [densinet] cloud.
#'
#' @param x An [anansi_netset] (divergent networks are excluded).
#' @param p Majority threshold for the consensus tree (default 0.5).
#' @param ret_support_min Keep only reticulation events with `freq >=` this
#'   (default 0, keep all).
#' @return An `anansi_consensus`: `tree` (the consensus `phylo`), `support`
#'   (named by internal node), `reticulations` (data.frame), `N`, `p`, `taxa`.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
#' cons <- consensus_network(anansi_netset(list(a, b)))
#' cons$reticulations
#' @export
consensus_network <- function(x, p = 0.5, ret_support_min = 0) {
  ns <- as_consistent_netset(x)
  trees <- backbone_multiphylo(ns)
  N <- length(trees)
  cons <- ape::consensus(trees, p = p, rooted = TRUE)
  if (is.null(cons$node.label)) cons$node.label <- rep("", cons$Nnode)
  # Per-internal-node support (fraction of networks containing that clade).
  ntip <- length(cons$tip.label)
  support <- stats::setNames(ape::prop.clades(cons, trees) / N,
                             (ntip + 1L):(ntip + cons$Nnode))
  rets <- reticulation_frequencies(ns)
  if (nrow(rets)) rets <- rets[rets$freq >= ret_support_min, , drop = FALSE]
  heights <- consensus_node_heights(cons, trees)
  structure(list(tree = cons, support = support, heights = heights,
                 N = N, p = p, reticulations = rets, taxa = ns$taxa),
            class = "anansi_consensus")
}

# Drawable segments for a consensus network on the shared, [0,1]-scaled frame.
# Returns list(backbone = df(x,y,xend,yend,support),
#              reticulations = df(x,y,xend,yend,freq,mean_gamma) | NULL). Internal.
consensus_segments <- function(cons, tip_order, mode = "cladogram", ret_min = 0) {
  tr <- cons$tree
  # y is mode-independent; compute via cladogram (the consensus tree has no
  # branch lengths, so phylogram x comes from cons$heights instead).
  nc <- compute_node_coords(tr, tip_order, "cladogram")
  xraw <- if (mode == "phylogram" && !is.null(cons$heights)) {
    cons$heights[nc$node]
  } else {
    nc$x
  }
  xr <- range(xraw, na.rm = TRUE)
  rg <- if (diff(xr) > 0) diff(xr) else 1
  nx <- stats::setNames((xraw - xr[1]) / rg, nc$node)
  ny <- stats::setNames(nc$y, nc$node)
  ntip <- length(tr$tip.label)

  E <- tr$edge
  childsup <- vapply(E[, 2], function(ch) {
    if (ch <= ntip) 1 else unname(cons$support[as.character(ch)])
  }, numeric(1))
  backbone <- data.frame(
    x = nx[as.character(E[, 1])], y = ny[as.character(E[, 1])],
    xend = nx[as.character(E[, 2])], yend = ny[as.character(E[, 2])],
    support = childsup, stringsAsFactors = FALSE)

  retseg <- NULL
  rets <- cons$reticulations
  if (nrow(rets)) rets <- rets[rets$freq >= ret_min, , drop = FALSE]
  if (nrow(rets)) {
    mrca <- function(tipstr) {
      tips <- strsplit(tipstr, ",", fixed = TRUE)[[1]]
      if (length(tips) == 1L) return(match(tips, tr$tip.label))
      m <- ape::getMRCA(tr, tips)
      if (is.null(m)) NA_integer_ else as.integer(m)
    }
    dn <- vapply(rets$donor, mrca, integer(1))
    rc <- vapply(rets$recipient, mrca, integer(1))
    ok <- !is.na(dn) & !is.na(rc)
    if (any(ok)) {
      retseg <- data.frame(
        x = nx[as.character(dn[ok])], y = ny[as.character(dn[ok])],
        xend = nx[as.character(rc[ok])], yend = ny[as.character(rc[ok])],
        freq = rets$freq[ok], mean_gamma = rets$mean_gamma[ok],
        stringsAsFactors = FALSE)
    }
  }
  list(backbone = backbone, reticulations = retseg)
}

#' Keep the top-N networks of a set by a metadata column
#'
#' @param x An [anansi_netset].
#' @param n Number of networks to keep.
#' @param by Metadata column to rank by; if NULL, the first available of
#'   `log_pseudo_likelihood`, `log_probability`, `network_rank`, `rank` is used.
#'   Likelihood columns are ranked high-to-low; `*rank` columns low-to-high.
#' @return An [anansi_netset] with the top-N networks.
#' @examples
#' ns <- anansi_netset(
#'   list(parse_network("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);"),
#'        parse_network("((A:1,(B:1)#H1:1::0.6):1,(C:1,#H1:1::0.4):1);")),
#'   meta = data.frame(score = c(-10, -5)))
#' top_networks(ns, 1, by = "score")
#' @export
top_networks <- function(x, n, by = NULL) {
  if (!inherits(x, "anansi_netset")) rlang::abort("Expected an anansi_netset.")
  if (is.null(x$meta)) rlang::abort("No metadata to rank by.")
  if (is.null(by)) {
    cand <- intersect(c("log_pseudo_likelihood", "log_probability",
                        "network_rank", "rank"), names(x$meta))
    if (!length(cand)) rlang::abort("Could not auto-pick a ranking column; set `by`.")
    by <- cand[1]
  }
  if (!by %in% names(x$meta)) rlang::abort(sprintf("No metadata column '%s'.", by))
  decreasing <- !grepl("rank", by, ignore.case = TRUE)
  ord <- order(x$meta[[by]], decreasing = decreasing)
  x[ord[seq_len(min(n, length(ord)))]]
}

#' @export
print.anansi_consensus <- function(x, ...) {
  cat("<anansi_consensus>\n")
  cat("  from networks:", x$N, "\n")
  cat("  consensus p:  ", x$p, "\n")
  cat("  taxa:         ", length(x$taxa), "\n")
  cat("  clades:       ", x$tree$Nnode, " (support ",
      sprintf("%.2f-%.2f", min(x$support), max(x$support)), ")\n", sep = "")
  cat("  reticulation events:", nrow(x$reticulations), "\n")
  invisible(x)
}
