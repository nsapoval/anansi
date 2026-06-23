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

# Backbone trees of the consistent networks, as a multiPhylo. Singleton
# (degree-2) nodes -- left behind at hybrid nodes once the reticulation edge is
# stripped -- are collapsed, since they do not change any clade but make ape's
# C-level consensus/prop.part routines segfault. Internal.
backbone_multiphylo <- function(ns) {
  trees <- lapply(ns$networks, function(n)
    ape::collapse.singles(backbone_tree(ensure_evonet(n))))
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
#' With `directed = FALSE`, events are instead keyed by the recipient clade plus
#' the *unordered* pair of its two parent lineages (the donor clade and the
#' backbone-sibling clade). This collapses reticulations that differ only by
#' which parent is the major (backbone) vs minor (reticulation) edge -- i.e.
#' direction-unstable events such as `((A,(B)#H1),(C,#H1))` and
#' `((C,(B)#H1),(A,#H1))` -- into a single event. Used by the undirected hybrid
#' layout (see [densinet]).
#'
#' @param x An [anansi_netset] (divergent networks are excluded).
#' @param directed If TRUE (default), key by directed donor -> recipient. If
#'   FALSE, key by recipient + unordered parent pair (direction-agnostic).
#' @return A data.frame ordered by decreasing frequency, with `key`, `count`,
#'   `freq`, `mean_gamma` and, for `directed = TRUE`, `donor`/`recipient`, or for
#'   `directed = FALSE`, `recipient`/`parent1`/`parent2`. Zero rows if there are
#'   no reticulations.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
#' reticulation_frequencies(anansi_netset(list(a, b)))
#' @export
reticulation_frequencies <- function(x, directed = TRUE) {
  ns <- as_consistent_netset(x)
  N <- length(ns$networks)
  empty <- if (directed) {
    data.frame(key = character(), donor = character(), recipient = character(),
               count = integer(), freq = numeric(), mean_gamma = numeric(),
               stringsAsFactors = FALSE)
  } else {
    data.frame(key = character(), recipient = character(), parent1 = character(),
               parent2 = character(), count = integer(), freq = numeric(),
               mean_gamma = numeric(), stringsAsFactors = FALSE)
  }
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
      if (directed) {
        rows[[length(rows) + 1L]] <- data.frame(
          key = paste(donor, "->", recip), donor = donor, recipient = recip,
          gamma = gamma, stringsAsFactors = FALSE)
      } else {
        # Backbone-sibling clade = tips below the hybrid's tree-parent, minus
        # the recipient itself. The hybrid's two parent lineages are then the
        # donor clade and this sibling clade; key by their unordered pair.
        tp <- ev$edge[ev$edge[, 2] == to, 1]
        sib <- if (length(tp)) paste(sort(setdiff(dt[[tp[1]]], dt[[to]])),
                                     collapse = ",") else ""
        pars <- sort(c(donor, sib))
        rows[[length(rows) + 1L]] <- data.frame(
          key = paste0(recip, " : ", paste(pars, collapse = " + ")),
          recipient = recip, parent1 = pars[1], parent2 = pars[2],
          gamma = gamma, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(empty)
  df <- do.call(rbind, rows)
  out <- do.call(rbind, lapply(split(df, df$key), function(d) {
    base <- data.frame(key = d$key[1], stringsAsFactors = FALSE)
    if (directed) {
      base$donor <- d$donor[1]; base$recipient <- d$recipient[1]
    } else {
      base$recipient <- d$recipient[1]
      base$parent1 <- d$parent1[1]; base$parent2 <- d$parent2[1]
    }
    base$count <- nrow(d)
    base$freq <- nrow(d) / N
    base$mean_gamma <- mean(d$gamma, na.rm = TRUE)
    base
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
#' @param directed Reticulation-event identity: `TRUE` (default) keys by directed
#'   donor -> recipient; `FALSE` keys direction-agnostically (recipient +
#'   unordered parent pair). See [reticulation_frequencies].
#' @return An `anansi_consensus`: `tree` (the consensus `phylo`), `support`
#'   (named by internal node), `reticulations` (data.frame), `directed`, `N`,
#'   `p`, `taxa`.
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
#' cons <- consensus_network(anansi_netset(list(a, b)))
#' cons$reticulations
#' @export
consensus_network <- function(x, p = 0.5, ret_support_min = 0, directed = TRUE) {
  ns <- as_consistent_netset(x)
  trees <- backbone_multiphylo(ns)
  N <- length(trees)
  cons <- ape::consensus(trees, p = p, rooted = TRUE)
  if (is.null(cons$node.label)) cons$node.label <- rep("", cons$Nnode)
  # Per-internal-node support (fraction of networks containing that clade).
  # rooted = TRUE matches the rooted consensus and avoids prop.clades returning
  # NA for clades it otherwise fails to match.
  ntip <- length(cons$tip.label)
  support <- stats::setNames(ape::prop.clades(cons, trees, rooted = TRUE) / N,
                             (ntip + 1L):(ntip + cons$Nnode))
  rets <- reticulation_frequencies(ns, directed = directed)
  if (nrow(rets)) rets <- rets[rets$freq >= ret_support_min, , drop = FALSE]
  heights <- consensus_node_heights(cons, trees)
  structure(list(tree = cons, support = support, heights = heights,
                 N = N, p = p, reticulations = rets, directed = directed,
                 taxa = ns$taxa),
            class = "anansi_consensus")
}

# Drawable segments for a consensus network on the shared, [0,1]-scaled frame.
# Returns list(backbone = df(x,y,xend,yend,support),
#              reticulations = df(x,y,xend,yend,freq,mean_gamma) | NULL).
# With a directed consensus (cons$directed), the reticulations are donor ->
# recipient edges (one per event) and the recipient keeps its solid backbone
# in-edge (the major/tree edge). With an undirected consensus they are two edges
# per event (recipient -> parent1, recipient -> parent2): this is a true
# consensus NETWORK, so the backbone in-edge into each hybrid recipient is
# omitted (the hybrid hangs from its two dashed edges only, like the single-
# network hybrid view) and the recipient end of those edges is anchored AT the
# node rather than nudged toward the now-removed parent. `ret_edge_frac` anchors
# each remaining endpoint a fraction of the way from its node toward the node's
# parent (0 = at the node, the legacy behavior). Internal.
consensus_segments <- function(cons, tip_order, mode = "cladogram",
                               ret_min = 0, ret_edge_frac = 0) {
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

  # Reticulation events (frequency-filtered) and recipient node ids are needed
  # both to draw the dashed edges below and -- for an undirected (hybrid)
  # consensus -- to drop the backbone in-edge into each hybrid recipient.
  rets <- cons$reticulations
  if (!is.null(rets) && nrow(rets)) rets <- rets[rets$freq >= ret_min, , drop = FALSE]
  has_rets <- !is.null(rets) && nrow(rets)
  directed <- is.null(cons$directed) || isTRUE(cons$directed)
  mrca <- function(tipstr) {
    tips <- strsplit(tipstr, ",", fixed = TRUE)[[1]]
    if (length(tips) == 1L) return(match(tips, tr$tip.label))
    m <- ape::getMRCA(tr, tips)
    if (is.null(m)) NA_integer_ else as.integer(m)
  }
  rc <- if (has_rets) vapply(rets$recipient, mrca, integer(1)) else integer(0)

  # A true consensus network: in undirected mode the hybrid node has only its
  # two dashed in-edges, so drop any backbone edge flowing INTO a recipient.
  keep_row <- rep(TRUE, nrow(E))
  if (!directed && length(rc)) {
    rc_drop <- rc[!is.na(rc)]
    if (length(rc_drop)) keep_row <- !(E[, 2] %in% rc_drop)
  }
  Ek <- E[keep_row, , drop = FALSE]
  childsup <- vapply(Ek[, 2], function(ch) {
    if (ch <= ntip) 1 else unname(cons$support[as.character(ch)])
  }, numeric(1))
  backbone <- data.frame(
    x = nx[as.character(Ek[, 1])], y = ny[as.character(Ek[, 1])],
    xend = nx[as.character(Ek[, 2])], yend = ny[as.character(Ek[, 2])],
    support = childsup, stringsAsFactors = FALSE)

  # Node coordinate, nudged `ret_edge_frac` toward the node's parent.
  anchor <- function(node) {
    if (is.na(node)) return(c(NA_real_, NA_real_))
    ax <- unname(nx[as.character(node)]); ay <- unname(ny[as.character(node)])
    pp <- E[E[, 2] == node, 1]
    if (!length(pp) || ret_edge_frac <= 0) return(c(ax, ay))
    px <- unname(nx[as.character(pp[1])]); py <- unname(ny[as.character(pp[1])])
    c(ax + ret_edge_frac * (px - ax), ay + ret_edge_frac * (py - ay))
  }
  # Raw node coordinate (no nudge); used for the recipient end of undirected
  # hybrid edges, whose backbone in-edge has been removed.
  node_xy <- function(node) {
    if (is.na(node)) return(c(NA_real_, NA_real_))
    c(unname(nx[as.character(node)]), unname(ny[as.character(node)]))
  }

  retseg <- NULL
  if (has_rets) {
    segs <- list()
    if (directed) {
      dn <- vapply(rets$donor, mrca, integer(1))
      for (i in seq_len(nrow(rets))) {
        if (is.na(dn[i]) || is.na(rc[i])) next
        a0 <- anchor(dn[i]); a1 <- anchor(rc[i])
        segs[[length(segs) + 1L]] <- data.frame(
          x = a0[1], y = a0[2], xend = a1[1], yend = a1[2],
          freq = rets$freq[i], mean_gamma = rets$mean_gamma[i],
          stringsAsFactors = FALSE)
      }
    } else {
      p1 <- vapply(rets$parent1, mrca, integer(1))
      p2 <- vapply(rets$parent2, mrca, integer(1))
      for (i in seq_len(nrow(rets))) {
        if (is.na(rc[i])) next
        ra <- node_xy(rc[i])
        for (pp in c(p1[i], p2[i])) {
          if (is.na(pp)) next
          pa <- anchor(pp)
          segs[[length(segs) + 1L]] <- data.frame(
            x = ra[1], y = ra[2], xend = pa[1], yend = pa[2],
            freq = rets$freq[i], mean_gamma = rets$mean_gamma[i],
            stringsAsFactors = FALSE)
        }
      }
    }
    if (length(segs)) retseg <- do.call(rbind, segs)
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
