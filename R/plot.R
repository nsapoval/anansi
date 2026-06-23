# plot.R --------------------------------------------------------------------
#
# Visualization. Phase 1: plot_network() (single network, parity with the
# reference plot_network.R). Phase 3: densinet() overlay. Phase 4: consensus
# encoding. See docs/DESIGN.md (section 7) and docs/WORKPLAN.md.

# Build a default title from per-network metadata. Internal.
default_title <- function(meta) {
  if (length(meta) == 0L) return(NULL)
  parts <- character(0)
  if (!is.null(meta$dataset)) parts <- c(parts, paste0("Dataset: ", meta$dataset))
  rc <- meta$reticulations %||% meta$num_reticulations
  if (!is.null(rc)) parts <- c(parts, paste0("Reticulations: ", rc))
  lp <- meta$log_probability %||% meta$log_pseudo_likelihood
  if (!is.null(lp)) parts <- c(parts, paste0("Log-prob: ", round(as.numeric(lp), 2)))
  if (length(parts) == 0L) return(NULL)
  paste(parts, collapse = "  |  ")
}

# "gamma = 0.78" with a real Greek gamma, built at runtime so the source stays
# ASCII (portable-package requirement). Internal.
gamma_label <- function(value) {
  sprintf("%s = %.2f", intToUtf8(0x03B3L), value)
}

#' Plot a single phylogenetic network
#'
#' Renders one network with `tanggle::ggevonet()` and (for an [anansi_network])
#' annotates each hybrid edge with its gamma inheritance probability: the minor
#' (reticulation) edge in blue and the major (tree) edge in red. This is the
#' framework version of the reference `plot_network.R`.
#'
#' @param x An [anansi_network] (preferred, enables gamma labels) or an ape
#'   `evonet`/`phylo`.
#' @param cladogram If TRUE (default) drop branch lengths so node depths are
#'   uniform (topology-only view).
#' @param layout `ggevonet` layout, e.g. "slanted" (default) or "rectangular".
#' @param tip_size Tip-label text size.
#' @param show_gamma If TRUE (default) annotate gamma values (needs an
#'   `anansi_network`).
#' @param title Plot title; if NULL a title is derived from metadata.
#' @return A `ggplot` object.
#' @examples
#' net <- parse_network("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);")
#' p <- plot_network(net)   # a ggplot object; print(p) to draw it
#' @export
plot_network <- function(x, cladogram = TRUE, layout = "slanted",
                         tip_size = 3, show_gamma = TRUE, title = NULL) {
  if (inherits(x, "anansi_network")) {
    net <- x$evonet; gammas <- x$gammas
    hybrid_nodes <- x$hybrid_nodes; meta <- x$meta
  } else if (inherits(x, c("evonet", "phylo"))) {
    net <- x; gammas <- NULL; hybrid_nodes <- NULL; meta <- list()
  } else {
    rlang::abort("plot_network() needs an anansi_network or an ape evonet/phylo.")
  }

  if (cladogram) net$edge.length <- NULL

  p <- tanggle::ggevonet(net, layout = layout) +
    ggtree::geom_tiplab(size = tip_size)

  ttl <- title %||% default_title(meta)
  if (!is.null(ttl)) p <- p + ggplot2::ggtitle(ttl)

  ret <- net$reticulation
  have_gamma <- show_gamma && nrow_or0(ret) > 0L &&
    !is.null(gammas) && nrow(gammas) > 0L && length(hybrid_nodes) > 0L

  if (have_gamma) {
    d <- p$data
    for (i in seq_len(nrow(ret))) {
      from <- ret[i, 1]; to <- ret[i, 2]
      lbl <- names(hybrid_nodes)[match(to, hybrid_nodes)]
      g <- gammas[gammas$label == lbl, , drop = FALSE]
      if (nrow(g) == 0L) next

      # Minor (reticulation) edge: donor -> hybrid (blue).
      xm <- (d$x[d$node == from] + d$x[d$node == to]) / 2
      ym <- (d$y[d$node == from] + d$y[d$node == to]) / 2
      if (length(xm) && !is.na(g$ret_gamma)) {
        p <- p + ggplot2::annotate("text", x = xm, y = ym + 0.3,
          label = gamma_label(g$ret_gamma),
          size = 2.5, color = "blue", fontface = "bold")
      }

      # Major (tree) edge: tree-parent -> hybrid (red).
      tp <- net$edge[net$edge[, 2] == to, 1]
      if (length(tp) && !is.na(g$tree_gamma)) {
        xm2 <- (d$x[d$node == tp[1]] + d$x[d$node == to]) / 2
        ym2 <- (d$y[d$node == tp[1]] + d$y[d$node == to]) / 2
        p <- p + ggplot2::annotate("text", x = xm2, y = ym2 + 0.3,
          label = gamma_label(g$tree_gamma),
          size = 2.5, color = "red", fontface = "bold")
      }
    }
    xmax <- max(d$x, na.rm = TRUE)
    p <- p + ggplot2::xlim(NA, xmax * 1.4)
  }

  p
}

# Expand straight (slanted) tree segments into right-angle elbows. Internal.
to_rectangular <- function(tree) {
  vert <- data.frame(x = tree$x, y = tree$y, xend = tree$x, yend = tree$yend,
                     .net = tree$.net)
  horiz <- data.frame(x = tree$x, y = tree$yend, xend = tree$xend, yend = tree$yend,
                      .net = tree$.net)
  rbind(vert, horiz)
}

#' DensiTree-style consensus/discrepancy overlay of a set of networks
#'
#' Overlays every network in a set on one shared tip order with alpha
#' transparency: the backbone trees form a "consensus cloud" (agreement is dense,
#' conflict is diffuse) and the reticulation edges form a distinct cloud. On top,
#' a solid consensus ("root-canal") backbone is drawn -- colored by clade support
#' -- together with the frequent reticulation events (donor -> recipient arrows),
#' so the dominant signal is unambiguous regardless of `alpha`. This is anansi's
#' primary view (DensiTree adapted to networks; see docs/DESIGN.md section 7 and
#' D4/D5). Returns a composable `ggplot` object.
#'
#' @param x An [anansi_netset] (see [read_networks_csv]).
#' @param tip_order Optional explicit tip order; otherwise computed via
#'   [consensus_tip_order] using `method`.
#' @param method Tip-ordering method (see [consensus_tip_order]); default `"mode"`.
#' @param outgroup Optional taxa to pin to one end of the figure (see
#'   [consensus_tip_order]).
#' @param outgroup_position `"top"` (default) or `"bottom"` placement for
#'   `outgroup`.
#' @param snap_to_consensus If TRUE, snap the shared tip order to the consensus
#'   topology so the consensus backbone has no crossings (see
#'   [consensus_tip_order]). Default FALSE.
#' @param mode `"cladogram"` (default) or `"phylogram"`.
#' @param layout `"slanted"` (default) or `"rectangular"` for the cloud backbone.
#' @param keep Optional character vector of taxa; if given, all networks are
#'   restricted to this subset via [restrict_taxa] before plotting.
#' @param top_n,top_by If `top_n` is set, restrict to the top-N networks by the
#'   metadata column `top_by` (see [top_networks]) before plotting.
#' @param tree_color,ret_color Colors for the backbone and reticulation clouds.
#' @param alpha,ret_alpha Per-edge cloud opacities; default to roughly `1/N` and
#'   `2/N` so consensus accumulates while conflict stays faint.
#' @param linewidth Cloud edge width.
#' @param jitter Max per-network x-offset to separate overlapping edges
#'   (DensiTree-style); 0 (default) disables it.
#' @param ret_linetype Line type for reticulation edges (default `"dashed"`).
#' @param consensus If TRUE (default), draw the consensus backbone + reticulations.
#' @param consensus_p Majority threshold for the consensus tree (default 0.5).
#' @param color_by_support If TRUE (default), color the consensus backbone by
#'   clade support (grey -> `consensus_color`).
#' @param consensus_color,consensus_linewidth Color/width of the consensus backbone.
#' @param consensus_ret If TRUE (default), draw consensus reticulation arrows.
#' @param consensus_ret_min Minimum event frequency to draw a consensus
#'   reticulation arrow (default 0.1).
#' @param consensus_ret_color Color of consensus reticulation arrows/edges.
#' @param ret_edge_frac Fraction (0-1) of the way from each reticulation
#'   endpoint's node toward its parent at which to anchor the edge, so
#'   reticulations attach along a branch rather than at the tip/MRCA node.
#'   Default 0.1; 0 reproduces at-the-node anchoring.
#' @param reticulation_style `"arrow"` (default) draws each consensus
#'   reticulation as a directed donor -> recipient arrow; `"hybrid"` keys events
#'   direction-agnostically and draws two undirected dotted edges into the hybrid
#'   node (so direction-unstable reticulations display identically).
#' @param tip_labels Whether to draw tip labels.
#' @param tip_size,tip_offset Tip-label text size and x-offset past the tips.
#' @param consistent_only If TRUE (default), drop networks with divergent taxa.
#' @param title Plot title; a default is generated if NULL.
#' @return A `ggplot` object.
#' @seealso [layout_netset], [consensus_network], [consensus_tip_order],
#'   [top_networks]
#' @examples
#' a <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
#' b <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
#' p <- densinet(anansi_netset(list(a, b)), method = "first")
#' @export
densinet <- function(x, tip_order = NULL, method = "mode",
                     outgroup = NULL, outgroup_position = c("top", "bottom"),
                     snap_to_consensus = FALSE,
                     mode = c("cladogram", "phylogram"),
                     layout = c("slanted", "rectangular"),
                     keep = NULL, top_n = NULL, top_by = NULL,
                     tree_color = "steelblue", ret_color = "firebrick",
                     alpha = NULL, ret_alpha = NULL, linewidth = 0.3,
                     jitter = 0, ret_linetype = "dashed",
                     consensus = TRUE, consensus_p = 0.5,
                     color_by_support = TRUE, consensus_color = "black",
                     consensus_linewidth = 0.7,
                     consensus_ret = TRUE, consensus_ret_min = 0.1,
                     consensus_ret_color = "darkred", ret_edge_frac = 0.1,
                     reticulation_style = c("arrow", "hybrid"),
                     tip_labels = TRUE, tip_size = 3, tip_offset = 0.02,
                     consistent_only = TRUE, title = NULL) {
  if (!inherits(x, "anansi_netset")) {
    rlang::abort("densinet() needs an anansi_netset (see read_networks_csv()).")
  }
  mode <- match.arg(mode)
  layout <- match.arg(layout)
  outgroup_position <- match.arg(outgroup_position)
  reticulation_style <- match.arg(reticulation_style)
  if (!is.null(keep)) x <- restrict_taxa(x, keep)
  if (!is.null(top_n)) x <- top_networks(x, top_n, by = top_by)

  seg <- layout_netset(x, tip_order = tip_order, method = method, mode = mode,
                       scale_x = TRUE, jitter = jitter,
                       consistent_only = consistent_only, outgroup = outgroup,
                       outgroup_position = outgroup_position,
                       snap_to_consensus = snap_to_consensus,
                       consensus_p = consensus_p)
  ord <- attr(seg, "tip_order")
  N <- length(unique(seg$.net))
  if (is.null(alpha)) alpha <- max(1 / N, 0.005)
  if (is.null(ret_alpha)) ret_alpha <- max(2 / N, 0.03)

  tree <- seg[seg$kind == "tree", ]
  ret  <- seg[seg$kind == "reticulation", ]
  if (layout == "rectangular") tree <- to_rectangular(tree)

  # --- cloud ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = tree,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = .net),
      color = tree_color, alpha = alpha, linewidth = linewidth)
  if (nrow(ret)) {
    p <- p + ggplot2::geom_segment(
      data = ret,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = .net),
      color = ret_color, alpha = ret_alpha, linewidth = linewidth,
      linetype = ret_linetype)
  }

  # --- consensus overlay (drawn slanted on top, regardless of cloud layout) ---
  if (consensus) {
    cons <- consensus_network(x, p = consensus_p,
                              directed = reticulation_style == "arrow")
    cs <- consensus_segments(cons, tip_order = ord, mode = mode,
                             ret_min = consensus_ret_min,
                             ret_edge_frac = ret_edge_frac)
    if (color_by_support) {
      p <- p +
        ggplot2::geom_segment(
          data = cs$backbone,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = support),
          linewidth = consensus_linewidth) +
        ggplot2::scale_color_gradient(low = "grey80", high = consensus_color,
                                      limits = c(0, 1), name = "clade\nsupport")
    } else {
      p <- p + ggplot2::geom_segment(
        data = cs$backbone,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = consensus_color, linewidth = consensus_linewidth)
    }
    if (consensus_ret && !is.null(cs$reticulations) && nrow(cs$reticulations)) {
      if (reticulation_style == "hybrid") {
        # Undirected: two dotted edges into the hybrid node, no arrowhead, so
        # direction-unstable reticulations read (and aggregate) identically.
        p <- p + ggplot2::geom_segment(
          data = cs$reticulations,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          color = consensus_ret_color, linewidth = consensus_linewidth,
          linetype = ret_linetype)
      } else {
        p <- p + ggplot2::geom_segment(
          data = cs$reticulations,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          color = consensus_ret_color, linewidth = consensus_linewidth,
          linetype = ret_linetype,
          arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed"))
      }
    }
  }

  # --- tips + theme ---
  if (tip_labels) {
    tlab <- data.frame(label = ord, y = seq_along(ord))
    p <- p +
      ggplot2::geom_text(data = tlab,
                         ggplot2::aes(x = 1 + tip_offset, y = y, label = label),
                         hjust = 0, size = tip_size) +
      # NA lower bound so a negative `jitter` offset is not clipped; the upper
      # bound reserves room for tip labels.
      ggplot2::scale_x_continuous(limits = c(NA, 1 + tip_offset + 0.5))
  }

  ttl <- title %||% sprintf("anansi densinet: %d networks%s", N,
    if (!is.null(top_n)) sprintf(" (top %d)", top_n) else "")
  p +
    ggplot2::ggtitle(ttl) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA))
}
