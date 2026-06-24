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
  if (!is.null(tree$.value)) {        # carry the gradient value through the split
    vert$.value <- tree$.value
    horiz$.value <- tree$.value
  }
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
#'   `2/N` so consensus accumulates while conflict stays faint. When `color_by`
#'   is set, `alpha` instead defaults to `max(3/N, 0.25)` so the gradient reads.
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
#' @param color_by Optional name of a numeric metadata column (see
#'   [read_networks_csv]); when set, each tree in the cloud is colored by its
#'   value along a gradient. To keep a single color scale, the scaffold is then
#'   drawn solid (consensus `color_by_support` is auto-disabled).
#' @param color_palette Gradient for `color_by`: one of the viridis options
#'   `"viridis"` (default), `"magma"`, `"plasma"`, `"inferno"`, `"cividis"`, or
#'   `"gradient"` for a two-color `color_low` -> `color_high` ramp.
#' @param color_direction Viridis direction (`1` or `-1`); ignored for
#'   `"gradient"`.
#' @param color_low,color_high Endpoints of the `"gradient"` palette.
#' @param color_legend Legend title for `color_by` (defaults to the column name).
#' @param color_trim Strong-outlier handling for `color_by`. A non-negative IQR
#'   multiplier `k`: per-network values outside the Tukey fences
#'   `[Q1 - k*IQR, Q3 + k*IQR]` are treated as outliers (`3` ~ "strong", `1.5` ~
#'   "mild"). `NULL` (default) disables trimming.
#' @param color_trim_action What to do with `color_trim` outliers: `"clamp"`
#'   (default) keeps every tree but pins the color scale to the robust range so
#'   outliers saturate at the end color; `"drop"` removes the outlier networks
#'   from the overlay.
#' @param backbone Optional explicit scaffold drawn solid in place of the
#'   consensus backbone: a Newick / extended-Newick string, [anansi_network],
#'   `evonet`, or `phylo` over (a subset of) the figure's taxa. If it is a
#'   network, its reticulation edges are drawn too.
#' @param backbone_color,backbone_ret_color,backbone_linewidth Colors of the
#'   explicit `backbone`'s tree and reticulation edges, and their width (width
#'   defaults to `consensus_linewidth`).
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
                     color_by = NULL, color_palette = "viridis",
                     color_direction = 1, color_low = "grey80",
                     color_high = "firebrick", color_legend = NULL,
                     color_trim = NULL, color_trim_action = c("clamp", "drop"),
                     backbone = NULL, backbone_color = "black",
                     backbone_ret_color = "darkred", backbone_linewidth = NULL,
                     consistent_only = TRUE, title = NULL) {
  if (!inherits(x, "anansi_netset")) {
    rlang::abort("densinet() needs an anansi_netset (see read_networks_csv()).")
  }
  mode <- match.arg(mode)
  layout <- match.arg(layout)
  outgroup_position <- match.arg(outgroup_position)
  reticulation_style <- match.arg(reticulation_style)
  color_palette <- match.arg(color_palette,
    c("viridis", "magma", "plasma", "inferno", "cividis", "gradient"))
  if (!is.null(keep)) x <- restrict_taxa(x, keep)
  if (!is.null(top_n)) x <- top_networks(x, top_n, by = top_by)

  seg <- layout_netset(x, tip_order = tip_order, method = method, mode = mode,
                       scale_x = TRUE, jitter = jitter,
                       consistent_only = consistent_only, outgroup = outgroup,
                       outgroup_position = outgroup_position,
                       snap_to_consensus = snap_to_consensus,
                       consensus_p = consensus_p, color_by = color_by)
  ord <- attr(seg, "tip_order")
  use_value <- !is.null(color_by)
  color_trim_action <- match.arg(color_trim_action)

  # Strong-outlier handling for the value coloring: a few extreme trees can
  # otherwise blow out the gradient so the bulk reads as one flat color. Tukey
  # fences on the per-network values (Q1 - k*IQR, Q3 + k*IQR; k = color_trim,
  # 3 ~ "strong", 1.5 ~ "mild") define the robust range. "clamp" keeps every
  # tree but pins the scale to that range (outliers saturate at the end color);
  # "drop" removes the outlier networks from the overlay.
  color_limits <- NULL
  if (use_value && !is.null(color_trim) && is.finite(color_trim) && color_trim >= 0) {
    pv <- tapply(seg$.value, seg$.net, function(v) v[1L])   # one value per network
    qs <- stats::quantile(pv, c(0.25, 0.75), names = FALSE, na.rm = TRUE)
    iqr <- qs[2L] - qs[1L]
    if (iqr > 0) {
      lo <- qs[1L] - color_trim * iqr
      hi <- qs[2L] + color_trim * iqr
      if (color_trim_action == "drop") {
        drop_nets <- as.integer(names(pv))[!is.na(pv) & (pv < lo | pv > hi)]
        if (length(drop_nets)) {
          seg <- seg[!seg$.net %in% drop_nets, , drop = FALSE]
          rlang::inform(sprintf(
            "color_trim: dropped %d outlier network(s) outside [%.3g, %.3g] (IQR x %g).",
            length(drop_nets), lo, hi, color_trim))
        }
      } else {
        lo <- max(lo, min(pv, na.rm = TRUE))   # don't pad the legend past the data
        hi <- min(hi, max(pv, na.rm = TRUE))
        if (hi > lo) {
          seg$.value <- pmin(pmax(seg$.value, lo), hi)
          color_limits <- c(lo, hi)
        }
      }
    }
  }

  N <- length(unique(seg$.net))
  # When coloring by value the gradient must read, so the cloud uses a higher
  # opacity than the consensus/discrepancy default (which favors faintness). The
  # 3/N term keeps small sets vivid; the 0.05 floor avoids saturating large ones.
  if (is.null(alpha)) alpha <- if (use_value) max(3 / N, 0.05) else max(1 / N, 0.005)
  if (is.null(ret_alpha)) ret_alpha <- max(2 / N, 0.03)

  tree <- seg[seg$kind == "tree", ]
  ret  <- seg[seg$kind == "reticulation", ]
  if (layout == "rectangular") tree <- to_rectangular(tree)

  # --- cloud ---
  p <- ggplot2::ggplot()
  if (use_value) {
    # Each tree shaded by its value; overlapping edges blend through alpha.
    p <- p + ggplot2::geom_segment(
      data = tree,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   color = .value, group = .net),
      alpha = alpha, linewidth = linewidth)
    legend_name <- color_legend %||% color_by
    p <- p + if (color_palette == "gradient") {
      ggplot2::scale_color_gradient(low = color_low, high = color_high,
                                    name = legend_name, limits = color_limits)
    } else {
      ggplot2::scale_color_viridis_c(option = color_palette,
                                     direction = color_direction,
                                     name = legend_name, limits = color_limits)
    }
  } else {
    p <- p + ggplot2::geom_segment(
      data = tree,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = .net),
      color = tree_color, alpha = alpha, linewidth = linewidth)
  }
  if (nrow(ret)) {
    p <- p + ggplot2::geom_segment(
      data = ret,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = .net),
      color = ret_color, alpha = ret_alpha, linewidth = linewidth,
      linetype = ret_linetype)
  }

  # --- scaffold (drawn slanted on top, regardless of cloud layout) ---
  # An explicit backbone replaces the consensus scaffold when supplied.
  if (!is.null(backbone)) {
    bb_net <- if (is.character(backbone)) parse_network(backbone) else backbone
    bb_ev <- ensure_evonet(bb_net)
    extra <- setdiff(bb_ev$tip.label, ord)
    if (length(extra)) {
      rlang::abort(sprintf("backbone has taxa not in the figure: %s",
                           paste(extra, collapse = ", ")))
    }
    bb <- layout_network(bb_net, tip_order = ord, mode = mode)
    s <- scale_x01(bb$x, bb$xend)        # scale tree + reticulation edges together
    bb$x <- s$x; bb$xend <- s$xend
    bb_lwd <- backbone_linewidth %||% consensus_linewidth
    p <- p + ggplot2::geom_segment(
      data = bb[bb$kind == "tree", , drop = FALSE],
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = backbone_color, linewidth = bb_lwd)
    bb_ret <- bb[bb$kind == "reticulation", , drop = FALSE]
    if (nrow(bb_ret)) {       # an explicit network backbone keeps its hybrid edges
      bb_arrow <- if (reticulation_style == "arrow")
        grid::arrow(length = grid::unit(0.1, "inches"), type = "closed") else NULL
      p <- p + ggplot2::geom_segment(
        data = bb_ret,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = backbone_ret_color, linewidth = bb_lwd,
        linetype = ret_linetype, arrow = bb_arrow)
    }
  } else if (consensus) {
    cons <- consensus_network(x, p = consensus_p,
                              directed = reticulation_style == "arrow")
    cs <- consensus_segments(cons, tip_order = ord, mode = mode,
                             ret_min = consensus_ret_min,
                             ret_edge_frac = ret_edge_frac)
    # The support gradient and the cloud value gradient would be two continuous
    # color scales on one plot; drop support coloring so the cloud scale wins.
    if (color_by_support && use_value) {
      rlang::inform(paste0(
        "color_by is set, so the consensus backbone is drawn solid ",
        "(color_by_support disabled to keep a single color scale)."),
        .frequency = "once",
        .frequency_id = "anansi_color_by_support_off")
    }
    if (color_by_support && !use_value) {
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
