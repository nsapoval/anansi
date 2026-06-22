library(ape)
library(tanggle)
library(ggtree)

# Read the results CSV
df <- read.csv("data/d10k_x100_results.csv", stringsAsFactors = FALSE)

# Parse the network as an evonet object
net_id <- 1
nwk <- df$network[net_id]
net <- read.evonet(text = nwk)

# --- Extract inheritance probabilities from extended Newick ---
# read.evonet does not populate net$prob, so parse gammas from the string.
# Each #Hx appears twice: 1st occurrence = tree (major) edge, 2nd = reticulation (minor) edge.
h_matches <- gregexpr("#H\\d+:[^,()]+::[0-9.eE+-]+", nwk, perl = TRUE)
h_strings <- regmatches(nwk, h_matches)[[1]]
h_labels <- sub("^(#H\\d+):.*", "\\1", h_strings)
h_gammas <- as.numeric(sub("^#H\\d+:[^:]+::", "", h_strings))

# For each hybrid label, split into tree-edge gamma (1st) and reticulation gamma (2nd)
gamma_table <- do.call(rbind, lapply(unique(h_labels), function(lbl) {
  g <- h_gammas[h_labels == lbl]
  data.frame(label = lbl, tree_gamma = g[1], ret_gamma = g[2],
             stringsAsFactors = FALSE)
}))

# Map hybrid labels to internal node numbers
ntip <- length(net$tip.label)
nl <- net$node.label
hybrid_node_map <- setNames(ntip + which(nl != ""), nl[nl != ""])

# --- Plot as cladogram (topology only) ---
# Remove edge lengths so ggevonet computes uniform node depths via node_depth_evonet
net$edge.length <- NULL
p <- ggevonet(net, layout = "slanted") +
  geom_tiplab(size = 3) +
  ggtitle(paste0("Dataset: ", df$dataset[net_id],
                 "  |  Reticulations: ", df$reticulations[net_id],
                 "  |  Log-prob: ", round(df$log_probability[net_id], 2)))

d <- p$data
ret <- net$reticulation

# Annotate each pair of edges to hybrid nodes with their gamma values
for (i in seq_len(nrow(ret))) {
  hybrid_node <- ret[i, 2]
  lbl <- names(hybrid_node_map)[hybrid_node_map == hybrid_node]
  g <- gamma_table[gamma_table$label == lbl, ]

  # Reticulation (minor) edge
  from <- ret[i, 1]
  to   <- ret[i, 2]
  x_mid <- (d$x[d$node == from] + d$x[d$node == to]) / 2
  y_mid <- (d$y[d$node == from] + d$y[d$node == to]) / 2
  p <- p + annotate("text", x = x_mid, y = y_mid + 0.3,
                     label = sprintf("\u03b3 = %.2f", g$ret_gamma),
                     size = 2.5, color = "blue", fontface = "bold")

  # Tree (major) edge to same hybrid node
  tree_parent <- net$edge[net$edge[, 2] == to, 1]
  if (length(tree_parent) > 0) {
    x_mid2 <- (d$x[d$node == tree_parent[1]] + d$x[d$node == to]) / 2
    y_mid2 <- (d$y[d$node == tree_parent[1]] + d$y[d$node == to]) / 2
    p <- p + annotate("text", x = x_mid2, y = y_mid2 + 0.3,
                       label = sprintf("\u03b3 = %.2f", g$tree_gamma),
                       size = 2.5, color = "red", fontface = "bold")
  }
}

x_max <- max(d$x, na.rm = TRUE)
p <- p + xlim(NA, x_max * 1.4)

png(paste0("network_", as.character(net_id), ".png"), width = 1200, height = 900, res = 150)
print(p)
dev.off()

cat(paste0("network_", as.character(net_id), ".png\n"))
