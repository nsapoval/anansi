# dev/make_figures.R --------------------------------------------------------
# Render a spread of densinet() outputs for visual inspection.
# Run from the package root:  Rscript dev/make_figures.R
# Outputs go to figures/ (git/Rbuild-ignored).

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(ggplot2))

dir.create("figures", showWarnings = FALSE)
save <- function(p, name, w = 10, h = 7, dpi = 150) {
  ggsave(file.path("figures", name), p, width = w, height = h, dpi = dpi,
         bg = "white")
  cat("  figures/", name, "\n", sep = "")
}

lac <- suppressWarnings(read_networks_csv("data/Lacerta_agilis_inferred_networks.csv"))
zoo <- suppressWarnings(read_networks_csv("data/Zootoca_vivipara_inferred_networks.csv"))
d10 <- suppressWarnings(read_networks_csv("data/d10k_mostlen_inferred_networks.csv"))

cat("Rendering figures...\n")

# 1. Tip-ordering comparison on Lacerta (which heuristic reads cleanest?)
for (m in c("first", "mode", "mds", "closest_leaf")) {
  save(densinet(lac, method = m, tip_size = 3.2,
                title = sprintf("Lacerta_agilis (230) - tip order: %s", m)),
       sprintf("01_lacerta_order_%s.png", m))
}

# 2. Alpha sweep on Lacerta (mds): contrast of consensus vs conflict
for (a in c(0.005, 0.02, 0.06)) {
  save(densinet(lac, method = "mds", alpha = a, ret_alpha = max(a * 3, 0.05),
                tip_size = 3.2,
                title = sprintf("Lacerta_agilis (mds) - alpha = %.3f", a)),
       sprintf("02_lacerta_alpha_%03.0f.png", a * 1000))
}

# 3. slanted vs rectangular backbone
save(densinet(lac, method = "mds", layout = "slanted", tip_size = 3.2,
              title = "Lacerta_agilis - slanted"), "03_lacerta_slanted.png")
save(densinet(lac, method = "mds", layout = "rectangular", tip_size = 3.2,
              title = "Lacerta_agilis - rectangular"), "03_lacerta_rectangular.png")

# 4. Dataset comparison (closest_leaf order)
save(densinet(zoo, method = "closest_leaf", tip_size = 3.2,
              title = "Zootoca_vivipara (230)"), "04_zootoca_closest_leaf.png")
save(densinet(d10, method = "closest_leaf", tip_size = 3.2,
              title = "d10k_mostlen (consistent subset, 1-3 reticulations)"),
     "04_d10k_closest_leaf.png")

# 5. Top-N by likelihood vs all (does focusing on the best networks read cleaner?)
ord <- order(lac$meta$log_pseudo_likelihood, decreasing = TRUE)
lac_top <- lac[ord[seq_len(min(50, length(ord)))]]
save(densinet(lac_top, method = "mds", tip_size = 3.2,
              title = "Lacerta_agilis - top 50 by pseudo-likelihood (mds)"),
     "05_lacerta_top50_mds.png")

# 6. Single-network reference (gamma labels)
save(plot_network(lac$networks[[1]]), "06_single_network.png", w = 9, h = 6)

# 7. Consensus overlay (Phase 4): support-colored backbone + consensus reticulations
save(densinet(lac, method = "mds",
              title = "Lacerta_agilis (230): cloud + consensus overlay"),
     "07_lacerta_consensus.png")
save(densinet(lac, method = "mds", top_n = 50,
              title = "Lacerta_agilis: top 50 + consensus overlay"),
     "07_lacerta_top50_consensus.png")
save(densinet(d10, method = "closest_leaf",
              title = "d10k_mostlen: cloud + consensus (reticulation events >=10%)"),
     "07_d10k_consensus.png")

# 8. Taxon subsetting (Phase 5): restrict all networks to a subset, then overlay
keep9 <- c("Lacerta_agilis", "Zootoca_vivipara", "Podarcis_muralis",
           "Podarcis_siculus", "Podarcis_liolepis", "Podarcis_bocagei",
           "Podarcis_vaucheri", "Podarcis_cretensis", "Podarcis_erhardii")
save(densinet(lac, method = "mds", keep = keep9,
              title = "Lacerta_agilis restricted to 9 taxa (consensus overlay)"),
     "08_lacerta_restricted.png")

# 9. Phylogram mode and jitter (Phase 6)
save(densinet(lac, method = "mds", mode = "phylogram", top_n = 80,
              title = "Lacerta_agilis (top 80): phylogram (branch-length scaled)"),
     "09_lacerta_phylogram.png")
save(densinet(lac, method = "mds", jitter = 0.015, top_n = 80,
              title = "Lacerta_agilis (top 80): cladogram with jitter"),
     "09_lacerta_jitter.png")

cat("Done.\n")
