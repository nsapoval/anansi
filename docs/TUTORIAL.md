# anansi — Tutorial

A practical walkthrough of the consensus/discrepancy workflow. Every code block
runs against the bundled example dataset, so you can copy-paste it as-is.

> A formal R-Markdown vignette (`vignette("anansi")`) is deferred until a pandoc
> toolchain is available in the build environment; this tutorial covers the same
> ground and is kept in sync with the code.

## 1. Load the package and a dataset

```r
library(anansi)            # during development: devtools::load_all(".")

# A bundled example: 40 inferred networks over 15 lizard taxa.
nets <- read_networks_csv(anansi_example())
nets
#> <anansi_netset>
#>   networks: 40
#>   taxa:     15 (modal/reference set)
#>   reticulations/network: 1:32  2:8
#>   taxa_ok: 40/40 networks match the reference taxon set
#>   meta columns: index, num_reticulations, network_rank, log_pseudo_likelihood, ...
```

Your own data is a CSV with one column of extended-Newick strings (named
`network`, `newick`, or `enewick` — auto-detected) and any number of metadata
columns:

```r
nets <- read_networks_csv("path/to/your_networks.csv")
```

Networks whose taxon set diverges from the majority (e.g. `ape::read.evonet`
artifacts on nested reticulations) are kept but flagged; keep the clean ones with
`nets[nets$taxa_ok]`. See [DATA.md](DATA.md).

## 2. The consensus/discrepancy overlay

`densinet()` is the main view: all networks overlaid on a shared tip order, with a
solid consensus backbone (colored by clade support) and the frequent reticulation
events drawn on top. It returns a `ggplot`, so you can add layers/themes.

```r
densinet(nets)                                  # cladogram, tip order via "mode"
densinet(nets, method = "mds")                  # MDS tip ordering
densinet(nets, method = "closest_leaf")         # DensiTree's greedy ordering
```

Tune the look:

```r
densinet(nets,
         alpha = 0.05,                # cloud opacity (default ~1/N)
         consensus_ret_min = 0.2,     # only draw reticulations seen in >=20%
         tree_color = "grey50",
         consensus_color = "navy")
```

Turn the consensus overlay off to see the raw cloud, or hide tip labels:

```r
densinet(nets, consensus = FALSE)
densinet(nets, tip_labels = FALSE)
```

## 3. Focus: top networks and taxon subsets

Restrict to the best-scoring networks (by a metadata column) and/or to a subset of
taxa (reticulation-aware pruning):

```r
# Top 20 networks by pseudo-likelihood
densinet(nets, top_n = 20)

# Restrict every network to 6 taxa, then overlay
keep <- c("Lacerta_agilis", "Zootoca_vivipara", "Podarcis_muralis",
          "Podarcis_siculus", "Podarcis_liolepis", "Podarcis_bocagei")
densinet(nets, keep = keep)

# Combine both
densinet(nets, top_n = 20, keep = keep)
```

The same operations are available as standalone transforms:

```r
sub  <- restrict_taxa(nets, keep)               # an anansi_netset on `keep`
best <- top_networks(nets, 20)                  # an anansi_netset of 20 networks
```

## 4. Layout options

```r
densinet(nets, mode = "phylogram")              # branch-length scaled (vs cladogram)
densinet(nets, layout = "rectangular")          # right-angle backbone edges
densinet(nets, jitter = 0.01)                   # nudge overlapping edges apart
```

## 5. Consensus summaries (the numbers behind the picture)

```r
clade_frequencies(nets)          # how often each clade appears
reticulation_frequencies(nets)   # reticulation events (donor->recipient) + mean gamma
consensus_network(nets, p = 0.5) # the majority-rule "root-canal" network
```

`reticulation_frequencies()` keys each reticulation by the pair of clades it
connects, so equivalent reticulations across networks are tallied together —
the basis for the consensus arrows in `densinet()`.

## 6. A single network

```r
plot_network(nets$networks[[1]])   # one network, with gamma inheritance labels
```

## 7. Saving figures

`densinet()` returns a `ggplot`; save it with a white background:

```r
p <- densinet(nets, method = "mds")
ggplot2::ggsave("consensus.png", p, width = 10, height = 6.5, dpi = 150, bg = "white")
```

See `dev/make_figures.R` for the full set of example figures.
