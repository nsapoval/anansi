# anansi

**Consensus & discrepancy visualization for sets of phylogenetic networks.**

anansi takes a *set* of explicit reticulate phylogenetic networks (extended Newick /
`ape::evonet`) over a fixed taxon set and overlays them DensiTree-style: networks are
drawn on a **shared tip ordering** with alpha transparency, so **agreement appears
dense** and **conflict appears diffuse**. Both backbone topology and reticulation
(hybrid-edge) placement are shown together. You can also restrict every network to a
**subset of taxa** and re-view.

It builds on [`ape`](https://cran.r-project.org/package=ape),
[`tanggle`](https://bioconductor.org/packages/tanggle) and
[`ggtree`](https://bioconductor.org/packages/ggtree); output is a composable `ggplot`
object. The name nods to *Anansi* the spider — fitting for webs of reticulate edges.

> **Status: complete (phases 0-6 + W1).** Reading network sets (via a native
> extended-Newick parser that handles nested/stacked reticulations), the shared
> layout engine, the `densinet()` consensus/discrepancy overlay, consensus/support
> encoding, reticulation-aware taxon subsetting, and phylogram/jitter polish all
> work; all sample files parse 100%. See [`docs/WORKPLAN.md`](docs/WORKPLAN.md) for
> the feature tracker.

## Why

A single network plot can't tell you how stable an inference is. Inference tools emit
*ranked sets* of candidate networks; anansi lets you see, across the whole set, where
the networks agree and where they don't — in topology **and** in where reticulations
are placed.

## Install

Requires R ≥ 4.1. Two dependencies (`tanggle`, `ggtree`) are on Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install(c("ggtree", "tanggle"))      # Bioconductor deps
install.packages(c("ape", "ggplot2", "tidytree", "dplyr", "rlang"))

# development install of anansi (from a local clone):
# install.packages("devtools"); devtools::install(".")
devtools::load_all(".")                            # for development
```

## Usage

```r
library(anansi)   # or devtools::load_all(".")

# 1. Read a set of networks from a CSV column of extended-Newick strings
nets <- read_networks_csv("data/Lacerta_agilis_inferred_networks.csv")
nets                                                # summary: networks, taxa, taxa_ok

# 2. Consensus/discrepancy overlay (cladogram; tip order via MDS), a ggplot object
densinet(nets, method = "mds")                      # backbone cloud + reticulation cloud

# 3. Consensus overlay focused on the best networks, restricted to a taxon subset
densinet(nets, top_n = 50,
         keep = c("Lacerta_agilis", "Zootoca_vivipara", "Podarcis_muralis",
                  "Podarcis_siculus", "Podarcis_liolepis"))

# 4. Plot a single network (gamma inheritance labels)
plot_network(nets$networks[[1]])

# Consensus summaries
clade_frequencies(nets)
reticulation_frequencies(nets)
consensus_network(nets, p = 0.5)
```

Working now: `read_networks_csv()`, `parse_network()`, `consensus_tip_order()`,
`layout_network()`/`layout_netset()`, `densinet()` (with `consensus`/`top_n`/`keep`),
`plot_network()`, `clade_frequencies()`, `reticulation_frequencies()`,
`consensus_network()`, `restrict_taxa()`, `top_networks()`.

Layout extras: `densinet(mode = "phylogram")`, `densinet(jitter = 0.01)`,
`densinet(layout = "rectangular")`.

Deferred: an R-Markdown vignette (pending pandoc; see [`docs/TUTORIAL.md`](docs/TUTORIAL.md))
and DensiTree-style top-N topology coloring — see the workplan.

## Documentation

- [`docs/TUTORIAL.md`](docs/TUTORIAL.md) — hands-on walkthrough on the bundled
  example dataset (`anansi_example()`).
- [`docs/DESIGN.md`](docs/DESIGN.md) — design decisions & rationale (the DensiTree →
  network adaptation, coordinate system, consensus definitions, open questions).
- [`docs/WORKPLAN.md`](docs/WORKPLAN.md) — phased roadmap + feature-status tracker.
- [`docs/DATA.md`](docs/DATA.md) — input CSV/eNewick formats & the taxon set.
- [`docs/REFERENCES.md`](docs/REFERENCES.md) — papers & R packages we build on.

## Sample data

`data/` holds four sets of inferred networks (60–460 each) over the same 15 lizard
taxa — see [`docs/DATA.md`](docs/DATA.md).

## License

AGPL-3 (see `LICENSE`).
