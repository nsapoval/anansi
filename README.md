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

### Web app

`run_anansi_app()` launches an interactive Shiny interface: upload a CSV (or
paste/upload raw extended-Newick, or load the bundled example), tune the key
`densinet()` controls — with the rest behind an **Advanced options** toggle —
preview the figure, and download it as PNG/PDF (SVG when `svglite` is installed).

```r
# install the (suggested) app dependencies once:
install.packages(c("shiny", "bslib", "colourpicker", "svglite"))  # only `shiny` is required
run_anansi_app()
```

The same input → figure mapping is available programmatically via
`build_densinet(netset, params)` and `netset_from_enewick(text)`.

#### Live app & deployment

A hosted copy runs on shinyapps.io: **https://&lt;account&gt;.shinyapps.io/anansi/**
(replace `<account>` with the shinyapps.io account name).

Deployment is automated: every push to `main` that touches `R/`, `inst/`,
`DESCRIPTION`, `NAMESPACE`, or `deploy/` triggers the
[`Deploy to shinyapps.io`](.github/workflows/deploy-shinyapps.yaml) GitHub Action,
which installs `anansi` from GitHub plus its dependency closure and runs
[`deploy/deploy.R`](deploy/deploy.R). Because the app is a thin UI over the
package, the public repo must stay **public** so the build server can install
`anansi` via `install_github`.

One-time setup for the auto-deploy:

1. Create a token at shinyapps.io → **Account → Tokens → Add Token → Show secret**.
2. Add three GitHub repository secrets (repo **Settings → Secrets and variables →
   Actions**): `SHINYAPPS_NAME` (account name), `SHINYAPPS_TOKEN`, `SHINYAPPS_SECRET`.

To deploy manually from a local clone instead:

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
Sys.setenv(SHINYAPPS_NAME = "<account>")
source("deploy/deploy.R")   # deploys inst/shiny as app "anansi"
```

Working now: `read_networks_csv()`, `parse_network()`, `consensus_tip_order()`,
`layout_network()`/`layout_netset()`, `densinet()` (with `consensus`/`top_n`/`keep`),
`plot_network()`, `clade_frequencies()`, `reticulation_frequencies()`,
`consensus_network()`, `restrict_taxa()`, `top_networks()`,
`run_anansi_app()`, `build_densinet()`, `netset_from_enewick()`.

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
