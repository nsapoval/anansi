# anansi — References

Curated background for the design. Grouped by topic.

## DensiTree (the inspiration)

- Bouckaert, R. (2010). *DensiTree: making sense of sets of phylogenetic trees.*
  Bioinformatics 26(10):1372–1373.
  https://academic.oup.com/bioinformatics/article/26/10/1372/192963
- Bouckaert, R. & Heled, J. (2014). *DensiTree 2: Seeing Trees Through the Forest.*
  bioRxiv 012401. https://www.biorxiv.org/content/10.1101/012401v1.full
- DensiTree site & manual (Remco Bouckaert):
  https://www.cs.auckland.ac.nz/~remco/DensiTree/ ·
  https://www.cs.auckland.ac.nz/~remco/DensiTree/DensiTreeManual.v2.2.pdf
- Taxon ordering for tree visualization: *Taxon ordering in phylogenetic trees: a
  workbench test.* BMC Bioinformatics 12:58 (2011).
  https://link.springer.com/article/10.1186/1471-2105-12-58

## DensiTree-style overlay in R

- ggtree `ggdensitree()` (overlay of multiple **trees**; tip-order strategies
  "mode"/"mds"): https://rdrr.io/bioc/ggtree/man/ggdensitree.html ·
  source: https://github.com/YuLab-SMU/ggtree/blob/devel/R/ggdensitree.R
- phangorn `densiTree()`: https://rdrr.io/cran/phangorn/man/densiTree.html

## Phylogenetic networks: representation & file formats

- `ape::evonet` (explicit reticulate networks; `$reticulation` matrix;
  `read.evonet`/`write.evonet`): https://search.r-project.org/CRAN/refmans/ape/html/evonet.html
- Cardona, Rosselló, Valiente (2008). *Extended Newick: it is time for a standard
  representation of phylogenetic networks.* BMC Bioinformatics.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2621367/
- SiPhyNetwork (evonet + inheritance probabilities; eNewick I/O):
  https://jjustison.github.io/SiPhyNetwork/articles/introduction.html

## Network plotting on top of ggtree

- `tanggle` (`ggevonet`, `ggsplitnet`, `minimize_overlap`, `fortify.evonet`):
  https://klausvigo.github.io/tanggle/articles/tanggle_vignette.html ·
  https://github.com/KlausVigo/tanggle
- `ggret` (alternative reticulate-network plotting; `read_enewick`,
  `read_beast_retnet`, curved/"snake" reticulation edges):
  https://grdspcht.github.io/ggret/articles/intro_to_ggret.html
- ggtree book (layouts, `fortify`, coordinate data frame):
  https://yulab-smu.top/treedata-book/

## Split / consensus networks (related but out of scope)

- `phangorn::consensusNet()` / `networx` (implicit/split networks):
  https://klausvigo.github.io/phangorn/articles/Networx.html

## Network drawing theory (layout / crossing minimization)

- *Drawing tree-based phylogenetic networks with minimum number of crossings.*
  arXiv:2008.08960. https://arxiv.org/pdf/2008.08960
- *Algorithms for visualizing phylogenetic networks.* arXiv:1609.00755.
  https://arxiv.org/abs/1609.00755
- *Sketch, capture and layout phylogenies* (trees vs networks layout), bioRxiv 2025.
  https://www.biorxiv.org/content/10.1101/2025.04.01.646633v3.full

## In-repo reference

- `../plot_network.R` — the original single-network script (γ parsing + `ggevonet`).
