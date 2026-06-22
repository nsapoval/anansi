#' anansi: Consensus/Discrepancy Visualization for Sets of Phylogenetic Networks
#'
#' anansi overlays a *set* of explicit reticulate phylogenetic networks
#' (extended Newick / [ape::evonet]) over a fixed taxon set in a DensiTree-style
#' consensus/discrepancy view: networks are drawn on a shared tip ordering with
#' alpha transparency so that agreement appears dense and conflict appears
#' diffuse. Backbone-tree topology and reticulation-edge placement are shown
#' together, and all networks can be restricted to a subset of taxa.
#'
#' The package is organized as a small framework:
#' \itemize{
#'   \item I/O & data model (`io.R`, `netset.R`) -- read networks and build a
#'         validated `anansi_netset`.
#'   \item Layout engine (`layout.R`) -- compute a shared coordinate frame.
#'   \item Consensus (`consensus.R`) -- clade and reticulation-event frequencies.
#'   \item Plotting (`plot.R`, `discrepancy.R`) -- the `densinet()` overlay.
#'   \item Taxa (`taxa.R`) -- restrict networks to a subset of taxa.
#' }
#'
#' See `docs/DESIGN.md` for design decisions and `docs/WORKPLAN.md` for the
#' phased roadmap and feature-status tracker.
#'
#' @keywords internal
"_PACKAGE"
