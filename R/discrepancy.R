# discrepancy.R -------------------------------------------------------------
#
# Discrepancy encodings: mapping consensus support onto visual channels.
#
# Planned API (Phase 4):
#   scale_by_support(...)        -> map clade/reticulation frequency to
#                                   opacity / color / line width
#   highlight_topologies(p, netset, top_n = 3)  -> color the top-N topologies
#   highlight_reticulation_conflict(p, netset)  -> emphasize conflicting rets
#
# Design (see docs/DESIGN.md):
#   * "Emphasize both signals equally": backbone clade support AND reticulation-
#     event support both drive the discrepancy encoding.
#   * Top-N topology coloring follows DensiTree (dominant topology darkest /
#     primary color; alternatives in secondary colors).
#
# (Implemented in Phase 4.)
