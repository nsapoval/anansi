# anansi — Design

> Living design document. Records *why* we make each choice. Pairs with
> [`WORKPLAN.md`](WORKPLAN.md) (what we build and in what order, with status),
> [`DATA.md`](DATA.md) (input formats), and [`REFERENCES.md`](REFERENCES.md).

## 1. Problem & goal

Phylogenetic network inference produces **sets** of candidate networks (e.g. ranked
by pseudo-likelihood). A single network plot tells you nothing about how stable that
estimate is. We want to look at a whole set at once and read off:

- **Consensus** — where do the networks agree (backbone topology *and* reticulation
  placement)?
- **Discrepancy** — where do they conflict, and how badly?

The taxon set is assumed **fixed** across all networks plotted together. We also want
to **restrict** all networks to a subset of taxa and re-view.

The inspiration is **DensiTree** (Bouckaert 2010, 2014), which does exactly this for
sets of *trees*: it overlays all trees with alpha-transparency on one shared tip
ordering, so consensus shows as dark/dense regions and conflict as a diffuse "web".
Our job is to adapt that idea to **explicit reticulate networks** (DAGs with
hybridization edges), which breaks several of DensiTree's assumptions (see §8).

## 2. Scope

**In scope:** explicit reticulate networks given as extended Newick (eNewick),
represented as `ape::evonet`; a fixed taxon set; cladogram-first overlay; consensus
& discrepancy encoding for both backbone and reticulations; taxon subsetting.

**Out of scope (for now):** network *inference*; split / implicit networks
(`phangorn::networx`); interactive GUI; guaranteed crossing-minimal layout
(NP-hard — we use heuristics). See [`WORKPLAN.md`](WORKPLAN.md) "Non-goals".

## 3. The DensiTree playbook (and what we borrow)

DensiTree's method, distilled:

1. **Shared 1-D tip order.** Pick one ordering of taxa once; reorder every tree's
   tips to match. This is what lets all trees overlay on a common (time × tip) plane.
   The order is chosen by a heuristic (e.g. closest-leaf greedy, or DFS of a
   consensus/backbone tree) to minimize crossings/clutter — order quality dominates
   readability.
2. **Alpha overlay.** Draw all trees with alpha ≈ 1/N. Agreement accumulates to
   opaque; disagreement stays faint. (Density = consensus.)
3. **Consensus / "root canal".** Group trees by topology; for each topology draw a
   representative with opacity ∝ its frequency, and/or draw a single consensus tree
   on top. Top-N topologies get distinct colors.
4. **Cladogram vs phylogram.** Cladogram = topology only (uniform depths); phylogram
   = branch-length/time scaled with tips aligned and node heights averaged.

We borrow **all four**. The adaptation to networks is in §4–§6.

## 4. Core design decisions

### D1 — Backbone-tree + reticulation-overlay decomposition *(the key idea)*
An `evonet` already separates the two parts we need: `$edge` is a **tree** (the
backbone), and `$reticulation` is a separate matrix of hybrid edges. We treat every
network as **backbone tree + reticulation edges**:

- The **backbone tree** drives tip ordering and the cladogram/phylogram coordinates —
  this is where DensiTree's tree machinery applies directly.
- The **reticulation edges** are a *separate visual layer* drawn on the same
  coordinate frame.

This sidesteps the fact that a DAG has no unique root-to-tip path: the tip order is
defined by the (tree) backbone, and reticulations are overlaid as arcs/segments.
*Caveat:* the backbone is the tree implied by `$edge`; the choice of which parent edge
is "the tree edge" vs. "the reticulation edge" comes from the eNewick (the `#Hk` minor
edge is the reticulation). We adopt that convention and document it.

### D2 — Shared coordinate frame pinned to a consensus tip order
All networks are laid out on **one** tip order, so tips share y-positions and the
networks superimpose. Internal-node positions are recomputed per network on that
order. The default layout is **cladogram** (uniform node depths, matching the
reference `plot_network.R`); phylogram is a later mode (§6, D6).

Coordinate spec (cladogram, the default):
- **y(tip)** = rank of the tip in the shared tip order (1..n), identical across all
  networks. **y(internal)** = mean of its children's y (slanted/“V” style), recomputed
  per network.
- **x(node)** = depth-based. Tips are **aligned** at a common x (right side) by
  default; root at the left. (Uniform depth = "node_depth"-style, as the reference
  script obtains by setting `edge.length <- NULL` before `ggevonet`.)
- A small optional **x-jitter** per network (à la DensiTree) separates exactly
  overlapping edges.

### D3 — Consensus tip ordering
`consensus_tip_order(netset, method=...)` over the backbone trees, methods:
- `"mode"` — most frequent leaf order (ggdensitree's default idea);
- `"mds"` — MDS on a leaf-distance matrix;
- `"closest_leaf"` — DensiTree's greedy closest-leaf heuristic on the averaged
  pairwise leaf distance (clade-span distance);
- `"consensus"` — DFS order of a majority-rule consensus backbone tree;
- `"manual"` — user-supplied vector.

We reuse `ggtree::ggdensitree` ordering ideas and `ape`/`phangorn` consensus helpers
rather than reimplementing from scratch where possible.

### D4 — Reticulation event identity & consensus
To tally "how often does this reticulation occur across the set", reticulations need a
**canonical identity**. We key a reticulation event by the **pair of backbone clades**
it connects: the clade below the *recipient* (hybrid) node and the clade below the
*donor* node, each represented as a taxon-set bipartition. Two reticulations in
different networks are "the same event" iff their (donor-clade, recipient-clade) keys
match. Event frequency drives color/opacity. The ambiguities here are real and tracked
in §9 (Open questions).

### D5 — Emphasize backbone and reticulations equally
The discrepancy encoding maps support onto visual channels for **both** layers:
backbone clade/split frequency *and* reticulation-event frequency drive opacity/color.
Top-N backbone topologies can be colored distinctly (DensiTree-style). Conflicting
reticulations (events present in only a minority) are highlightable.

### D6 — Output is a composable ggplot/ggtree object
`densinet()` returns a `ggplot` object so users can `+ geom_tiplab()`, `+ theme_*()`,
etc., consistent with `ggtree`/`tanggle` and the existing script.

### D7 — Reticulation-aware taxon pruning
`restrict_taxa(netset, keep)` (Phase 5): prune the backbone with `ape::keep.tip`, then
resolve each reticulation — keep it iff both endpoints survive (suppressing resulting
degree-2 nodes), otherwise drop it. Document edge cases (a hybrid that loses one
parent becomes a normal node; a reticulation that collapses to a parallel edge is
dropped). Pruning happens on the data model, before layout, so every view (single or
overlay) can be computed on the restricted set.

### D8 — Reuse, don't reinvent
Use `ape` (`read.evonet`, `as.phylo`, `keep.tip`, clade/bipartition utilities),
`tanggle` (`fortify.evonet`, `ggevonet`, `minimize_overlap`), `ggtree` (geoms,
`ggdensitree` ordering), optionally `phangorn` (consensus). New code only for: the
multi-network overlay, reticulation consensus, and reticulation-aware pruning.

## 5. Data model

```
anansi_network            # one network
  $ evonet  : ape::evonet # backbone in $edge, reticulations in $reticulation
  $ gammas  : data.frame  # label, tree_gamma, ret_gamma (parsed from eNewick)
  $ meta    : list        # log-likelihood, rank, num_reticulations, source row, ...

anansi_netset             # a set of anansi_network over a fixed taxon set
  $ networks : list<anansi_network>
  $ taxa     : character  # the shared, validated taxon set
  $ meta     : data.frame # one row per network (the source CSV columns)
```

Rationale: keep the trusted `ape::evonet` as the substrate (so all of `ape`/`tanggle`
keeps working) and attach the things `read.evonet()` drops (gammas) plus provenance.

## 6. Layout modes

- **Cladogram (default).** Topology only; robust to the noisy inferred branch lengths
  (which range ~0.1–10 in the sample data). Best for reading topological consensus.
- **Phylogram (Phase 6).** Branch-length/time scaled. Requires cross-network
  normalization (e.g. scale each network to common root height, or rank-normalize
  node heights) because absolute lengths are not comparable across networks. Node-
  height uncertainty then shows as a blur, as in DensiTree.

## 7. Rendering layers (the overlay)

1. **Backbone cloud** — all backbone trees, alpha-blended, one color.
2. **Reticulation cloud** — all reticulation edges, alpha-blended, a distinct color;
   drawn as curved/“snake” segments to distinguish from tree edges.
3. **Consensus overlay (optional)** — the consensus backbone ("root canal") and/or the
   most-frequent reticulation events, opaque and colored by support.
4. **Tip labels** — once, in the shared order.

## 8. Why networks break DensiTree's assumptions (and our responses)

| DensiTree assumption | Why networks break it | anansi response |
|---|---|---|
| Unique root→tip path per taxon | Reticulation nodes have ≥2 parents | Order from the **backbone tree**; overlay reticulations separately (D1) |
| One 1-D tip order minimizes crossings | Reticulation edges connect distant clades; a single order can't make them all short | Accept some long reticulation arcs; minimize via backbone ordering + `minimize_overlap`; document as a known limitation |
| Majority-rule consensus tree is well-defined | Clades aren't uniquely defined in a DAG | Compute consensus on the **backbone**; tally reticulations as separate events (D4) |
| Branch lengths comparable across trees | Network branch lengths vary wildly / units differ | Cladogram default; phylogram only with explicit normalization (D6) |

## 9. Open questions (revisit as we build)

- **Reticulation identity (D4).** Is (donor-clade, recipient-clade) the right key? How
  to handle near-misses (clades that differ by one taxon)? Should we cluster events
  rather than require exact key matches?
- **γ direction / donor vs recipient.** eNewick marks the minor edge; is the donor the
  minor-edge parent? Confirm against the inference tool's convention for this data.
- **Level-k overlap.** Networks with 2–3 reticulations may have interacting events;
  does per-event tallying mislead when reticulations are nested?
- **Pruning semantics (D7).** Exact rule when a hybrid loses a parent, or when both a
  reticulation's endpoints collapse onto the same backbone edge.
- **Phylogram normalization (D6).** Which normalization preserves interpretable node-
  height uncertainty without letting one network dominate?
- **Parser robustness (W1) — RESOLVED.** `ape::read.evonet` mis-parses
  nested/stacked (level-k) reticulations (phantom tips) and fails on
  zero-reticulation trees. We now use a **native recursive-descent extended-Newick
  parser** (`R/enewick.R`, `enewick_to_anansi()`) as the primary parser, with
  `read.evonet` as a fallback. It merges hybrid placeholders into reticulation
  edges and prunes nodes with no sampled descendants (ghost donors / leaf hybrids).
  This revises D8: we own the eNewick parser rather than reusing `read.evonet`.
  Remaining nuance: leaf-hybrid reticulations (recipient with no sampled
  descendants) are dropped as invisible over the taxon set — reasonable for
  visualization, but worth noting if exact reticulation counts matter.

## 10. Dependencies & environment

R ≥ 4.1. `Imports`: ape, tanggle, ggtree, ggplot2, tidytree, dplyr, rlang.
`tanggle` and `ggtree` are **Bioconductor** packages — install via `BiocManager`.
Dependency pinning (`renv`) is optional and may be added later; not blocking.
