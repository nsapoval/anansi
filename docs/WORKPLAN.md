# anansi — Workplan & Feature Tracker

> Living roadmap. Update the **status** column as features land. Design rationale
> lives in [`DESIGN.md`](DESIGN.md); this file is *what* and *when*.

Status legend: ✅ done · 🚧 in progress · ⬜ planned

## Phases

### Phase 0 — Scaffolding & docs 🚧
Objective: a loadable R-package skeleton plus living documentation.
- ✅ Package skeleton: `DESCRIPTION`, `NAMESPACE` (roxygen), `R/` module files,
  `.Rbuildignore`, `tests/testthat` harness.
- ✅ `devtools::load_all()` succeeds; smoke test passes.
- ✅ Docs authored: `DESIGN.md`, `WORKPLAN.md`, `REFERENCES.md`, `DATA.md`, `README.md`.

*Acceptance:* `devtools::load_all(".")` and `devtools::test()` run clean. **Met.**
(`R CMD check`: 0 errors / 0 warnings / 1 NOTE — "declared Imports not used", expected
while module files are stubs; resolves as Phase 1+ add code that uses each package.)

### Phase 1 — I/O & data model ✅
Objective: read network sets and build a validated data model.
- ✅ `read_networks_csv(path, network_col=NULL)` — auto-detects the network column
  (`network`/`newick`/`enewick`).
- ✅ `parse_network(enewick, meta)` — robust evonet read (with plain-tree fallback)
  + generalized γ parsing (N hybrids #H1..#Hk) + hybrid-node→label map + `issues`.
- ✅ `anansi_network` / `anansi_netset` constructors, `print`/`[`/`length`.
- ✅ `network_taxa()`, `validate_taxon_set()`, `taxa_consistency()` (modal-set check
  with per-network `taxa_ok` flagging — warns rather than erroring).
- ✅ `plot_network(net)` — single-network plot via the new API.

*Acceptance:* re-plotting network 1 of `d10k_x100_results.csv` through the API
reproduces the original `plot_network.R` PNG (verified, pixel-faithful); 46 unit
tests pass; `R CMD check` is clean (0/0/0). **Met.**

*Finding (see [`DATA.md`](DATA.md)):* `ape::read.evonet` mis-parses nested/stacked
reticulations (phantom empty-label tips) and fails on zero-reticulation trees. The
data model now flags these instead of crashing. `Lacerta`/`Zootoca` are 100% clean;
the `d10k` files carry a flagged minority. A robust parser is tracked below (W1).

### Phase 2 — Shared layout engine ✅
Objective: a shared coordinate frame so networks overlay.
- ✅ `consensus_tip_order(netset, method=c("mode","mds","closest_leaf","first","manual"))`
  (`"consensus"` deferred; the four heuristics + manual are implemented).
- ✅ `layout_network(net, tip_order, mode=c("cladogram","phylogram"))` → tidy
  segments (`from, to, x, y, xend, yend, kind ∈ {tree,reticulation}, is_tip, label`).
- ✅ `layout_netset(netset, tip_order, method, mode, scale_x, consistent_only)` →
  segments + `.net` id; x normalized to `[0,1]`; `taxa_ok=FALSE` networks dropped.

*Acceptance:* tip-y is identical per label across networks (tested + verified on 230
Lacerta networks); a single network's layout matches its own `ggevonet` ordering
(tested); 64 unit tests pass; engine runs in ~0.8s on 230 networks. Overlay preview
on real data shows a coherent consensus cloud + reticulation cloud. **Met.**

### Phase 3 — `densinet()` overlay MVP ✅
Objective: the DensiTree-for-networks overlay.
- ✅ `densinet(x, tip_order, method, mode, layout, tree_color, ret_color, alpha,
  ret_alpha, ret_linetype, tip_labels, ...)`.
- ✅ backbone cloud (alpha) + reticulation cloud (distinct color, dashed) + tip
  labels on the shared order; auto alpha (~1/N, ~3/N).
- ✅ `"slanted"` and `"rectangular"` backbone layouts; returns a composable `ggplot`.

*Acceptance:* `densinet()` renders a readable consensus cloud + reticulation cloud on
`Lacerta_agilis` (clean consensus) and on the `d10k_mostlen` consistent subset (dense
reticulation cloud, 1-3 reticulations). 71 tests pass; `R CMD check` clean. **Met.**

*Note:* the layout engine was reimplemented on `ape` directly (no `fortify.evonet`)
because that S3 method is not reliably found when tanggle is imported-but-not-attached
(under `R CMD check`). Coordinates are computed from the backbone tree.

### Phase 4 — Consensus & discrepancy encoding ✅
Objective: make agreement vs conflict legible.
- ✅ `clade_frequencies()`, `reticulation_frequencies()` (events keyed by
  donor/recipient clades, with mean gamma; DESIGN D4).
- ✅ `consensus_network(netset, p)` -> `anansi_consensus` ("root-canal" backbone
  with per-clade support + frequent reticulation events).
- ✅ `densinet(consensus=TRUE)` overlay: solid consensus backbone colored by clade
  support + consensus reticulation arrows (donor->recipient, frequency-thresholded
  via `consensus_ret_min`).
- ✅ `top_networks()` / `densinet(top_n=, top_by=)` to focus on the best networks.
- ⬜ Top-N *topology* coloring (DensiTree multi-color of the cloud) -> Phase 6.

*Acceptance:* the dominant topology (dark, support-colored consensus backbone) and
frequent reticulations (red arrows) are clearly distinct from the faint cloud on
`Lacerta_agilis` and the `d10k_mostlen` subset (see figures/07_*). 89 tests pass;
`R CMD check` clean (0/0/0). **Met.**

### Phase 5 — Taxon subsetting ✅
Objective: restrict all networks to a subset of taxa.
- ✅ `restrict_taxa(x, keep)` (set or single network) — reticulation-aware
  pruning: backbone `ape::keep.tip` + degree-2 suppression, then each
  reticulation is re-located by the MRCA of its surviving donor/recipient tips
  and dropped if an endpoint vanishes, the endpoints collapse, or one clade
  nests in the other (DESIGN D7). Reticulation identity is tip-set based, so no
  fragile node-number tracking through pruning.
- ✅ `densinet(x, keep = ...)` restricts before plotting.

*Acceptance:* restricting Lacerta to 9 taxa yields a valid netset (230/230
taxa_ok), correctly drops reticulations whose endpoints were removed (1-2/net ->
0-1/net), and produces a clean overlay (figures/08_lacerta_restricted.png). 104
tests pass; `R CMD check` clean (0/0/0). **Met.**

### Phase 6 — Polish ✅ (mostly; two items deferred with reasons)
Objective: production quality.
- ✅ Phylogram mode (`densinet(mode="phylogram")`): per-network branch-length x
  scaled to [0,1]; consensus drawn with mean clade heights across networks.
- ✅ Jitter (`densinet(jitter=)`, `layout_netset(jitter=)`).
- ✅ Performance: ~1s for 230 networks; layout/consensus scale to the 460-network
  files without trouble.
- ✅ Legends / themes / scales: clade-support gradient legend, white background,
  `theme_void`; auto alpha defaults.
- ✅ Bundled example data (`inst/extdata/lacerta_sample.csv`) + `anansi_example()`
  so examples/tutorial run in the installed package.
- ✅ Fuller test suite (114 tests); `R CMD check` clean (0/0/0).
- ◐ Walkthrough: `docs/TUTORIAL.md` (verified) instead of an Rmd vignette --
  **pandoc is unavailable** in this environment, so a formal
  `vignette("anansi")` is deferred until it is (the tutorial covers the same
  ground).
- ⬜ DensiTree-style top-N *topology* coloring of the cloud: deferred -- it needs a
  discrete colour scale that conflicts with the consensus support gradient (would
  require `ggnewscale`); the support-colored consensus backbone already conveys
  the dominant topology.

## Feature tracker (cross-cutting)

| Feature | Status | Phase | Notes |
|---|---|---|---|
| Read eNewick network sets from CSV | ✅ | 1 | auto-detect `network`/`newick`/`enewick` |
| Parse γ inheritance probabilities (N hybrids) | ✅ | 1 | generalized `plot_network.R` regex |
| `anansi_network` / `anansi_netset` model | ✅ | 1 | + `issues`, `taxa_ok` |
| Fixed-taxon-set validation | ✅ | 1 | modal-set check, warn + flag |
| Single-network plot (parity with reference) | ✅ | 1 | verified pixel-faithful |
| **Robust nested-reticulation parser (W1)** | ✅ | — | native parser; all files 100% |
| Consensus tip ordering | ✅ | 2 | mode/mds/closest_leaf/first/manual |
| Shared cladogram layout engine | ✅ | 2 | + phylogram path; x scaled to [0,1] |
| `densinet()` overlay (backbone + reticulation cloud) | ✅ | 3 | slanted + rectangular |
| x-jitter | ⬜ | 3 | deferred to Phase 6 polish |
| Clade & reticulation frequencies | ✅ | 4 | event identity = donor/recipient clades |
| Consensus / root-canal overlay | ✅ | 4 | support-colored backbone + ret arrows |
| Top-N network filtering (`top_n`) | ✅ | 4 | by likelihood/rank metadata |
| Top-N topology coloring (cloud) | ⬜ | 6 | DensiTree multi-color, deferred |
| Reticulation-conflict highlight | ✅ | 4 | freq-thresholded consensus arrows |
| Restrict to taxon subset (pruning) | ✅ | 5 | reticulation-aware (MRCA re-location) |
| Phylogram mode + normalization | ✅ | 6 | per-network [0,1]; consensus mean heights |
| Jitter | ✅ | 6 | `densinet(jitter=)` |
| Bundled example data + `anansi_example()` | ✅ | 6 | inst/extdata |
| Walkthrough (docs/TUTORIAL.md) | ✅ | 6 | Rmd vignette deferred (no pandoc) |

## Backlog / decisions

**W1 — Robust nested-reticulation eNewick parsing ✅ (done).** Implemented a
self-contained recursive-descent extended-Newick parser (`R/enewick.R`,
`enewick_to_anansi()`), now the **primary** parser in `parse_network()`
(`ape::read.evonet` kept only as a fallback). It merges each hybrid's placeholder
into its defining occurrence as a reticulation edge, then prunes nodes with no
sampled descendants (ghost donors -> re-point reticulation to parent; leaf hybrids
-> drop the reticulation). Resolves nested/stacked (level-k) reticulations and
hybrid-free trees. Chosen over `ggret` (GitHub-only dep) and over repairing
`read.evonet` output (fragile). Revises design D8 (we now own the parser).
*Result:* all sample files parse 100% (`taxa_ok` = 60/60, 460/460, 230/230,
230/230); 140 tests; `R CMD check` clean.

## Non-goals (current)

- Network inference. · Split/implicit (`networx`) networks. · Interactive GUI. ·
  Provably crossing-minimal layout (NP-hard; heuristics + documented limits).
