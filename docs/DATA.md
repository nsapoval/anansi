# anansi — Input Data

The sample data in [`../data/`](../data) are sets of inferred phylogenetic networks
over a **fixed 15-taxon** lizard set (genera *Podarcis*, *Lacerta*, *Zootoca*), each
network stored as an extended-Newick (eNewick) string inside one column of a CSV.

## Files & schemas

The CSV column layouts differ between files (so `read_networks_csv()` must auto-detect
the network column). Verified schemas:

| File | Rows (networks) | Network column | Other columns |
|---|---|---|---|
| `d10k_x100_results.csv` | 60 | `network` | `dataset`, `reticulations`, `rank`, `log_probability` |
| `Lacerta_agilis_inferred_networks.csv` | 230 | `newick` | `index`, `num_reticulations`, `network_rank`, `log_pseudo_likelihood`, `starting_network` |
| `Zootoca_vivipara_inferred_networks.csv` | 230 | `newick` | (same as Lacerta) |
| `d10k_mostlen_inferred_networks.csv` | 460 | `newick` | `index`, `species`, `search_r`, `network_rank`, `num_reticulations`, `log_pseudo_likelihood`, `starting_network` |

Network-column name seen so far: `network` or `newick`. The per-network metadata
columns (rank, (pseudo-)likelihood, reticulation count) are carried into
`anansi_netset$meta` (one row per network).

Reticulation counts vary within a file (the discrepancy we want to surface). Example —
`d10k_mostlen_inferred_networks.csv`: 0 reticulations × 3, 1 × 165, 2 × 228, 3 × 64.

## Taxon set (fixed; 15 taxa)

All four files share the same 15 taxa (verified):

```
Lacerta_agilis           Podarcis_liolepis        Podarcis_raffonei
Podarcis_bocagei         Podarcis_melisellensis   Podarcis_siculus
Podarcis_cretensis       Podarcis_muralis         Podarcis_tiliguerta
Podarcis_erhardii        Podarcis_pityusensis     Podarcis_vaucheri
Podarcis_filfolensis     Podarcis_gaigeae         Zootoca_vivipara
```

`validate_taxon_set()` (Phase 1) confirms this holds across *every* network in a set,
not just the first row.

## Extended Newick (eNewick) format

Standard Newick describes the backbone tree (`TaxonName:branch_length`, nested
parentheses). Reticulations add **hybrid nodes** that appear **twice** in the string,
tagged `#Hk`:

- **1st occurrence** — attached after a closing `)`, on the **major (tree) edge**:
  `)#H1:<edge_len>::<major_gamma>`
- **2nd occurrence** — as a leaf-like clade on the **minor (reticulation) edge**:
  `#H1:<edge_len>::<minor_gamma>`

The two values after the double colon are **γ inheritance probabilities** (major +
minor ≈ 1.0). Example hybrid pair (from `d10k_x100_results.csv`, row 1):

```
)#H1:0.48679011719881643::0.7801372331743996      # major edge, γ≈0.78
#H1:2.052600271490827::0.21986276682560038        # minor edge, γ≈0.22
```

### Parsing notes (→ `R/io.R`)

- `ape::read.evonet(text = ...)` parses the **structure** into an `evonet`
  (`$edge` = backbone, `$reticulation` = hybrid edges, `$node.label` carries `#Hk`).
- `read.evonet()` does **not** populate γ. We extract γ from the raw string with a
  regex generalizing the reference `plot_network.R`:
  - find all `#H\d+:[^,()]+::[0-9.eE+-]+` matches,
  - for each label, 1st match = `tree_gamma` (major), 2nd = `ret_gamma` (minor),
  - map `#Hk` labels to internal node numbers via `$node.label`.
- Networks may have multiple hybrids (`#H1..#H3` in this data); the parser handles
  arbitrary `k`.

## Known parsing issues (`ape::read.evonet`)

Validating *every* network (not just the first row) surfaced two real
`read.evonet` limitations. anansi loads such files without crashing, flags the
affected networks, and lets you filter them.

1. **Phantom tips on nested/stacked reticulations.** When a hybrid placeholder is
   nested inside another hybrid's defining clade (level-k structures), `read.evonet`
   leaves an unlabeled placeholder node that surfaces as a tip with an empty (`""`)
   label — so the parsed network has 16–17 "taxa" instead of 15, and a node that is
   a multi-hybrid donor. The real 15-taxon set is unchanged. These are flagged via
   `anansi_network$issues` and the netset's `taxa_ok` vector. Observed incidence
   (it only hits some 2–3 reticulation networks):
   - `d10k_x100_results.csv`: 4 / 60
   - `d10k_mostlen_inferred_networks.csv`: 52 / 460
   - `Lacerta_agilis_*` and `Zootoca_vivipara_*`: **0** (fully clean; ≤2 reticulations)

   Keep only the consistent networks with `ns[ns$taxa_ok]`.

2. **Zero-reticulation networks fail `read.evonet`.** Plain trees (no `#H`) make
   `read.evonet` error/return a degenerate object (3 such networks in
   `d10k_mostlen`). `parse_network()` falls back to `ape::read.tree()` and attaches
   an empty reticulation matrix, so these parse cleanly.

A robust nested-reticulation parser (custom Rich-Newick parser, `ggret::read_enewick`,
or a `read.evonet` repair pass) is a **tracked roadmap item** — see
[`WORKPLAN.md`](WORKPLAN.md). Until then, `Lacerta`/`Zootoca` work end-to-end with no
flagged networks, and the `d10k` files work on their `taxa_ok` subset.

## Provenance

These are inference outputs (ranked candidate networks with (pseudo-)likelihood
scores). anansi only **visualizes** them; it does not infer networks.
