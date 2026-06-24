# Phase 2: shared layout engine.

# Two 4-taxon networks (A,B,C,D), each with one reticulation, different topology.
na <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
nb <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
ns <- anansi_netset(list(na, nb))

test_that("consensus_tip_order returns a permutation of the taxa", {
  for (m in c("mode", "mds", "closest_leaf", "first")) {
    ord <- consensus_tip_order(ns, method = m)
    expect_setequal(ord, c("A", "B", "C", "D"))
  }
  expect_equal(consensus_tip_order(ns, method = "manual",
                                   order = c("D", "C", "B", "A")),
               c("D", "C", "B", "A"))
})

test_that("layout_network emits tree and reticulation segments", {
  ord <- c("A", "B", "C", "D")
  seg <- layout_network(na, tip_order = ord)
  expect_s3_class(seg, "data.frame")
  expect_setequal(unique(seg$kind), c("tree", "reticulation"))
  # one tree segment per backbone edge, one segment per reticulation
  expect_equal(sum(seg$kind == "tree"), nrow(na$evonet$edge))
  expect_equal(sum(seg$kind == "reticulation"), nrow(na$evonet$reticulation))
})

test_that("tips are pinned to their rank in the shared tip order", {
  ord <- c("C", "A", "D", "B")
  seg <- layout_network(na, tip_order = ord)
  tips <- seg[seg$is_tip, ]
  # a tip's y equals its position in tip_order
  expect_equal(tips$yend[match("A", tips$label)], match("A", ord))
  expect_equal(tips$yend[match("D", tips$label)], match("D", ord))
})

test_that("layout matches the network's own ggevonet ordering", {
  ord <- network_tip_order(na)            # the network's intrinsic order
  seg <- layout_network(na, tip_order = ord)
  tips <- seg[seg$is_tip, ]
  expect_equal(tips$label[order(tips$yend)], ord)
})

test_that("layout_netset shares tip y-positions across networks", {
  seg <- layout_netset(ns, method = "first")
  expect_true(".net" %in% names(seg))
  expect_length(unique(seg$.net), 2L)

  tips <- seg[seg$is_tip, c(".net", "label", "yend")]
  # same label -> same y in both networks
  wide <- tapply(tips$yend, list(tips$label, tips$.net), unique)
  expect_true(all(wide[, 1] == wide[, 2]))

  # x normalized to [0, 1] per network
  expect_equal(min(c(seg$x, seg$xend)), 0)
  expect_equal(max(c(seg$x, seg$xend)), 1)
})

test_that("color_by attaches a per-network .value aligned by .net", {
  nsv <- anansi_netset(list(na, nb), meta = data.frame(value = c(10, 20)))
  seg <- layout_netset(nsv, method = "first", color_by = "value")
  expect_true(".value" %in% names(seg))
  # every segment of network i carries that network's value
  by_net <- tapply(seg$.value, seg$.net, unique)
  expect_equal(as.numeric(by_net), c(10, 20))
})

test_that("color_by validates the column exists and is numeric", {
  nsv <- anansi_netset(list(na, nb),
                       meta = data.frame(value = c(1, 2), label = c("x", "y")))
  expect_error(layout_netset(nsv, method = "first", color_by = "missing"),
               "not a metadata column")
  expect_error(layout_netset(nsv, method = "first", color_by = "label"),
               "must be numeric")
})

test_that("divergent (taxa_ok = FALSE) networks are dropped with a warning", {
  od <- parse_network("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")  # 5 taxa -> divergent
  mixed <- suppressWarnings(anansi_netset(list(na, nb, od)))
  expect_warning(layout_netset(mixed, method = "first"), "Dropping 1 network")
})

test_that("outgroup pinning moves outgroup taxa to the chosen end", {
  top <- consensus_tip_order(ns, method = "mds", outgroup = "D",
                             outgroup_position = "top")
  expect_setequal(top, c("A", "B", "C", "D"))
  expect_equal(top[length(top)], "D")               # top = end of the vector
  bot <- consensus_tip_order(ns, method = "mds", outgroup = "D",
                             outgroup_position = "bottom")
  expect_equal(bot[1], "D")                          # bottom = start
  og <- consensus_tip_order(ns, method = "mds", outgroup = c("C", "D"),
                            outgroup_position = "top")
  expect_setequal(og[(length(og) - 1L):length(og)], c("C", "D"))  # contiguous
})

test_that("outgroup pinning errors on taxa not in the set", {
  expect_error(consensus_tip_order(ns, method = "first", outgroup = "Z"),
               "not in the taxon set")
})

test_that("snap_to_consensus yields a consensus-contiguous (crossing-free) order", {
  t1 <- parse_network("(((A:1,B:1):1,C:1):1,D:1);")
  t2 <- parse_network("(((A:1,B:1):1,C:1):1,D:1);")
  t3 <- parse_network("(((A:1,C:1):1,B:1):1,D:1);")
  nss <- anansi_netset(list(t1, t2, t3))
  ord <- consensus_tip_order(nss, method = "mds", snap_to_consensus = TRUE)
  expect_setequal(ord, c("A", "B", "C", "D"))
  tr <- consensus_network(nss)$tree
  dt <- descendant_tips(tr)
  ntip <- length(tr$tip.label)
  for (nd in (ntip + 1L):(ntip + tr$Nnode)) {
    pos <- sort(match(dt[[nd]], ord))
    # every consensus clade occupies a contiguous run -> no crossings
    expect_equal(pos, seq.int(pos[1], pos[length(pos)]))
  }
})
