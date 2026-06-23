# Phase 4: consensus & discrepancy encoding.

# Two networks, same backbone topology and same reticulation event (C is the
# hybrid, donor lineage = D); only the gammas differ.
n1 <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
n2 <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
ns <- anansi_netset(list(n1, n2))

# A third with a conflicting backbone (A,C grouped) for support gradients.
n3 <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
ns3 <- anansi_netset(list(n1, n2, n3))

# A direction-unstable pair: same hybrid B, but its two parent lineages (A and C)
# swap roles of backbone-sibling vs donor. Identical undirected, opposite directed.
u1 <- parse_network("((A,(B)#H1),(C,#H1));")
u2 <- parse_network("((C,(B)#H1),(A,#H1));")
uns <- anansi_netset(list(u1, u2))

test_that("clade_frequencies returns valid frequencies", {
  cf <- clade_frequencies(ns)
  expect_s3_class(cf, "data.frame")
  expect_true(all(cf$freq > 0 & cf$freq <= 1))
  expect_true("A,B" %in% cf$clade)          # the (A,B) clade exists
  expect_equal(cf$freq[cf$clade == "A,B"], 1)  # in both networks

  cf3 <- clade_frequencies(ns3)
  expect_equal(cf3$freq[cf3$clade == "A,B"], 2 / 3)  # only 2 of 3 networks
})

test_that("reticulation_frequencies keys events by donor/recipient clades", {
  rf <- reticulation_frequencies(ns)
  expect_equal(nrow(rf), 1L)                 # one shared event
  expect_equal(rf$count, 2L)
  expect_equal(rf$freq, 1)
  expect_equal(rf$mean_gamma, mean(c(0.4, 0.3)))
})

test_that("consensus_network returns a tree with support and reticulations", {
  cons <- consensus_network(ns3, p = 0.5)
  expect_s3_class(cons, "anansi_consensus")
  expect_s3_class(cons$tree, "phylo")
  expect_true(all(cons$support >= 0 & cons$support <= 1))
  expect_true(nrow(cons$reticulations) >= 1L)
})

test_that("top_networks selects by a metadata column", {
  df <- data.frame(
    score = c(-10, -1, -5),
    newick = c("((A:1,B:1):1,C:1);", "((A:1,B:1):1,C:1);", "((A:1,B:1):1,C:1);"),
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".csv"); utils::write.csv(df, f, row.names = FALSE)
  nset <- read_networks_csv(f)
  top <- top_networks(nset, 1, by = "score")   # highest score = -1 (row 2)
  expect_equal(length(top), 1L)
  expect_equal(top$meta$score, -1)
})

test_that("densinet draws the consensus overlay", {
  p <- densinet(ns3, method = "first", consensus = TRUE)
  expect_s3_class(p, "ggplot")
  f <- tempfile(fileext = ".png")
  ggplot2::ggsave(f, p, width = 5, height = 4, dpi = 72, bg = "white")
  expect_true(file.exists(f))

  # consensus = FALSE still works
  expect_s3_class(densinet(ns3, method = "first", consensus = FALSE), "ggplot")
})

test_that("directed = FALSE merges direction-swapped reticulation events", {
  d <- reticulation_frequencies(uns, directed = TRUE)
  expect_equal(nrow(d), 2L)                 # two distinct directed events (A->B, C->B)
  u <- reticulation_frequencies(uns, directed = FALSE)
  expect_equal(nrow(u), 1L)                 # collapsed into one
  expect_equal(u$count, 2L)
  expect_equal(u$freq, 1)
  expect_equal(u$recipient, "B")
  expect_setequal(c(u$parent1, u$parent2), c("A", "C"))
})

test_that("ret_edge_frac anchors reticulations off the node", {
  ord <- consensus_tip_order(ns, method = "first")
  cons <- consensus_network(ns)
  cs0 <- consensus_segments(cons, ord, ret_edge_frac = 0)
  csf <- consensus_segments(cons, ord, ret_edge_frac = 0.3)
  cols <- c("x", "y", "xend", "yend")
  expect_false(isTRUE(all.equal(cs0$reticulations[, cols],
                                csf$reticulations[, cols])))
})

test_that("undirected hybrid consensus_segments returns two edges per event", {
  ord <- consensus_tip_order(uns, method = "first")
  cs <- consensus_segments(consensus_network(uns, directed = FALSE), ord)
  expect_equal(nrow(cs$reticulations), 2L)            # one event -> two dotted edges
  # both edges meet at the recipient (shared start point)
  starts <- unique(paste(cs$reticulations$x, cs$reticulations$y))
  expect_length(starts, 1L)
})
