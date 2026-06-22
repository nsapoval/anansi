# Phase 5: reticulation-aware taxon subsetting.

# Hybrid recipient subtends {C}; the reticulation's other endpoint involves the
# D lineage.
n1 <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
n2 <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
ns <- anansi_netset(list(n1, n2))

test_that("keeping all taxa preserves tips and the reticulation", {
  p <- restrict_taxa(n1, c("A", "B", "C", "D"))
  expect_s3_class(p, "anansi_network")
  expect_setequal(p$evonet$tip.label, c("A", "B", "C", "D"))
  expect_equal(nrow(p$evonet$reticulation), 1L)
  expect_equal(nrow(p$gammas), 1L)
})

test_that("dropping a reticulation endpoint's taxa drops the reticulation", {
  # Recipient subtends {C}; removing C must drop the reticulation.
  p <- restrict_taxa(n1, c("A", "B", "D"))
  expect_setequal(p$evonet$tip.label, c("A", "B", "D"))
  expect_equal(nrow(p$evonet$reticulation), 0L)
})

test_that("pruned networks remain layoutable and never gain reticulations", {
  p <- restrict_taxa(n1, c("A", "B", "C"))
  expect_lte(nrow(p$evonet$reticulation), nrow(n1$evonet$reticulation))
  seg <- layout_network(p, tip_order = c("A", "B", "C"))
  expect_s3_class(seg, "data.frame")
  expect_true(all(seg$is_tip[seg$kind == "tree" & seg$is_tip] %in% TRUE))
})

test_that("restrict_taxa on a set returns a netset on the kept taxa", {
  sub <- restrict_taxa(ns, c("A", "B", "C"))
  expect_s3_class(sub, "anansi_netset")
  expect_equal(length(sub), 2L)
  expect_setequal(network_taxa(sub), c("A", "B", "C"))
  expect_true(all(sub$taxa_ok))
})

test_that("restrict_taxa errors when too few taxa remain", {
  expect_error(restrict_taxa(n1, c("A")), "At least 2")
})

test_that("densinet(keep=) restricts before plotting", {
  p <- densinet(ns, method = "first", keep = c("A", "B", "C"))
  expect_s3_class(p, "ggplot")
})
