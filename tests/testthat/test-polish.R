# Phase 6: bundled data, phylogram mode, jitter.

test_that("anansi_example ships a readable dataset", {
  path <- anansi_example()
  expect_true(file.exists(path))
  ns <- read_networks_csv(path)
  expect_s3_class(ns, "anansi_netset")
  expect_gt(length(ns), 1L)
  expect_equal(length(network_taxa(ns)), 15L)
})

test_that("phylogram mode lays out and plots (incl. consensus overlay)", {
  ns <- read_networks_csv(anansi_example())
  seg <- layout_netset(ns, method = "first", mode = "phylogram")
  expect_true(all(c("x", "xend") %in% names(seg)))
  expect_true(all(is.finite(c(seg$x, seg$xend))))
  p <- densinet(ns, method = "first", mode = "phylogram")  # consensus = TRUE
  expect_s3_class(p, "ggplot")
})

test_that("jitter offsets networks horizontally", {
  ns <- read_networks_csv(anansi_example())
  s0 <- layout_netset(ns, method = "first", jitter = 0)
  sj <- layout_netset(ns, method = "first", jitter = 0.02)
  expect_false(isTRUE(all.equal(s0$x, sj$x)))
})

test_that("consensus_network carries per-node heights for phylogram", {
  ns <- read_networks_csv(anansi_example())
  cons <- consensus_network(ns)
  expect_length(cons$heights, length(cons$tree$tip.label) + cons$tree$Nnode)
  expect_true(any(is.finite(cons$heights)))
})
