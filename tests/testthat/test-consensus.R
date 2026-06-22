# Phase 4: consensus & discrepancy encoding.

# Two networks, same backbone topology and same reticulation event (C is the
# hybrid, donor lineage = D); only the gammas differ.
n1 <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
n2 <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.7):1,(D:1,#H1:1::0.3):1);")
ns <- anansi_netset(list(n1, n2))

# A third with a conflicting backbone (A,C grouped) for support gradients.
n3 <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
ns3 <- anansi_netset(list(n1, n2, n3))

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
