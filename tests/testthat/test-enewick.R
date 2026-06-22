# W1: native extended-Newick parser (nested/stacked reticulations).

test_that("native parser reads a clean network with gammas", {
  net <- parse_network("((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);")
  expect_setequal(net$evonet$tip.label, c("A", "B", "C"))
  expect_equal(nrow(net$evonet$reticulation), 1L)
  expect_equal(net$gammas$tree_gamma, 0.7)   # major (defining) edge
  expect_equal(net$gammas$ret_gamma, 0.3)    # minor (placeholder) edge
  expect_length(net$issues, 0L)              # native parse, no fallback note
  # the recipient (to) is the hybrid defining node
  expect_true(unname(net$hybrid_nodes) %in% net$evonet$reticulation[, 2])
})

test_that("native parser handles hybrid-free trees (read.evonet fails on these)", {
  net <- parse_network("((A:1,B:1):1,C:1);")
  expect_s3_class(net$evonet, "evonet")
  expect_equal(nrow(net$evonet$reticulation), 0L)
  expect_setequal(net$evonet$tip.label, c("A", "B", "C"))
})

test_that("native parser resolves nested/stacked reticulations without phantom tips", {
  # The real failure cases that leave read.evonet with phantom 16th tips.
  fix <- read_networks_csv(anansi_example("nested_reticulations.csv"),
                           network_col = "network")
  for (net in fix$networks) {
    tl <- net$evonet$tip.label
    expect_false(any(is.na(tl) | !nzchar(tl)))     # no phantom/empty tips
    expect_false(any(grepl("^#", tl)))             # no leftover #H tips
    expect_equal(length(tl), 15L)                  # exactly the 15 taxa
    # backbone must be a valid ape tree
    expect_silent(ape::reorder.phylo(
      structure(`[[<-`(net$evonet, "reticulation", NULL), class = "phylo"),
      "postorder"))
  }
  expect_true(all(fix$taxa_ok))                    # all consistent now
})
