# Phase 0 smoke tests: the package loads and the eNewick toolchain is wired up.
# Replaced/expanded by feature tests from Phase 1 onward.

test_that("the eNewick / evonet toolchain is available", {
  # A minimal extended-Newick network: B is a hybrid (#H1) with two parents.
  nwk <- "((A,(B)#H1),(C,#H1));"
  net <- ape::read.evonet(text = nwk)

  expect_s3_class(net, "evonet")
  expect_setequal(net$tip.label, c("A", "B", "C"))
  expect_equal(nrow(net$reticulation), 1L)
})
