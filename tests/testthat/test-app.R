# Web app helpers (run_anansi_app, build_densinet, netset_from_enewick).

test_that("run_anansi_app errors clearly when shiny is unavailable", {
  skip_if(requireNamespace("shiny", quietly = TRUE),
          "shiny is installed; cannot test the missing-shiny path")
  expect_error(run_anansi_app(), "shiny")
})

test_that("the bundled Shiny app.R parses", {
  app <- system.file("shiny", "app.R", package = "anansi")
  skip_if(!nzchar(app), "bundled app not found (run devtools::load_all/install)")
  expect_silent(parse(app))   # syntactic validity, without executing the app
})

test_that("build_densinet maps params to a ggplot from the example", {
  ns <- read_networks_csv(anansi_example())
  p <- build_densinet(ns, list(method = "first",
                               reticulation_style = "hybrid",
                               consensus_p = 0.5))
  expect_s3_class(p, "ggplot")
})

test_that("build_densinet forces outgroup into keep (overlap handling)", {
  ns <- read_networks_csv(anansi_example())
  tx <- network_taxa(ns)
  # outgroup deliberately NOT in keep: build_densinet must union it back in so
  # the outgroup is not pruned away (which would break tip-order pinning).
  p <- build_densinet(ns, list(method = "first",
                               outgroup = tx[1], keep = tx[2:5]))
  expect_s3_class(p, "ggplot")
})

test_that("build_densinet rejects a non-netset", {
  expect_error(build_densinet(list()), "anansi_netset")
})

test_that("netset_from_enewick parses pasted strings into a netset", {
  txt <- paste(
    "((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);",
    "# a comment line is ignored",
    "((A:1,(C:1)#H1:1::0.6):1,(B:1,#H1:1::0.4):1);",
    sep = "\n")
  ns <- netset_from_enewick(txt)
  expect_s3_class(ns, "anansi_netset")
  expect_equal(length(ns), 2L)
  expect_setequal(network_taxa(ns), c("A", "B", "C"))
  expect_s3_class(build_densinet(ns, list(method = "first")), "ggplot")
})

test_that("netset_from_enewick errors on empty input", {
  expect_error(netset_from_enewick("\n  \n# only comments"), "No extended-Newick")
})
