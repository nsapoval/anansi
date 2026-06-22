# Phase 1: I/O & data model.

nwk1 <- "((A:1,(B:1)#H1:1::0.7):1,(C:1,#H1:1::0.3):1);"          # taxa A,B,C; 1 hybrid
nwk2 <- "(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);"  # taxa A,B,C,D; 1 hybrid

test_that("parse_gammas extracts major/minor gammas per hybrid", {
  g <- parse_gammas(nwk1)
  expect_equal(nrow(g), 1L)
  expect_equal(g$label, "#H1")
  expect_equal(g$tree_gamma, 0.7)
  expect_equal(g$ret_gamma, 0.3)
})

test_that("parse_gammas returns zero rows when there are no hybrids", {
  g <- parse_gammas("((A:1,B:1):1,C:1);")
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 0L)
})

test_that("parse_network builds an anansi_network with a hybrid-node map", {
  net <- parse_network(nwk1, meta = list(rank = 1))
  expect_s3_class(net, "anansi_network")
  expect_s3_class(net$evonet, "evonet")
  expect_setequal(net$evonet$tip.label, c("A", "B", "C"))
  expect_equal(nrow(net$evonet$reticulation), 1L)
  expect_equal(length(net$hybrid_nodes), 1L)
  expect_equal(names(net$hybrid_nodes), "#H1")
  # the mapped node must be the recipient node in $reticulation
  expect_true(unname(net$hybrid_nodes) %in% net$evonet$reticulation[, 2])
})

test_that("read_networks_csv auto-detects the network column and builds a netset", {
  for (col in c("network", "newick", "enewick")) {
    df <- data.frame(idx = 1:2, score = c(-1, -2), stringsAsFactors = FALSE)
    df[[col]] <- c(nwk1, nwk1)
    f <- tempfile(fileext = ".csv")
    utils::write.csv(df, f, row.names = FALSE)

    ns <- read_networks_csv(f)
    expect_s3_class(ns, "anansi_netset")
    expect_equal(length(ns), 2L)
    expect_setequal(network_taxa(ns), c("A", "B", "C"))
    # metadata (non-network columns) is carried through
    expect_true(all(c("idx", "score") %in% names(ns$meta)))
    expect_false(col %in% names(ns$meta))
  }
})

test_that("read_networks_csv errors when no network column is detectable", {
  df <- data.frame(a = 1, b = 2)
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE)
  expect_error(read_networks_csv(f), "auto-detect")
})

test_that("validate_taxon_set accepts a fixed set; mismatch warns and flags", {
  same <- anansi_netset(list(parse_network(nwk1), parse_network(nwk1)))
  expect_true(validate_taxon_set(same))
  expect_true(all(same$taxa_ok))

  # A divergent taxon set is kept but flagged (warn, not error).
  expect_warning(
    mixed <- anansi_netset(list(parse_network(nwk1), parse_network(nwk1),
                                parse_network(nwk2))),
    "do not match the modal taxon set"
  )
  expect_equal(sum(mixed$taxa_ok), 2L)        # the two nwk1 are the modal set
  expect_false(mixed$taxa_ok[3])

  # Strict mode is available on demand.
  expect_error(validate_taxon_set(mixed, error = TRUE), "differ from the modal")
})

test_that("zero-reticulation (plain tree) networks parse via fallback", {
  net <- parse_network("((A:1,B:1):1,C:1);")
  expect_s3_class(net, "anansi_network")
  expect_s3_class(net$evonet, "evonet")
  expect_equal(nrow(net$evonet$reticulation), 0L)
  expect_length(net$issues, 0L)
  expect_setequal(net$evonet$tip.label, c("A", "B", "C"))
})

test_that("netset subsetting keeps class and aligns metadata", {
  df <- data.frame(idx = 1:3, newick = rep(nwk1, 3), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE)
  ns <- read_networks_csv(f)

  sub <- ns[2:3]
  expect_s3_class(sub, "anansi_netset")
  expect_equal(length(sub), 2L)
  expect_equal(sub$meta$idx, c(2L, 3L))
})
