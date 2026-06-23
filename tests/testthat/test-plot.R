# Phase 3: densinet() overlay.

na <- parse_network("(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
nb <- parse_network("(((A:1,C:1):1,(B:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);")
ns <- anansi_netset(list(na, nb))

test_that("densinet returns a composable ggplot with the expected layers", {
  p <- densinet(ns, method = "first")
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(L) class(L$geom)[1], character(1))
  expect_true("GeomSegment" %in% geoms)  # tree + reticulation clouds
  expect_true("GeomText" %in% geoms)     # tip labels
})

test_that("densinet builds for both layouts and without tip labels", {
  expect_s3_class(densinet(ns, method = "first", layout = "rectangular"), "ggplot")
  p <- densinet(ns, method = "first", tip_labels = FALSE)
  geoms <- vapply(p$layers, function(L) class(L$geom)[1], character(1))
  expect_false("GeomText" %in% geoms)
})

test_that("densinet renders to a file without error", {
  p <- densinet(ns, method = "first")
  f <- tempfile(fileext = ".png")
  ggplot2::ggsave(f, p, width = 5, height = 4, dpi = 72)
  expect_true(file.exists(f))
})

test_that("densinet rejects non-netset input", {
  expect_error(densinet(na), "anansi_netset")
})

test_that("densinet supports outgroup, snap, and hybrid reticulation style", {
  expect_s3_class(densinet(ns, method = "mds", outgroup = "D",
                           outgroup_position = "top"), "ggplot")
  expect_s3_class(densinet(ns, method = "mds", snap_to_consensus = TRUE), "ggplot")
  p <- densinet(ns, method = "first", reticulation_style = "hybrid")
  expect_s3_class(p, "ggplot")
  # hybrid style draws dotted edges, no arrowhead layer
  arrows <- vapply(p$layers, function(L) !is.null(L$geom_params$arrow), logical(1))
  expect_false(any(arrows))
})
