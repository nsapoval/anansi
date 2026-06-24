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

test_that("densinet colors the cloud by a per-tree value (viridis + gradient)", {
  trees <- lapply(c("((A,B),C);", "((A,C),B);", "((B,C),A);"), parse_network)
  nsv <- anansi_netset(trees, meta = data.frame(value = c(0, 1, 2)))

  p <- densinet(nsv, color_by = "value")
  expect_s3_class(p, "ggplot")
  # a continuous color scale is mapped (the value gradient)
  scales <- vapply(p$scales$scales, function(s) paste(s$aesthetics, collapse = ","),
                   character(1))
  expect_true(any(grepl("colour|color", scales)))

  expect_s3_class(densinet(nsv, color_by = "value", color_palette = "gradient"),
                  "ggplot")
})

test_that("color_trim clamps the scale (keeps trees) and drops outliers", {
  vals <- c(1:10, 1000)                       # one strong high outlier
  trees <- lapply(rep("((A,B),C);", length(vals)), parse_network)
  nsv <- anansi_netset(trees, meta = data.frame(value = vals))
  q <- stats::quantile(vals, c(0.25, 0.75), names = FALSE); iqr <- q[2] - q[1]
  hi <- min(q[2] + 3 * iqr, max(vals))

  cloud_data <- function(p) {
    L <- Filter(function(l) ".value" %in% names(l$data), p$layers)
    L[[1]]$data
  }

  # clamp: all 11 networks kept, but the outlier value is squished to the fence
  pc <- densinet(nsv, color_by = "value", color_trim = 3,
                 color_trim_action = "clamp")
  dc <- cloud_data(pc)
  expect_length(unique(dc$.net), 11L)
  expect_equal(max(dc$.value), hi)

  # drop: the outlier network is removed from the overlay
  pd <- densinet(nsv, color_by = "value", color_trim = 3,
                 color_trim_action = "drop")
  dd <- cloud_data(pd)
  expect_length(unique(dd$.net), 10L)
  expect_lt(max(dd$.value), 1000)
})

test_that("densinet draws an explicit backbone and validates its taxa", {
  trees <- lapply(c("((A,B),C);", "((A,C),B);"), parse_network)
  nsv <- anansi_netset(trees, meta = data.frame(value = c(0, 1)))
  expect_s3_class(densinet(nsv, backbone = "((A,B),C);"), "ggplot")
  expect_s3_class(densinet(nsv, color_by = "value", backbone = "((A,B),C);"),
                  "ggplot")
  expect_error(densinet(nsv, backbone = "((A,B),Z);"),
               "taxa not in the figure")
})

test_that("a reticulate backbone draws its hybrid edge as an extra layer", {
  # 4-taxon cloud + a backbone network with one reticulation over the same taxa.
  cloud <- anansi_netset(list(na, nb))
  net_bb <- "(((A:1,B:1):1,(C:1)#H1:1::0.6):1,(D:1,#H1:1::0.4):1);"
  p_tree <- densinet(cloud, method = "first", backbone = "(((A,B),C),D);")
  p_net  <- densinet(cloud, method = "first", backbone = net_bb)
  seg_layers <- function(p)
    sum(vapply(p$layers, function(L) inherits(L$geom, "GeomSegment"), logical(1)))
  # the network backbone adds one more GeomSegment layer (the reticulation edge)
  expect_gt(seg_layers(p_net), seg_layers(p_tree))
})

test_that("build_densinet threads color_by through to a ggplot", {
  trees <- lapply(c("((A,B),C);", "((A,C),B);", "((B,C),A);"), parse_network)
  nsv <- anansi_netset(trees, meta = data.frame(value = c(0, 1, 2)))
  expect_s3_class(build_densinet(nsv, list(color_by = "value")), "ggplot")
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
