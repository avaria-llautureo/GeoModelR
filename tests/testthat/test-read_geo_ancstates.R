test_that("example AncStates file is parsed correctly", {
  geo <- read_geo_ancstates(example_files()$anc, verbose = FALSE)

  expect_s3_class(geo, "geo_ancstates")
  expect_equal(geo$n_nodes, 1783L)
  expect_equal(geo$n_iterations, 50L)
  expect_equal(sum(geo$nodes$n_taxa == 1), 894L)
  expect_equal(nrow(geo$coords), 1783L * 50L)
  expect_true(all(c("iteration", "node_name", "longitude", "latitude",
                    "branch_length_scaled", "height", "is_tip") %in%
                    names(geo$coords)))
  expect_true(all(abs(geo$coords$longitude) <= 180))
  expect_true(all(abs(geo$coords$latitude) <= 90))
  expect_equal(max(geo$nodes$height), 74.046823, tolerance = 1e-6)
})

test_that("print and summary methods run", {
  geo <- example_geo()
  expect_output(print(geo), "BayesTraits Geo Model")
  expect_output(summary(geo), "geo_ancstates summary")
})

test_that("a file without the expected blocks gives a clear error", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("not", "an AncStates file"), f)
  expect_error(read_geo_ancstates(f, verbose = FALSE), "block headers")
})

test_that("utility helpers work on the example data", {
  geo <- example_geo()
  ages <- get_node_ages(geo)
  expect_equal(ages, sort(unique(round(geo$nodes$height))))
  expect_true("Lemur_catta" %in% find_tips(geo, "^Lemur_"))
  expect_length(find_tips(geo, "no_such_taxon"), 0)
})
