test_that("branch distances cover every edge and iteration", {
  geo <- example_geo()
  phy <- example_phy()
  bd  <- get_dispersal_distance(geo, phy, type = "branch")

  expect_s3_class(bd, "data.frame")
  expect_equal(nrow(bd), nrow(phy$edge) * geo$n_iterations)
  expect_false(anyNA(bd$distance_km))
  expect_true(all(bd$distance_km >= 0))
  expect_identical(attr(bd, "frame"), "palaeo")
  expect_null(attr(bd, "reconstruction"))

  # distance matches an independent Haversine calculation
  r <- bd[123, ]
  expect_equal(r$distance_km,
               geosphere::distHaversine(c(r$from_lon, r$from_lat),
                                        c(r$to_lon, r$to_lat)) / 1000)

  # fixed branch time comes from the tree
  e <- which(phy$edge[, 1] == r$parent_ape_node & phy$edge[, 2] == r$child_ape_node)
  expect_equal(r$branch_time_myr, phy$edge.length[e])
})

test_that("nodes = 'tips' keeps only terminal branches", {
  geo <- example_geo()
  phy <- example_phy()
  bd  <- get_dispersal_distance(geo, phy, type = "branch", nodes = "tips")
  expect_true(all(bd$child_ape_node <= ape::Ntip(phy)))
  expect_equal(nrow(bd), ape::Ntip(phy) * geo$n_iterations)
  expect_equal(bd$tip_label, phy$tip.label[bd$child_ape_node])
})

test_that("path distance equals the sum of branch distances along the path", {
  geo <- example_geo()
  phy <- example_phy()
  both <- get_dispersal_distance(geo, phy, type = "both", nodes = "all")
  bd <- both$branch
  pd <- both$path
  it <- bd$iteration[1]

  tip  <- match("Lemur_catta", phy$tip.label)
  path <- ape::nodepath(phy, ape::Ntip(phy) + 1, tip)
  sub  <- bd[bd$iteration == it, ]
  steps <- vapply(seq_len(length(path) - 1), function(k)
    sub$distance_km[sub$parent_ape_node == path[k] &
                      sub$child_ape_node == path[k + 1]], numeric(1))

  got <- pd$path_distance_km[pd$iteration == it & pd$ape_node == tip]
  expect_equal(got, sum(steps))

  root <- pd[pd$iteration == it & pd$ape_node == ape::Ntip(phy) + 1, ]
  expect_equal(root$path_distance_km, 0)
  expect_equal(root$lineage_time_myr, 0)
})

test_that("missing ape_node mapping gives a clear error", {
  geo <- read_geo_ancstates(example_files()$anc, verbose = FALSE)
  expect_error(get_dispersal_distance(geo, example_phy()), "ape_node mapping")
})

test_that("frame = 'plate' uses reconstructed coordinates and reports failures", {
  geo <- example_geo()
  phy <- example_phy()

  # Fake reconstruction: shift longitudes by 1 degree per Ma, and fail all
  # points at 50 Ma plus the first point at 10 Ma.
  fake <- function(xy, age, model, chunk = 2000L) {
    out <- xy
    out[, 1] <- ((xy[, 1] + age + 180) %% 360) - 180
    err <- NULL
    if (age == 50) { out[] <- NA; err <- "simulated failure" }
    if (age == 10) out[1, ] <- NA
    list(xy = out, error = err)
  }
  local_mocked_bindings(.reconstruct_points = fake)

  expect_warning(
    bd <- get_dispersal_distance(geo, phy, frame = "plate", verbose = FALSE),
    "could not be reconstructed")

  expect_identical(attr(bd, "frame"), "plate")
  rec <- attr(bd, "reconstruction")
  expect_s3_class(rec, "data.frame")
  expect_equal(sum(rec$n_samples), nrow(geo$coords))
  at50 <- rec$reconstruction_age == 50
  expect_true(all(rec$n_failed[at50] == rec$n_samples[at50]))
  expect_equal(sum(rec$n_failed[rec$reconstruction_age == 10]), 1)
  expect_true(all(rec$n_failed[!rec$reconstruction_age %in% c(10, 50)] == 0))
  expect_identical(attr(rec, "errors")[["50"]], "simulated failure")

  # failed samples are kept (same number of rows), with NA distances
  expect_equal(nrow(bd), nrow(phy$edge) * geo$n_iterations)
  expect_true(anyNA(bd$distance_km))

  # an edge between two age-0 tips' parents is unaffected only if both ends
  # are at age 0; check one fully reconstructed edge instead
  ok <- bd[!is.na(bd$distance_km), ][1, ]
  expect_equal(ok$distance_km,
               geosphere::distHaversine(c(ok$from_lon, ok$from_lat),
                                        c(ok$to_lon, ok$to_lat)) / 1000)
})

test_that("frame = 'plate' with no failures gives no warning", {
  local_mocked_bindings(.reconstruct_points = function(xy, age, model, chunk = 2000L)
    list(xy = xy, error = NULL))
  geo <- example_geo()
  phy <- example_phy()
  expect_no_warning(
    bd <- get_dispersal_distance(geo, phy, frame = "plate", verbose = FALSE))
  ref <- get_dispersal_distance(geo, phy)
  # identity reconstruction: same distances as the palaeo frame
  expect_equal(bd$distance_km, ref$distance_km)
  expect_true(all(attr(bd, "reconstruction")$n_failed == 0))
})

test_that(".reconstruct_points never errors and returns NA on failure", {
  local_mocked_bindings(reconstruct = function(...) stop("no connection"),
                        .package = "rgplates")
  res <- GeoModelR:::.reconstruct_points(cbind(c(1, 2, 3), c(4, 5, 6)), age = 10,
                                         model = "PALEOMAP")
  expect_true(all(is.na(res$xy)))
  expect_match(res$error, "no connection")
})
