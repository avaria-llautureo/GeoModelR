fake_maps <- function(ages, failed = integer(0)) {
  sq <- sf::st_sf(id = 1, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(-10, -10), c(10, -10), c(10, 10),
                              c(-10, 10), c(-10, -10)))), crs = 4326))
  cl <- lapply(ages, function(a) if (a %in% failed) NULL else sq)
  names(cl) <- as.character(ages)
  structure(list(coastlines = cl, model = "PALEOMAP", max_age = max(ages),
                 ages = ages, n_ages = length(ages) - length(failed),
                 failed = failed), class = "geo_paleomaps")
}

test_that("write_buildmaps_file() writes contiguous time intervals", {
  f <- tempfile(fileext = ".txt")
  write_buildmaps_file(fake_maps(0:3), out_file = f, verbose = FALSE)
  out <- readLines(f)
  expect_equal(out[1], "BinaryMap.bin\t250\t150")
  expect_equal(out[-1], c("mask_0.csv\t0.0\t0.5", "mask_1.csv\t0.5\t1.5",
                          "mask_2.csv\t1.5\t2.5", "mask_3.csv\t2.5\t3.5"))
})

test_that("write_buildmaps_file() skips failed ages with a warning", {
  f <- tempfile(fileext = ".txt")
  expect_warning(write_buildmaps_file(fake_maps(0:3, failed = 2L), out_file = f,
                                      verbose = FALSE), "skipped")
  out <- readLines(f)
  expect_false(any(grepl("mask_2.csv", out)))
  expect_equal(out[3], "mask_1.csv\t0.5\t2.0")
  expect_error(write_buildmaps_file(fake_maps(0:1), bin_file = "a b.bin",
                                    out_file = f), "spaces")
})

test_that("get_paleomaps() keeps a NULL slot for failed ages", {
  sq <- fake_maps(0)$coastlines[[1]]
  local_mocked_bindings(
    reconstruct = function(x, age, model, ...) {
      if (age == 1) stop("simulated GWS failure")
      sq
    },
    .package = "rgplates")
  expect_warning(maps <- suppressMessages(get_paleomaps(ages = 0:3)), "failed")
  expect_identical(names(maps$coastlines), c("0", "1", "2", "3"))
  expect_null(maps$coastlines[["1"]])
  expect_false(is.null(maps$coastlines[["3"]]))
  expect_identical(maps$failed, 1L)
  expect_output(print(maps), "Failed")

  # fill_missing_maps() retries the failed age
  local_mocked_bindings(reconstruct = function(x, age, model, ...) sq,
                        .package = "rgplates")
  maps2 <- suppressMessages(fill_missing_maps(maps))
  expect_length(maps2$failed, 0)
  expect_false(is.null(maps2$coastlines[["1"]]))
})

test_that("sf_to_csv() rasterises a maps object to a binary grid", {
  df <- sf_to_csv(age = 0, maps = fake_maps(0), resolution = 5, verbose = FALSE)
  expect_named(df, c("Lon", "Lat", "mask"))
  expect_setequal(unique(df$mask), c(0L, 1L))
  expect_equal(df$mask[df$Lon == 2.5 & df$Lat == 2.5], 1L)
  expect_equal(df$mask[df$Lon == 102.5 & df$Lat == 42.5], 0L)

  d <- file.path(tempdir(), "masks")
  sf_to_csv(age = 0:1, maps = fake_maps(0:1), resolution = 10, dir = d,
            verbose = FALSE)
  expect_true(all(file.exists(file.path(d, c("mask_0.csv", "mask_1.csv")))))
})
