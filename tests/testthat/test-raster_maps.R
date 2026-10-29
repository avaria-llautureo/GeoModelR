# Tests for raster palaeomaps (load_raster_maps(), paleomap_to_csv() with
# source = "raster", export_restriction_maps(), write_buildmaps_file(),
# fetch_etopo2022()). All rasters are small synthetic GeoTIFFs written here.

# One synthetic "elevation" raster: value = latitude * 100 + shift, so cells
# north of the equator are positive (land) and south negative (sea).
make_rast <- function(res = 10, shift = 0, name = "z", crs = "EPSG:4326") {
  r <- terra::rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                   resolution = res, crs = crs)
  terra::values(r) <- terra::yFromCell(r, seq_len(terra::ncell(r))) * 100 + shift
  names(r) <- name
  r
}

# A folder with one Scotese-style file per age: Map_<age>Ma.tif
make_series <- function(ages, res = 10) {
  d <- tempfile("series_"); dir.create(d)
  for (a in ages)
    terra::writeRaster(make_rast(res, shift = a),
                       file.path(d, sprintf("Map_%gMa.tif", a)))
  d
}

load_series <- function(d, ...)
  load_raster_maps(d, age_regex = "_([0-9]+)Ma\\.tif$", pattern = "\\.tif$",
                   verbose = FALSE, ...)

test_that("load_raster_maps() finds ages with age_regex, present to past", {
  d <- make_series(c(10, 0, 5))
  m <- load_series(d)
  expect_s3_class(m, "geo_paleomaps")
  expect_identical(m$type, "raster")
  expect_equal(m$ages, c(0, 5, 10))
  expect_equal(m$files$age, c(0, 5, 10))
  expect_equal(m$values_var, "z")
  expect_equal(m$native_res, c(10, 10))
  expect_output(print(m), "local raster")

  # 'ages' selects a time span and warns about ages without a file
  expect_warning(m2 <- load_series(d, ages = c(0, 5, 7)), "7")
  expect_equal(m2$ages, c(0, 5))

  expect_error(load_raster_maps(d, age_regex = "^nomatch([0-9]+)",
                                pattern = "\\.tif$"), "matched")
  expect_error(load_raster_maps(age_regex = "x"), "dir")
  expect_error(load_raster_maps(dir = tempfile()), "dir")
  expect_error(load_raster_maps(dir = tempfile(), age_regex = "x"), "not found")
})

test_that("load_raster_maps() stops on duplicate ages", {
  d <- make_series(c(0, 5))
  file.copy(file.path(d, "Map_5Ma.tif"), file.path(d, "Other_5Ma.tif"))
  expect_error(load_series(d), "same age")
})

test_that("load_raster_maps() checks the CRS and the variable", {
  d <- tempfile(); dir.create(d)
  r <- terra::rast(xmin = 0, xmax = 1e6, ymin = 0, ymax = 1e6, resolution = 1e5,
                   crs = "EPSG:3857")
  terra::values(r) <- 1
  terra::writeRaster(r, file.path(d, "Map_0Ma.tif"))
  expect_error(load_series(d), "EPSG:4326")
  expect_s3_class(load_series(d, check_crs = FALSE), "geo_paleomaps")

  # Two layers: values_var is required, and selects the layer
  d2 <- tempfile(); dir.create(d2)
  two <- c(make_rast(name = "elev"), make_rast(shift = 10000, name = "temp"))
  terra::writeRaster(two, file.path(d2, "Map_0Ma.tif"))
  expect_error(load_series(d2), "values_var")
  m <- load_series(d2, values_var = "temp")
  expect_equal(m$files$values_var, "temp")
  df <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(0, Inf))
  expect_true(all(df$mask == 1L))          # temp layer is positive everywhere
  expect_error(load_series(d2, values_var = "nope"), "nope")
})

test_that("mask_range bounds are inclusive and several ranges combine", {
  m  <- load_series(make_series(0))
  # values at cell centres are -8500, -7500, ..., 8500
  df <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(500, 1500))
  expect_equal(sort(unique(df$Lat[df$mask == 1L])), c(5, 15))

  df <- paleomap_to_csv(0, maps = m, source = "raster",
                        mask_range = list(c(-1500, -500), c(7500, Inf)))
  expect_equal(sort(unique(df$Lat[df$mask == 1L])), c(-15, -5, 75, 85))
  expect_true(all(df$mask %in% c(0L, 1L)))
  expect_equal(nrow(df), 36 * 18)

  # invert swaps 1 and 0
  di <- paleomap_to_csv(0, maps = m, source = "raster",
                        mask_range = list(c(-1500, -500), c(7500, Inf)),
                        invert = TRUE)
  expect_identical(di$mask, 1L - df$mask)

  expect_error(paleomap_to_csv(0, maps = m, source = "raster"), "mask_range")
  expect_error(paleomap_to_csv(0, maps = m, source = "raster",
                               mask_range = c(5, 1)), "min <= max")
  expect_error(paleomap_to_csv(0, maps = m, source = "raster",
                               mask_range = "land"), "numeric")
  expect_error(paleomap_to_csv(0, maps = m, mask_range = c(0, 1)),
               "load_raster_maps")
  expect_error(paleomap_to_csv(0, source = "raster", mask_range = c(0, 1)),
               "load_raster_maps")
})

test_that("cells with no value (NA) always become 0", {
  d <- tempfile(); dir.create(d)
  r <- make_rast()
  r[terra::cellFromXY(r, cbind(5, 45))] <- NA      # a land cell with no value
  terra::writeRaster(r, file.path(d, "Map_0Ma.tif"))
  m <- load_series(d)
  for (inv in c(FALSE, TRUE)) {
    df <- paleomap_to_csv(0, maps = m, source = "raster",
                          mask_range = c(0, Inf), invert = inv)
    expect_false(anyNA(df$mask))
    expect_equal(df$mask[df$Lon == 5 & df$Lat == 45], 0L)
  }
})

test_that("resolution resamples every age, with either method", {
  m  <- load_series(make_series(c(0, 5)))
  df <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(0, Inf),
                        resolution = 5)
  expect_equal(nrow(df), 72 * 36)
  expect_true(all(c(-177.5, 177.5) %in% df$Lon))
  dc <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(0, Inf),
                        resolution = 5, resample_method = "cubicspline")
  expect_equal(nrow(dc), 72 * 36)
  expect_true(all(dc$mask %in% c(0L, 1L)))

  # extent crops the output
  de <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(0, Inf),
                        extent = c(-20, 20, -10, 10))
  expect_true(all(abs(de$Lon) < 20 & abs(de$Lat) < 10))

  o <- tempfile()
  paleomap_to_csv(m$ages, maps = m, source = "raster", mask_range = c(0, Inf),
                  resolution = 5, dir = o, verbose = FALSE)
  expect_setequal(list.files(o), c("mask_0.csv", "mask_5.csv"))
  expect_equal(nrow(utils::read.csv(file.path(o, "mask_5.csv"))), 72 * 36)
})

test_that("age0_override replaces or adds the 0 Ma map", {
  d  <- make_series(c(0, 5))
  ov <- tempfile(fileext = ".tif")
  terra::writeRaster(make_rast(shift = 100000), ov)   # positive everywhere

  m <- suppressMessages(load_series(d, age0_override = ov))
  expect_equal(m$ages, c(0, 5))
  expect_equal(m$files$path[m$files$age == 0], ov)
  df0 <- paleomap_to_csv(0, maps = m, source = "raster", mask_range = c(0, Inf))
  expect_true(all(df0$mask == 1L))
  df5 <- paleomap_to_csv(5, maps = m, source = "raster", mask_range = c(0, Inf))
  expect_false(all(df5$mask == 1L))

  # series without a 0 Ma file: override is added as 0 Ma
  m2 <- load_series(make_series(c(5, 10)),
                    age0_override = list(path = ov, values_var = "z"))
  expect_equal(m2$ages, c(0, 5, 10))
  expect_equal(m2$files$path[1], ov)

  expect_error(load_series(d, age0_override = tempfile()), "existing file")
})

test_that("the size warning checks each map, including age0_override", {
  d  <- make_series(c(0, 5), res = 10)
  ov <- tempfile(fileext = ".tif")
  terra::writeRaster(make_rast(res = 0.25), ov)
  m  <- load_series(d, age0_override = ov)

  # coarse palaeo series alone: no warning
  expect_no_warning(paleomap_to_csv(5, maps = m, source = "raster",
                                    mask_range = c(0, Inf)))
  # fine present-day layer at native resolution: warning names 0 Ma
  expect_warning(paleomap_to_csv(c(0, 5), maps = m, source = "raster",
                                 mask_range = c(0, Inf), dir = tempfile(),
                                 verbose = FALSE),
                 "0 Ma \\(0.25 x 0.25 deg")
  # resolution set: no warning
  expect_no_warning(paleomap_to_csv(0, maps = m, source = "raster",
                                    mask_range = c(0, Inf), resolution = 1))
})

test_that("load_raster_maps(file = ) registers one present-day map", {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(make_rast(), f)
  p <- load_raster_maps(file = f, verbose = FALSE, label = "Present")
  expect_equal(p$ages, 0)
  expect_equal(p$files$path, f)
  expect_equal(p$model, "Present")
  df <- paleomap_to_csv(0, maps = p, source = "raster", mask_range = c(0, Inf))
  expect_equal(nrow(df), 36 * 18)

  expect_error(load_raster_maps(file = f, dir = tempdir(), age_regex = "x"),
               "not both")
  expect_error(load_raster_maps(file = f, ages = 0), "ages")
  expect_error(load_raster_maps(file = f, age0_override = f), "age0_override")
  expect_error(load_raster_maps(file = tempfile()), "not found")
})

test_that("write_buildmaps_file() handles 5 Ma spacing and max_age", {
  m <- load_series(make_series(c(0, 5, 10)))
  f <- tempfile(fileext = ".txt")
  write_buildmaps_file(m, out_file = f, verbose = FALSE)
  expect_equal(readLines(f)[-1], c("mask_0.csv\t0.0\t2.5",
                                   "mask_5.csv\t2.5\t7.5",
                                   "mask_10.csv\t7.5\t10.5"))

  # max_age extends the last map to cover the root
  write_buildmaps_file(m, out_file = f, max_age = 12.3, verbose = FALSE)
  expect_equal(readLines(f)[4], "mask_10.csv\t7.5\t13.5")
  expect_error(write_buildmaps_file(m, out_file = f, max_age = 8), "younger")
  expect_error(write_buildmaps_file(m, out_file = f, max_age = "a"), "max_age")

  # ages selects a subset
  write_buildmaps_file(m, out_file = f, ages = c(0, 5), verbose = FALSE)
  expect_equal(readLines(f)[-1], c("mask_0.csv\t0.0\t2.5", "mask_5.csv\t2.5\t5.5"))
})

test_that("a single present-day map can cover the whole tree", {
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(make_rast(), f)
  p   <- load_raster_maps(file = f, verbose = FALSE)
  out <- tempfile(fileext = ".txt")

  write_buildmaps_file(p, out_file = out, max_age = 1.03, verbose = FALSE)
  expect_equal(readLines(out)[-1], "mask_0.csv\t0.0\t2.5")

  expect_warning(write_buildmaps_file(p, out_file = out, verbose = FALSE),
                 "Only one map")

  # the same with a GPlates present-day map
  sq <- sf::st_sf(id = 1, geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(-10, -10), c(10, -10), c(10, 10),
                              c(-10, 10), c(-10, -10)))), crs = 4326))
  gp <- structure(list(type = "coastlines", coastlines = list(`0` = sq),
                       model = "PALEOMAP", max_age = 0L, ages = 0L, n_ages = 1L,
                       failed = integer(0)), class = "geo_paleomaps")
  write_buildmaps_file(gp, out_file = out, max_age = 0.8, verbose = FALSE)
  expect_equal(readLines(out)[-1], "mask_0.csv\t0.0\t1.5")
})

test_that("export_restriction_maps() traces the mask boundary for rasters", {
  m <- load_series(make_series(c(0, 5)))
  o <- tempfile()
  export_restriction_maps(m, dir = o, mask_range = c(0, Inf), verbose = FALSE)
  expect_setequal(list.files(o), c("mask_0.csv", "mask_5.csv"))
  xy <- utils::read.csv(file.path(o, "mask_0.csv"), header = FALSE)
  expect_equal(ncol(xy), 2)
  expect_true(all(abs(xy[[2]]) < 1e-8))     # boundary is the equator

  o2 <- tempfile()
  export_restriction_maps(m, dir = o2, ages = 5, mask_range = c(0, Inf),
                          verbose = FALSE)
  expect_equal(list.files(o2), "mask_5.csv")

  # a range no cell falls in: no boundary, no file
  o3 <- tempfile()
  export_restriction_maps(m, dir = o3, mask_range = c(1e6, Inf), verbose = FALSE)
  expect_length(list.files(o3), 0)

  expect_error(export_restriction_maps(m, dir = tempfile()), "mask_range")
})

test_that("raster maps are refused where GPlates maps are needed", {
  m <- load_series(make_series(0))
  expect_error(fill_missing_maps(m), "load_raster_maps")
  expect_error(.get_palaeomap(0, maps = m), "get_paleomaps")
})

test_that("fetch_etopo2022() uses the cache and the package's own name", {
  dest <- eval(formals(fetch_etopo2022)$dest, environment(fetch_etopo2022))
  expect_match(dest, utils::packageName(environment(fetch_etopo2022)),
               fixed = TRUE)

  f <- tempfile(fileext = ".nc")
  writeLines("cached", f)
  expect_message(out <- fetch_etopo2022(dest = f), "cached")
  expect_identical(out, f)
})
