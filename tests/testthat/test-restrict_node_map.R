# Tests for restrict_node_map() and node maps / offset in
# write_buildmaps_file(). The primate tests use the bundled example tree;
# maps are small synthetic grids (10 degree cells).

square <- function(xmin, xmax, ymin, ymax, crs = 4326)
  sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax),
    c(xmin, ymin)))), crs = crs))

# GPlates-like maps object for ages, with series CSVs written to a new folder
series <- function(ages) {
  cl <- rep(list(square(-10, 10, -10, 10)), length(ages))
  names(cl) <- as.character(ages)
  maps <- structure(list(type = "coastlines", coastlines = cl,
                         model = "PALEOMAP", max_age = max(ages), ages = ages,
                         n_ages = length(ages), failed = integer(0)),
                    class = "geo_paleomaps")
  d <- tempfile("series_")
  paleomap_to_csv(ages, maps = maps, resolution = 10, dir = d, verbose = FALSE)
  list(maps = maps, dir = d)
}

crown_tips <- function() {
  x <- readLines(example_files()$anc)
  strsplit(x[grepl("^Node-00008\t", x)], "\t")[[1]][-(1:4)]
}

# Small non-ultrametric tree: heights A,B,C = 0, D = 1;
# node (A,B) = 1, node (C,D) = 1.5, root = 2
small_tree <- function()
  ape::read.tree(text = "((A:1,B:1):1,(C:1.5,D:0.5):0.5);")

africa <- square(-20, 50, -35, 35)

test_that("restrict_node_map() finds crown primates and writes the node map", {
  s <- series(0:76)
  expect_message(
    m <- restrict_node_map(s$maps, map = africa, label = "AFA_64",
                           phy = example_phy(), tips = crown_tips(),
                           dir = s$dir),
    "Node 903.*age 63.909104 Ma")
  nm <- m$node_maps
  expect_equal(nrow(nm), 1)
  expect_equal(nm$node, 903L)
  expect_equal(nm$n_tips, length(crown_tips()))
  expect_equal(nm$age, 63.909104, tolerance = 1e-6)
  expect_equal(c(nm$start, nm$end), c(63.908104, 63.910104), tolerance = 1e-9)
  expect_output(print(m), "mask_AFA_64.csv : node 903")

  # same Lon/Lat points, in the same order, as the series
  node_csv <- utils::read.csv(file.path(s$dir, "mask_AFA_64.csv"))
  ser_csv  <- utils::read.csv(file.path(s$dir, "mask_0.csv"))
  expect_identical(node_csv[, 1:2], ser_csv[, 1:2])
  expect_equal(node_csv$mask[node_csv$Lon == 25 & node_csv$Lat == 25], 1L)
  expect_equal(node_csv$mask[node_csv$Lon == -105 & node_csv$Lat == 25], 0L)

  # the same node by ape number
  m2 <- restrict_node_map(s$maps, map = africa, label = "AFA_b",
                          phy = example_phy(), node = 903, dir = s$dir,
                          verbose = FALSE)
  expect_equal(m2$node_maps$age, nm$age)
})

test_that("write_buildmaps_file() splits the series around the node map", {
  s <- series(0:76)
  m <- restrict_node_map(s$maps, map = africa, label = "AFA_64",
                         phy = example_phy(), tips = crown_tips(),
                         dir = s$dir, verbose = FALSE)
  f <- tempfile(fileext = ".txt")
  write_buildmaps_file(m, out_file = f, max_age = 76.05, verbose = FALSE)
  out <- readLines(f)
  i <- grep("AFA_64", out)
  expect_equal(out[(i - 1):(i + 1)],
               c("mask_64.csv\t63.5\t63.908104",
                 "mask_AFA_64.csv\t63.908104\t63.910104",
                 "mask_64.csv\t63.910104\t64.5"))
  expect_equal(out[2], "mask_0.csv\t0.0\t0.5")
  expect_equal(out[length(out)], "mask_76.csv\t75.5\t77.5")

  # no gaps or overlaps
  tab <- utils::read.table(text = out[-1], sep = "\t")
  expect_equal(tab$V2[-1], tab$V3[-nrow(tab)])
})

test_that("restrict_node_map() warns when other nodes fall in the interval", {
  s <- series(60:70)
  expect_warning(
    restrict_node_map(s$maps, map = africa, label = "W", phy = example_phy(),
                      tips = crown_tips(), width = 0.01, dir = s$dir,
                      verbose = FALSE),
    "1 other internal node.*63.900068")
  expect_no_warning(
    restrict_node_map(s$maps, map = africa, label = "N", phy = example_phy(),
                      tips = crown_tips(), dir = s$dir, verbose = FALSE))
})

test_that("two node maps are both inserted, and overlaps stop", {
  s   <- series(0:3)
  phy <- small_tree()
  m <- restrict_node_map(s$maps, africa, "AB", phy, tips = c("A", "B"),
                         dir = s$dir, verbose = FALSE)
  m <- restrict_node_map(m, africa, "CD", phy, tips = c("C", "D"),
                         dir = s$dir, verbose = FALSE)
  f <- tempfile(fileext = ".txt")
  write_buildmaps_file(m, out_file = f, verbose = FALSE)
  expect_equal(readLines(f)[-1], c(
    "mask_0.csv\t0.0\t0.5",
    "mask_1.csv\t0.5\t0.999",
    "mask_AB.csv\t0.999\t1.001",
    "mask_1.csv\t1.001\t1.499",
    "mask_CD.csv\t1.499\t1.501",
    "mask_2.csv\t1.501\t2.5",
    "mask_3.csv\t2.5\t3.5"))

  # a node map that overlaps another one
  expect_warning(
    w <- restrict_node_map(m, africa, "wide", phy, tips = c("C", "D"),
                           width = 0.6, dir = s$dir, verbose = FALSE),
    "other internal node")
  expect_error(write_buildmaps_file(w, out_file = f, verbose = FALSE), "overlap")

  # a node older than the series
  s2 <- series(0)
  paleomap_to_csv(0, maps = s2$maps, resolution = 10, dir = s2$dir,
                  write = TRUE, verbose = FALSE)
  r  <- restrict_node_map(s2$maps, africa, "root", phy, node = 5,
                          dir = s2$dir, verbose = FALSE)
  expect_error(suppressWarnings(
    write_buildmaps_file(r, out_file = f, verbose = FALSE)), "max_age")
  write_buildmaps_file(r, out_file = f, max_age = 2, verbose = FALSE)
  expect_equal(readLines(f)[-1], c("mask_0.csv\t0.0\t1.999",
                                   "mask_root.csv\t1.999\t2.001",
                                   "mask_0.csv\t2.001\t2.5"))
})

test_that("raster node maps use mask_range, invert and the series grid", {
  s   <- series(0:3)
  phy <- small_tree()
  # 1/0 raster on a finer grid than the series: 1 in the northern hemisphere.
  # Resampling gives 0.875 / 0.125 next to the equator, hence c(0.5, 1).
  r <- terra::rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
                   resolution = 5, crs = "EPSG:4326")
  terra::values(r) <- as.integer(terra::yFromCell(r, seq_len(terra::ncell(r))) > 0)
  m <- restrict_node_map(s$maps, r, "north", phy, tips = c("A", "B"),
                         mask_range = c(0.5, 1), dir = s$dir, verbose = FALSE)
  csv <- utils::read.csv(file.path(s$dir, "mask_north.csv"))
  expect_identical(csv[, 1:2],
                   utils::read.csv(file.path(s$dir, "mask_0.csv"))[, 1:2])
  expect_true(all(csv$mask[csv$Lat > 0] == 1L))
  expect_true(all(csv$mask[csv$Lat < 0] == 0L))

  # invert, and NA cells become 0
  r2 <- r
  terra::values(r2)[terra::yFromCell(r2, seq_len(terra::ncell(r2))) < -40 &
                    terra::xFromCell(r2, seq_len(terra::ncell(r2))) > 0 &
                    terra::xFromCell(r2, seq_len(terra::ncell(r2))) < 10] <- NA
  restrict_node_map(s$maps, r2, "south", phy, tips = c("A", "B"),
                    mask_range = c(0.5, 1), invert = TRUE, dir = s$dir,
                    verbose = FALSE)
  csv2 <- utils::read.csv(file.path(s$dir, "mask_south.csv"))
  expect_false(anyNA(csv2$mask))
  expect_equal(csv2$mask[csv2$Lon == 5 & csv2$Lat == -45], 0L)
  expect_equal(csv2$mask[csv2$Lon == 5 & csv2$Lat == -25], 1L)

  # raster file with two layers needs values_var
  f2 <- tempfile(fileext = ".tif")
  two <- c(r, 1L - r); names(two) <- c("north", "south")
  terra::writeRaster(two, f2)
  expect_error(restrict_node_map(s$maps, f2, "x", phy, tips = c("A", "B"),
                                 mask_range = c(1, 1), dir = s$dir),
               "layers")
  restrict_node_map(s$maps, f2, "x", phy, tips = c("A", "B"),
                    mask_range = c(0.5, 1), values_var = "south", dir = s$dir,
                    verbose = FALSE)
  csv3 <- utils::read.csv(file.path(s$dir, "mask_x.csv"))
  expect_true(all(csv3$mask[csv3$Lat < 0] == 1L))

  expect_error(restrict_node_map(s$maps, r, "y", phy, tips = c("A", "B"),
                                 dir = s$dir), "mask_range")
})

test_that("vector node maps can be files, and invert gives the sea", {
  s   <- series(0:3)
  phy <- small_tree()
  g <- tempfile(fileext = ".gpkg")
  sf::st_write(africa, g, quiet = TRUE)
  restrict_node_map(s$maps, g, "afr_file", phy, tips = c("A", "B"),
                    invert = TRUE, dir = s$dir, verbose = FALSE)
  csv <- utils::read.csv(file.path(s$dir, "mask_afr_file.csv"))
  expect_equal(csv$mask[csv$Lon == 25 & csv$Lat == 25], 0L)
  expect_equal(csv$mask[csv$Lon == -105 & csv$Lat == 25], 1L)

  merc <- sf::st_transform(africa, 3857)
  expect_error(restrict_node_map(s$maps, merc, "m", phy, tips = c("A", "B"),
                                 dir = s$dir), "EPSG:4326")
  pts <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326))
  expect_error(restrict_node_map(s$maps, pts, "p", phy, tips = c("A", "B"),
                                 dir = s$dir), "no polygons")
})

test_that("restrict_node_map() checks its inputs", {
  s   <- series(0:3)
  phy <- small_tree()
  expect_error(restrict_node_map(s$maps, africa, "64", phy, tips = c("A", "B"),
                                 dir = s$dir), "must not be a number")
  expect_error(restrict_node_map(s$maps, africa, "a b", phy, tips = c("A", "B"),
                                 dir = s$dir), "without spaces")
  expect_error(restrict_node_map(s$maps, africa, "x", phy, tips = c("A", "Z"),
                                 dir = s$dir), "Z")
  expect_error(restrict_node_map(s$maps, africa, "x", phy, tips = c("A", "B"),
                                 node = 6, dir = s$dir), "either")
  expect_error(restrict_node_map(s$maps, africa, "x", phy, node = 2,
                                 dir = s$dir), "internal")
  expect_error(restrict_node_map(s$maps, africa, "x", phy, tips = c("A", "B"),
                                 dir = tempfile()), "paleomap_to_csv")
  expect_error(restrict_node_map(s$maps, africa, "x", phy, tips = c("A", "B"),
                                 width = 0, dir = s$dir), "width")
  m <- restrict_node_map(s$maps, africa, "x", phy, tips = c("A", "B"),
                         dir = s$dir, verbose = FALSE)
  expect_error(restrict_node_map(m, africa, "x", phy, tips = c("C", "D"),
                                 dir = s$dir), "already")
})

test_that("offset shifts the series and node maps into tree time", {
  s <- series(64:70)
  f <- tempfile(fileext = ".txt")
  write_buildmaps_file(s$maps, out_file = f, offset = 66, verbose = FALSE)
  expect_equal(readLines(f)[-1], c("mask_66.csv\t0.0\t0.5",
                                   "mask_67.csv\t0.5\t1.5",
                                   "mask_68.csv\t1.5\t2.5",
                                   "mask_69.csv\t2.5\t3.5",
                                   "mask_70.csv\t3.5\t4.5"))
  expect_error(write_buildmaps_file(s$maps, out_file = f, offset = 80,
                                    verbose = FALSE), "younger than offset")

  # node (C,D) has height 1.5 -> real age 67.5 with offset 66
  phy <- small_tree()
  expect_message(
    m <- restrict_node_map(s$maps, africa, "CD", phy, tips = c("C", "D"),
                           offset = 66, dir = s$dir),
    "age 67.5 Ma.*tree time 1.499 to 1.501")
  write_buildmaps_file(m, out_file = f, offset = 66, verbose = FALSE)
  expect_equal(readLines(f)[3:5], c("mask_67.csv\t0.5\t1.499",
                                    "mask_CD.csv\t1.499\t1.501",
                                    "mask_68.csv\t1.501\t2.5"))
  expect_error(write_buildmaps_file(m, out_file = f, verbose = FALSE),
               "same offset")
  expect_error(restrict_node_map(m, africa, "AB", phy, tips = c("A", "B"),
                                 dir = s$dir), "same offset")
})

test_that(".fmt_time() writes up to 6 decimals", {
  expect_equal(.fmt_time(c(0, 0.5, 63.9081, 63.9093, 77, 1e-6)),
               c("0.0", "0.5", "63.9081", "63.9093", "77.0", "0.000001"))
})
