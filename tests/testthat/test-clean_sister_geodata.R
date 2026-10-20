test_that("pairs include cherries, triples and polytomies", {
  pairs <- GeoModelR:::.get_sister_pairs(toy_tree(), include_triples = TRUE)
  key <- vapply(pairs, function(p) paste(sort(p), collapse = "-"), "")
  expect_setequal(key, c("A-B", "D-E", "C-D", "C-E"))

  pairs2 <- GeoModelR:::.get_sister_pairs(toy_tree(), include_triples = FALSE)
  expect_length(pairs2, 2)

  poly <- ape::read.tree(text = "((A,B,C),D);")
  key3 <- vapply(GeoModelR:::.get_sister_pairs(poly), function(p)
    paste(sort(p), collapse = "-"), "")
  expect_setequal(key3, c("A-B", "A-C", "B-C"))
})

test_that("no overlaps: data returned unchanged", {
  res <- clean_sister_geodata(toy_tree(), toy_single(), verbose = FALSE)
  expect_s3_class(res, "geo_cleaned")
  expect_equal(res$clean_singledata, toy_single())
  expect_equal(nrow(res$deletion_log), 0)
  expect_equal(nrow(res$clean_distdata), 0)
  expect_identical(res$clean_tree, toy_tree())
  expect_output(print(res), "Sister-taxon coordinate cleaning")
})

test_that("overlapping centroids stop with a list of pairs", {
  single <- toy_single()
  single$V2[2] <- -70.0002          # B about 19 m from A
  err <- tryCatch(clean_sister_geodata(toy_tree(), single, verbose = FALSE),
                  geo_centroid_conflict = function(e) e)
  expect_s3_class(err, "geo_centroid_conflict")
  expect_equal(err$conflicts$taxon_1, "A")
  expect_equal(err$conflicts$taxon_2, "B")
  expect_lt(err$conflicts$distance_km, 0.05)
  expect_match(conditionMessage(err), "delete = c\\(\"A\"\\)")

  # resolved by the user: B removed from data and tree
  res <- clean_sister_geodata(toy_tree(), single, delete = "B", verbose = FALSE)
  expect_false("B" %in% res$clean_singledata$V1)
  expect_false("B" %in% res$clean_tree$tip.label)
  expect_identical(res$deleted_taxa, "B")
  expect_equal(single[single$V1 == "A", "V2"],
               res$clean_singledata[res$clean_singledata$V1 == "A", "V2"])
})

test_that("threshold controls what counts as overlap (default 50 m)", {
  single <- toy_single()
  single$V2[2] <- -70.0007          # about 67 m from A at 30 S
  expect_silent(clean_sister_geodata(toy_tree(), single, verbose = FALSE))
  expect_error(clean_sister_geodata(toy_tree(), single,
                                    distance_threshold = 0.1, verbose = FALSE),
               class = "geo_centroid_conflict")
})

test_that("both distribution sets: points removed only from the larger set", {
  dist <- data.frame(V1 = c("D", "E"), V2 = "Linked",
                     V3 = c("10,10.5,11", "10.00001,12"),
                     V4 = c("45,45.5,46", "45,47"))
  res <- clean_sister_geodata(toy_tree(), toy_single(), dist, verbose = FALSE)

  log <- res$deletion_log
  expect_equal(nrow(log), 1)
  expect_identical(log$deleted_from, "D")
  expect_equal(c(log$lon, log$lat), c(10, 45))
  expect_match(log$reason, "more points")

  expect_equal(res$clean_distdata$V3[res$clean_distdata$V1 == "D"], "10.5,11")
  expect_equal(res$clean_distdata$V3[res$clean_distdata$V1 == "E"], "10.00001,12")
  expect_length(res$fallback_to_centroid, 0)
})

test_that("equal-sized sets: random choice, reproducible with seed, RNG restored", {
  dist <- data.frame(V1 = c("D", "E"), V2 = "Linked",
                     V3 = c("10,10.5,11", "10.00001,12,13"),
                     V4 = c("45,45.5,46", "45,47,48"))
  set.seed(123); before <- runif(1)
  set.seed(123)
  a <- clean_sister_geodata(toy_tree(), toy_single(), dist, seed = 7, verbose = FALSE)
  after <- runif(1)
  expect_identical(before, after)

  b <- clean_sister_geodata(toy_tree(), toy_single(), dist, seed = 7, verbose = FALSE)
  expect_identical(a$deletion_log, b$deletion_log)
  expect_equal(nrow(a$deletion_log), 1)
  expect_match(a$deletion_log$reason, "random")

  # across seeds, both taxa get chosen
  losers <- vapply(1:30, function(s)
    clean_sister_geodata(toy_tree(), toy_single(), dist, seed = s,
                         verbose = FALSE)$deletion_log$deleted_from, "")
  expect_setequal(unique(losers), c("D", "E"))
})

test_that("distribution points overlapping a sister centroid are removed", {
  # C is sister to (D, E); C uses its centroid, D has distribution data
  single <- toy_single()
  single[single$V1 == "C", c("V2", "V3")] <- c(20, 50)
  dist <- data.frame(V1 = "D", V2 = "Linked",
                     V3 = "20.00001,30,31", V4 = "50,50,51")
  res <- clean_sister_geodata(toy_tree(), single, dist, verbose = FALSE)
  expect_identical(res$deletion_log$deleted_from, "D")
  expect_equal(res$clean_distdata$V3, "30,31")
  expect_equal(res$clean_singledata, single)   # centroids untouched
})

test_that("taxa left with < 2 points fall back to their centroid", {
  dist <- data.frame(V1 = c("D", "E"), V2 = "Linked",
                     V3 = c("10,11,12", "10.00001,11.00001"),
                     V4 = c("45,46,47", "45,46"))
  res <- clean_sister_geodata(toy_tree(), toy_single(), dist, verbose = FALSE)
  expect_true("D" %in% res$fallback_to_centroid)
  expect_false("D" %in% res$clean_distdata$V1)
  expect_true("D" %in% res$clean_singledata$V1)
})

test_that("input data rules are enforced", {
  single <- toy_single()

  extra <- rbind(single, data.frame(V1 = "Z", V2 = 0, V3 = 0))
  expect_output(
    expect_message(res <- clean_sister_geodata(toy_tree(), extra), "not in the tree"),
    "Dropped \\(not in tree\\)")
  expect_identical(res$not_in_tree$singledata, "Z")

  expect_error(clean_sister_geodata(toy_tree(), single[-1, ], verbose = FALSE),
               "no centroid")
  expect_error(clean_sister_geodata(toy_tree(), rbind(single, single[1, ]),
                                    verbose = FALSE), "more than one row")

  bad <- single; bad$V3[1] <- 95
  expect_error(clean_sister_geodata(toy_tree(), bad, verbose = FALSE),
               "invalid centroid")

  dist <- data.frame(V1 = "D", V2 = "Linked", V3 = "10,11", V4 = "45")
  expect_error(clean_sister_geodata(toy_tree(), single, dist, verbose = FALSE),
               "Invalid distribution coordinates")

  expect_error(clean_sister_geodata(toy_tree(), single, delete = "Q",
                                    verbose = FALSE), "not found")
  expect_error(clean_sister_geodata(toy_tree(), single, include_distdata = TRUE,
                                    verbose = FALSE), "NULL")
})

test_that("ask = TRUE fills in 'delete' from the menu", {
  single <- toy_single()
  single$V2[2] <- -70.0002
  local_mocked_bindings(.is_interactive = function() TRUE,
                        menu = function(choices, title) 2L)   # "Delete B"
  expect_message(
    res <- clean_sister_geodata(toy_tree(), single, ask = TRUE, verbose = FALSE),
    "delete = c\\(\"B\"\\)")
  expect_identical(res$deleted_taxa, "B")

  # choosing "keep both" still stops
  local_mocked_bindings(menu = function(choices, title) 3L)
  expect_error(clean_sister_geodata(toy_tree(), single, ask = TRUE, verbose = FALSE),
               class = "geo_centroid_conflict")
})

test_that("works with multiPhylo input", {
  trees <- c(toy_tree(), toy_tree())
  single <- toy_single(); single$V2[2] <- -70.0002
  res <- clean_sister_geodata(trees, single, delete = "A", verbose = FALSE)
  expect_s3_class(res$clean_tree, "multiPhylo")
  expect_false("A" %in% res$clean_tree[[2]]$tip.label)
})
