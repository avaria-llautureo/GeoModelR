test_that("ape node numbers match the tree", {
  geo <- example_geo()
  phy <- example_phy()

  expect_false(anyNA(geo$nodes$ape_node))
  expect_false(anyNA(geo$coords$ape_node))

  tips <- geo$nodes[geo$nodes$n_taxa == 1, ]
  expect_equal(phy$tip.label[tips$ape_node], unlist(tips$taxa))

  # Internal nodes: the clade below each ape node has exactly the same tips
  inner <- geo$nodes[geo$nodes$n_taxa > 1, ]
  for (i in round(seq(1, nrow(inner), length.out = 25))) {
    desc <- ape::extract.clade(phy, inner$ape_node[i])$tip.label
    expect_setequal(desc, inner$taxa[[i]])
  }
})

test_that("add_ape_nodes() can be run twice without duplicating columns", {
  geo  <- example_geo()
  phy  <- example_phy()
  geo2 <- add_ape_nodes(add_ape_nodes(geo, phy, verbose = FALSE), phy,
                        verbose = FALSE)

  expect_identical(grep("^ape_node", names(geo2$coords), value = TRUE), "ape_node")
  expect_identical(grep("^ape_node", names(geo2$nodes), value = TRUE), "ape_node")
  expect_identical(geo2$coords$ape_node, geo$coords$ape_node)
  expect_identical(geo2$coords$node_name, geo$coords$node_name)
})

test_that("old ape_node.x / ape_node.y columns are removed", {
  geo <- example_geo()
  geo$coords$ape_node.x <- geo$coords$ape_node
  geo$coords$ape_node.y <- geo$coords$ape_node
  geo$coords$ape_node   <- NULL
  geo2 <- add_ape_nodes(geo, example_phy(), verbose = FALSE)
  expect_identical(grep("^ape_node", names(geo2$coords), value = TRUE), "ape_node")
})

test_that("a tree with a different topology is rejected", {
  geo <- example_geo()
  bad <- example_phy()
  i <- match(c("Lemur_catta", "Alouatta_caraya"), bad$tip.label)
  bad$tip.label[i] <- bad$tip.label[rev(i)]
  expect_error(add_ape_nodes(geo, bad, verbose = FALSE),
               "no clade in 'phy' with exactly the same tips")
})

test_that("a tree with different tips is rejected", {
  geo <- example_geo()
  expect_error(add_ape_nodes(geo, ape::drop.tip(example_phy(), "Lemur_catta"),
                             verbose = FALSE),
               "does not match")
  expect_error(add_ape_nodes(geo, "not a tree"), "phylo")
})
