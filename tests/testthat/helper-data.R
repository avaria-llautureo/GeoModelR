# Shared example data (bundled with the package). Read once per test run.
example_files <- function() {
  list(
    anc  = system.file("extdata", "Primate.AncStates.txt", package = "GeoModelR"),
    tree = system.file("extdata", "Median.trees", package = "GeoModelR")
  )
}

example_phy <- local({
  phy <- NULL
  function() {
    if (is.null(phy)) phy <<- ape::read.nexus(example_files()$tree)
    phy
  }
})

example_geo <- local({
  geo <- NULL
  function() {
    if (is.null(geo))
      geo <<- read_geo_ancstates(example_files()$anc, phy = example_phy(),
                                 verbose = FALSE)
    geo
  }
})

# Small toy tree and coordinates for clean_sister_geodata()
toy_tree <- function() ape::read.tree(text = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")

toy_single <- function() {
  data.frame(V1 = c("A", "B", "C", "D", "E"),
             V2 = c(-70, -60, -50, 10, 10.2),
             V3 = c(-30, -30, -20, 45, 45.2),
             stringsAsFactors = FALSE)
}
