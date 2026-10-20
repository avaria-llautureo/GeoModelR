# GeoModelR 0.3.0

## New features

* `clean_sister_geodata()` (with Andrew Meade) checks BayesTraits Geo model
  input coordinates for sister taxa closer than `distance_threshold`
  (default 50 m), which can compress branches to zero length under the
  variable-rates model.
  - Overlapping distribution points are deleted from the taxon with more
    points; ties are broken at random (reproducible with `seed`).
  - Overlapping single centroids are never changed automatically: the
    function stops and lists the pairs (error class
    `geo_centroid_conflict`). Resolve them with `delete = c(...)`, which
    removes taxa from the data and the tree (`$clean_tree`), or
    interactively with `ask = TRUE`.
  - Every deleted point is recorded in `$deletion_log`.
  - Tips in polytomies are also compared pairwise.

* `get_dispersal_distance()` gains `frame = c("palaeo", "plate")`.
  `"palaeo"` (default) keeps the previous behaviour: distances between the
  posterior palaeocoordinates. `"plate"` reconstructs coordinates to their
  present-day positions with rgplates before measuring, removing plate
  motion. Samples that cannot be reconstructed are kept with `NA`
  distances and counted per node in `attr(result, "reconstruction")`.
  New `model` and `verbose` arguments; results carry `attr(result, "frame")`.

* Package vignette: `vignette("GeoModelR")`, a walkthrough with the primate
  example data.

## Bug fixes

* `add_ape_nodes()` can be run more than once: it no longer creates
  `ape_node.x` / `ape_node.y` columns.
* `add_ape_nodes()` now requires an exact match between each node's tips
  and a clade in `phy`, and stops with a list of mismatched nodes. Before,
  `ape::getMRCA()` silently returned a node for a mismatched tree.
* `get_paleomaps()`: when a download failed, the empty slot was removed
  from the list, shifting later maps to the wrong ages. Failed ages now keep
  an empty (`NULL`) slot, as `fill_missing_maps()` expects.
* `get_dispersal_distance()` warned about unmatched edges whenever a
  coordinate was `NA`; it now warns only when an edge has no matching node.

## Other changes

* `geosphere` added to Imports (it was used but not declared); unused
  `httr2` removed. Minimum R version lowered from 4.5.0 to 4.1.0.
* Example data moved from `data-raw/` to `inst/extdata/`
  (`system.file("extdata", package = "GeoModelR")`).
* All examples run under `R CMD check`; examples that need the GPlates Web
  Service are wrapped in `\dontrun{}`.
* `fill_missing_maps()` is now documented; duplicated and missing argument
  documentation fixed in `plot_timeinterval()` and `plot_clade()`.
* testthat test suite and GitHub Actions `R-CMD-check` workflow added.
* `inst/CITATION` added (`citation("GeoModelR")`).
* Andrew Meade added as contributor.

# GeoModelR 0.2.0

* Previous release.
