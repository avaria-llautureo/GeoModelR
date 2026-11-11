# GeoModelR 0.3.2

## New features

* `restrict_node_map()` assigns a user map to one internal node, e.g. to
  restrict a clade's crown node to one palaeocontinent and compare
  alternative origins. The node is found from its tips (or ape node number)
  and its age from the tree (height above the most recent tip, as in
  BayesTraits); non-ultrametric trees are supported. The map can be polygons
  (`sf` object or file) or a raster (`SpatRaster` or file, thresholded with
  `mask_range`); `invert` swaps 1 and 0. It is written on the same grid as
  the series CSVs. The node gets the interval node age +/- `width` (default
  0.001 Ma), and a warning lists any other node inside that interval.
* `write_buildmaps_file()` inserts node maps and splits the series map
  around them, so intervals have no gaps or overlaps; overlapping node maps
  stop with an error. Times are now written with up to 6 decimals (existing
  lines are unchanged).
* `offset` in `write_buildmaps_file()` and `restrict_node_map()`: the age of
  the youngest tip, for all-extinct clades (BayesTraits measures node ages
  from the most recent tip). Intervals are shifted into tree time.

## Bug fixes

* `get_paleomaps()` could return an empty map for some ages and plate models
  (with the warning "x contains no geometries of specified type"). Repairing
  the coastlines can leave some features empty, which made every polygon of
  that age be dropped. Empty features are now removed first, so all polygons
  are kept. This also applies where coastlines are downloaded in
  `paleomap_to_csv()` and `fill_missing_maps()`, and to polygon maps in
  `restrict_node_map()`.
* `get_paleomaps()` now reports an age with no coastline polygons as failed
  (in `$failed` and the warning) instead of storing an empty map.
* `get_paleomaps()` start message shows the real age range (e.g. "1-10 Ma"
  for `ages = 1:10`, not "0-10 Ma").

# GeoModelR 0.3.1

## New features

* Restriction maps from local raster files. Geo model restrictions can now
  be built from continuous palaeo-elevation and bathymetry rasters (e.g. the
  Scotese PaleoDEMs), or any other gridded variable,
  not only from GPlates coastline polygons.
  - `load_raster_maps()` registers a folder of raster files, one per age
    (ages read from file names with `age_regex`; `values_var` names the
    variable), or one present-day map with `file =`. Only EPSG:4326 rasters
    are accepted (`check_crs = FALSE` overrides). `age0_override` replaces
    the series' own 0 Ma map with an observational present-day DEM.
  - `fetch_etopo2022()` downloads and caches the ETOPO 2022 60 arc-second
    bedrock relief model for use as a present-day map.
  - `paleomap_to_csv()` replaces `sf_to_csv()`. With `source = "raster"`,
    cells inside `mask_range` (one or several ranges, bounds inclusive)
    become 1 and all others 0; `invert` swaps them. Cells with no value
    become 0. `resolution` (optional) resamples every map with
    `resample_method` (`"bilinear"` or `"cubicspline"`); a warning is given
    for any map finer than 0.5 degrees written at native resolution.
  - `export_restriction_maps()` traces the `mask_range` boundary as a contour
    line for raster maps.
* `write_buildmaps_file()` gains `max_age`, which extends the last map to
  cover the root, so every node gets a restriction map. This also allows a
  single present-day map to be used for the whole tree.

## Deprecated

* `sf_to_csv()` is deprecated in favour of `paleomap_to_csv()`. It still
  works with the same arguments and defaults.

## Other changes

* Plotting functions give a clear error if given maps from
  `load_raster_maps()`; `fill_missing_maps()` likewise.
* Description updated to mention raster restriction maps.

# GeoModelR 0.3.0

## New features

* `clean_sister_geodata()` checks BayesTraits Geo model
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
* New Title and Description in `DESCRIPTION`.
* The package help page now also lists `export_restriction_maps()`.

# GeoModelR 0.2.0

* Previous release.
