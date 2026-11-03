# GeoModelR

<!-- badges: start -->
[![R-CMD-check](https://github.com/avaria-llautureo/GeoModelR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/avaria-llautureo/GeoModelR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

**GeoModelR** provides R tools for the **Geographical (Geo) model** in
[BayesTraits v5+](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html).
The Geo model reconstructs the posterior distribution of ancestral geographic
coordinates across the nodes of a phylogeny. The model works in continuous, spherical space, 
allowing dispersal rates to vary across branches, and optionally restricting inference to land or ocean using palaemaps.

The current package version covers three workflows.

### 1. Input data preparation *(before BayesTraits)*

| Function | Description |
|---|---|
| `clean_sister_geodata()` | Find sister taxa with overlapping coordinates (< 50 m by default) and remove overlapping distribution points to avoid zero-length scaled branches under the variable-rates model |

### 2. Palaeomap preparation *(before BayesTraits)*

| Function | Description |
|---|---|
| `get_paleomaps()` | Download all palaeomaps (using rgplates::reconstruct()) for a tree's age range |
| `fill_missing_maps()` | Re-download any ages that failed |
| `load_raster_maps()` | Register local raster maps: a time series (e.g. palaeo-elevation and bathymetry from the Scotese PaleoDEMs), or one present-day map |
| `fetch_etopo2022()` | Download and cache the ETOPO 2022 present-day elevation and bathymetry model |
| `paleomap_to_csv()` | Convert palaeomaps (GPlates or raster) to binary 1/0 CSV grids; for rasters, 1 = cells inside a chosen value range |
| `export_restriction_maps()` | Write restriction boundaries as CSVs (polygon vertices, or traced contours for rasters) |
| `restrict_node_map()` | Restrict one node to its own map (e.g. one palaeocontinent), with a narrow time interval around the node's age |
| `write_buildmaps_file()` | Write the BayesTraits `BuildMaps` input file |

`sf_to_csv()` is deprecated: use `paleomap_to_csv()`.

### 3. Output exploration *(after BayesTraits)*

| Function | Description |
|---|---|
| `read_geo_ancstates()` | Parse a `*.AncStates.txt` output file |
| `add_ape_nodes()` | Match BayesTraits node names to ape node numbers (exact clade match) |
| `plot_node()` | Posterior distribution for a single node |
| `plot_timeinterval()` | All nodes within a time window |
| `plot_clade()` | All nodes within a clade |
| `get_dispersal_distance()` | Posterior branch and root-to-tip dispersal distances, in the palaeo or plate reference frame |

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("avaria-llautureo/GeoModelR", build_vignettes = TRUE)
```

Optional packages for the tree panel in `plot_node()` and HPD contours:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages(c("deeptime", "patchwork", "ggdensity"))
```

Downloading palaeomaps and `frame = "plate"` need an internet connection
(GPlates Web Service, via [rgplates](https://CRAN.R-project.org/package=rgplates)).

---

## Quick start

A full walkthrough with the bundled primate data is in the vignette:

```r
vignette("GeoModelR")
```

### Clean input coordinates

```r
library(GeoModelR)

phy    <- ape::read.nexus("your_tree.nex")
single <- read.table("Single_Coord_Data.txt", sep = "\t")
dist   <- read.table("Distribution_Data.txt", sep = "\t")

res <- clean_sister_geodata(phy, single, dist, seed = 1)
# If sister centroids overlap, the function stops and lists the pairs.
# Choose which taxa to remove and run again:
res <- clean_sister_geodata(phy, single, dist, seed = 1,
                            delete = c("Species_a"))
res$deletion_log
res$clean_tree        # tree without the deleted taxa
```

### Prepare maps for BayesTraits

```r
maps <- get_paleomaps(max_age = ceiling(max(ape::node.depth.edgelength(phy))))
maps <- fill_missing_maps(maps)

# Inspect one map before writing them all
df <- paleomap_to_csv(age = 50, maps = maps, resolution = 2)

paleomap_to_csv(age = maps$ages, maps = maps, resolution = 1, dir = "maps_csv")
write_buildmaps_file(maps)
# Then follow the BayesTraits manual: BuildMaps / LoadMaps
```

Restrictions can also be built from local raster maps of palaeo-elevation and
bathymetry (or any other gridded variable). Cells inside `mask_range` become 1,
so the restriction can be land, ocean, or a band such as the continental shelf.
The raster files are not part of the package; download them first.

```r
scotese <- load_raster_maps(dir = "scotese_nc/", age_regex = "_([0-9]+)Ma\\.nc$",
                            values_var = "z",
                            age0_override = fetch_etopo2022())  # real present day

# Check one age, then write all ages (shelf band: -200 to 0 m)
df <- paleomap_to_csv(age = 100, maps = scotese, source = "raster",
                      mask_range = c(-200, 0), resolution = 1)
paleomap_to_csv(age = scotese$ages, maps = scotese, source = "raster",
                mask_range = c(-200, 0), resolution = 1, dir = "maps_csv")
write_buildmaps_file(scotese, max_age = 120)   # max_age: root age of the tree

# One present-day map for the whole tree (e.g. a very young clade)
present <- load_raster_maps(file = fetch_etopo2022())
paleomap_to_csv(age = 0, maps = present, source = "raster",
                mask_range = c(0, Inf), resolution = 0.5, write = TRUE)
write_buildmaps_file(present, max_age = 1.2)
```

A single node can be restricted to its own map, for example the crown node to
one palaeocontinent (your own polygon or raster map). The node map gets a
narrow time interval around the node's age, and the series is split around it:

```r
crown <- c("Homo_sapiens", "Lemur_catta")          # tips defining the node
maps_afr <- restrict_node_map(maps, map = "africa_64Ma.shp", label = "AFA_64",
                              phy = phy, tips = crown, dir = "maps_csv")
write_buildmaps_file(maps_afr, out_file = "BuildMaps_AFA.txt",
                     bin_file = "Maps_AFA.bin", max_age = 76)
```

### Explore Geo model output

```r
# Bundled example data
phy <- ape::read.nexus(system.file("extdata", "Median.trees", package = "GeoModelR"))
geo <- read_geo_ancstates(system.file("extdata", "Primate.AncStates.txt",
                                      package = "GeoModelR"), phy = phy)

maps <- get_paleomaps(geo)
plot_node(geo, "Node-00008", maps = maps)                 # crown primates
plot_node(geo, "Node-00008", phy = phy, maps = maps)      # with tree panel
plot_timeinterval(geo, age_min = 8, age_max = 10, maps = maps)
plot_clade(geo, tips = find_tips(geo, "^Homo"), maps = maps)

# Dispersal distances per branch and posterior sample
bd <- get_dispersal_distance(geo, phy, type = "branch")
bd_plate <- get_dispersal_distance(geo, phy, frame = "plate")  # plate motion removed
```

---

## Example data

`inst/extdata/` (available with `system.file("extdata", package = "GeoModelR")`)
contains a minimal example from the primate analysis in Avaria-Llautureo et
al. (2025):

- `Primate.AncStates.txt`: truncated AncStates output (50 posterior samples)
- `Median.trees`: one of the primate phylogenies analysed in the study

---

## Citation

```r
citation("GeoModelR")
```

> Avaria-Llautureo, J., Püschel, T. A., Meade, A., Baker, J., Nicholson, S. L.,
> & Venditti, C. (2025). The radiation and geographic expansion of primates
> through diverse climates. *PNAS*, 122(32), e2423833122.
> <https://doi.org/10.1073/pnas.2423833122>

---

## Feedback

Bug reports and feature requests are welcome via
[GitHub Issues](https://github.com/avaria-llautureo/GeoModelR/issues).

## License

GPL (>= 3)
