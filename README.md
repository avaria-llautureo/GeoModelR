# GeoModelR <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

---

## Overview

**GeoModelR** provides R tools for working with the **Geographical (Geo) model** implemented in [BayesTraits v5+](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.3/BayesTraitsV5.0.3.html). The Geo model reconstructs the posterior distribution of ancestral geographic coordinates across the nodes of a phylogenetic tree, optionally restricting inference to land or ocean areas using palaeomaps.

The package covers two workflows:

### Section 1 — Palaeomap Preparation *(run before BayesTraits)*
Download palaeocoastlines and generate the binary restriction map files that BayesTraits needs to run the Geo model with palaeogeographic constraints.

| Function | Description |
|---|---|
| `get_paleomaps()` | Download all palaeomaps for a tree's age range |
| `fill_missing_maps()` | Re-download any ages that failed |
| `sf_to_csv()` | Convert palaeomaps to binary land/ocean CSV grids |
| `export_restriction_maps()` | Write raw polygon CSVs |
| `write_buildmaps_file()` | Write the BayesTraits `BuildMaps` input file |

### Section 2 — Output Exploration *(run after BayesTraits)*
Read Geo model output and visualise the posterior distribution of ancestral locations on age-matched palaeomaps.

| Function | Description |
|---|---|
| `read_geo_ancstates()` | Parse a `*.AncStates.txt` output file |
| `add_ape_nodes()` | Match BayesTraits node names to ape node numbers |
| `plot_node()` | Posterior distribution for a single node |
| `plot_timeinterval()` | All nodes within a time window |
| `plot_clade()` | All nodes within a clade |

---

## Citation

If you use GeoModelR, please cite the paper that introduced the Geo model with palaeomap restrictions:

> Avaria-Llautureo, J., Püschel, H.P., et al. (2025). *The radiation and geographic expansion of primates through diverse climates*. **PNAS**, 122(18). <https://doi.org/10.1073/pnas.2423833122>

---

## Installation

GeoModelR is currently in beta. Install directly from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install GeoModelR
devtools::install_github("avaria-llautureo/GeoModelR")
```

### Required dependencies

```r
install.packages(c(
  "ape",       # phylogenetic trees
  "rgplates",  # GPlates Web Service
  "httr2",     # required by rgplates
  "sf",        # spatial features
  "terra",     # rasterisation (Section 1)
  "ggplot2"    # plotting
))

# For the tree panel in plot_node():
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages(c("deeptime", "patchwork"))

# Optional: HPD density contours
install.packages("ggdensity")
```

---

## Quick start

### Section 1 — Prepare maps for BayesTraits

```r
library(GeoModelR)
library(ape)

phy  <- read.nexus("your_tree.nex")

# Download all palaeomaps for the tree age range
maps <- get_paleomaps(max_age = ceiling(max(node.depth.edgelength(phy))))

# Fix any failed downloads
maps <- fill_missing_maps(maps)

# Write binary land/ocean CSV files
sf_to_csv(age = maps$ages, maps = maps, resolution = 1, dir = "maps_csv")

# Write the BayesTraits BuildMaps input file
write_buildmaps_file(maps)
# Then, follow the BayesTraits manual instructions to build the input map-file format for analysis
```

### Section 2 — Explore Geo model output

```r
library(GeoModelR)
library(ape)

# Read output file and link to tree
phy <- read.nexus("your_tree.nex")
geo <- read_geo_ancstates("your_analysis.AncStates.txt", phy = phy)
print(geo)

# Plot root node posterior on its age-matched palaeomap
plot_node(geo, "Node-00000") 

# With tree panel and regional zoom
plot_node(geo, "Node-00000",
                 phy  = phy,
                 bbox = c(-30, 60, -20, 60))

# All nodes in a time window
plot_timeinterval(geo, age_min = 0, age_max = 10)

# All nodes in a clade
plot_clade(geo, tips = find_tips(geo, "Homo")) # Homo clade
```

---

## Example dataset

The `data-raw/` folder contains a minimal example from the primate analysis in Avaria-Llautureo et al. (2025):
- `Primate.AncStates.txt` — truncated AncStates output (50 iterations)
- `Median.trees` — one of the hundreds of primate phylogenetic trees analysed in the study

---

## Feedback

This is a beta release. Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/avaria-llautureo/GeoModelR/issues).

---


## License

GPL (>= 3)
