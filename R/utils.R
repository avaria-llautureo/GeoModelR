# =============================================================================
# GeoModelR — Internal helpers
# These functions are used internally by prepare_maps.R and visualise.R.
# They are not exported to the user.
# =============================================================================


# =============================================================================
# Internal helpers
# =============================================================================

# Palaeomap session cache — avoids repeated GWS requests for the same age
.map_cache <- new.env(parent = emptyenv())

# Look up one age from a geo_paleomaps object, or fetch from GWS if maps = NULL
.get_palaeomap <- function(age, model = "PALEOMAP",
                            feature = "coastlines", maps = NULL) {
  age_int <- as.integer(round(age))

  if (!is.null(maps)) {
    if (!inherits(maps, "geo_paleomaps"))
      stop("'maps' must be a geo_paleomaps object from get_paleomaps().")
    key <- as.character(age_int)
    if (!key %in% names(maps$coastlines))
      stop(sprintf(
        "Age %d Ma not found in maps object.\n  Range: %d-%d Ma.\n  Rebuild with get_paleomaps(max_age = %d).",
        age_int,
        min(as.integer(names(maps$coastlines))),
        max(as.integer(names(maps$coastlines))),
        age_int
      ))
    return(maps$coastlines[[key]])
  }

  # No maps object — fetch directly from GWS
  if (!requireNamespace("rgplates", quietly = TRUE))
    stop("Package 'rgplates' is required. Install with: install.packages('rgplates')")
  rgplates::reconstruct(feature, age = age_int, model = model)
}

# Round a node height (Ma) to the nearest integer, minimum 0
.snap_age <- function(height) max(0L, as.integer(round(height)))

# Apply burn-in and return posterior coords for one node.
# node can be:
#   - A BayesTraits node name string: "Node-00005"
#   - An ape node number integer:     62  (requires ape_node column)
#
# Both paths resolve to a BayesTraits node_name first, then filter $coords
# by that name — ensuring both inputs produce identical results.
.get_node_coords <- function(geo, node, burnin = 0.1) {

  if (is.numeric(node) || is.integer(node)) {
    # --- Lookup by ape node number ---
    if (!"ape_node" %in% names(geo$nodes))
      stop("ape_node column not found in geo$nodes.\n",
           "  Run: geo <- add_ape_nodes(geo, phy)  to add it.")
    node_int <- as.integer(node)
    row_idx  <- which(!is.na(geo$nodes$ape_node) &
                        geo$nodes$ape_node == node_int)
    if (length(row_idx) == 0)
      stop("ape node number ", node_int, " not found in geo$nodes$ape_node.\n",
           "  Check with: geo$nodes[, c('node_name','ape_node')]")
    bt_name <- geo$nodes$node_name[row_idx[1]]

  } else {
    # --- Lookup by BayesTraits node name: exact match first, then partial ---
    exact <- which(geo$nodes$node_name == node)
    if (length(exact) == 1) {
      bt_name <- node
    } else {
      partial <- grep(node, geo$nodes$node_name, fixed = TRUE, value = TRUE)
      if (length(partial) == 0)
        stop("Node '", node, "' not found. Check geo$nodes$node_name.")
      if (length(partial) > 1)
        stop("'", node, "' matches multiple nodes: ",
             paste(partial, collapse = ", "), "\nBe more specific.")
      bt_name <- partial[1]
    }
  }

  # Both paths filter coords by the resolved BayesTraits node name
  coords <- geo$coords[geo$coords$node_name == bt_name, ]
  if (nrow(coords) == 0)
    stop("No coordinates found for node '", bt_name, "' in geo$coords.")

  iters <- sort(unique(coords$iteration))
  keep  <- iters[(floor(length(iters) * burnin) + 1):length(iters)]
  coords[coords$iteration %in% keep, ]
}

# Validate and unpack a bbox argument into xlim / ylim
.unpack_bbox <- function(bbox) {
  if (length(bbox) != 4 || !is.numeric(bbox))
    stop("bbox must be c(lon_min, lon_max, lat_min, lat_max).\n",
         "  Example: bbox = c(-80, -60, -45, -15)")
  list(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]))
}

# Build the ggplot2 palaeomap base layers (ocean background + land polygons)
.palaeomap_base <- function(age, model = "PALEOMAP",
                             xlim = c(-180, 180), ylim = c(-90, 90),
                             land_col   = "#d2c8a0",
                             ocean_fill = "#b8d4e8",
                             maps       = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.")

  age_int  <- .snap_age(age)
  coast_raw <- .get_palaeomap(age_int, model, "coastlines", maps = maps)

  coast_sf <- tryCatch({
    cs <- sf::st_make_valid(coast_raw)
    if (any(sf::st_geometry_type(cs) == "GEOMETRYCOLLECTION")) {
      cs <- tryCatch(
        sf::st_collection_extract(cs, "POLYGON"),
        error = function(e) cs
      )
    }
    cs
  }, error = function(e) {
    message("Note: coastline geometry repair skipped (", conditionMessage(e), ")")
    coast_raw
  })

  list(
    ggplot2::geom_sf(data = coast_sf, fill = land_col, colour = land_col,
                     linewidth = 0.2),
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                      datum = sf::st_crs(4326)),
    ggplot2::theme_bw(base_size = 11),
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = ocean_fill, colour = NA),
      panel.border     = ggplot2::element_rect(colour = "grey40", fill = NA),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 8),
      axis.ticks       = ggplot2::element_line(colour = "grey40")
    )
  )
}

# Infer a short clade label from a set of tip names (genus-level prefix)
.infer_clade_label <- function(tips) {
  if (length(tips) == 1) return(tips)
  prefixes <- unique(sub("_.*", "", tips))
  if (length(prefixes) == 1) prefixes else paste(prefixes[1], "et al.")
}

# Reconstruct palaeocoords to present (internal — see reconstruct_coords_to_present)
.reconstruct_to_present <- function(coords, model = "PALEOMAP") {

  if (!requireNamespace("rgplates", quietly = TRUE))
    stop("Package 'rgplates' is required. Install with: install.packages('rgplates')")

  coords$longitude_present <- coords$longitude
  coords$latitude_present  <- coords$latitude
  coords$.age_int <- vapply(coords$height, .snap_age, integer(1))
  ages <- sort(unique(coords$.age_int[coords$.age_int > 0]))

  if (length(ages) == 0) {
    coords$.age_int <- NULL
    return(coords)
  }

  message(sprintf(
    "[.reconstruct_to_present] Reconstructing %d age(s) to present using '%s': %s",
    length(ages), model, paste(ages, collapse = ", ")
  ))

  for (age in ages) {
    idx <- which(coords$.age_int == age)
    pts <- coords[idx, c("longitude", "latitude"), drop = FALSE]
    cache_key <- paste0("recon_present_", model, "_", age, "_",
                        nrow(pts), "_",
                        round(sum(pts$longitude), 4), "_",
                        round(sum(pts$latitude),  4))

    if (!exists(cache_key, envir = .map_cache)) {
      recon <- tryCatch(
        rgplates::reconstruct(x = pts, age = age, model = model, reverse = TRUE),
        error = function(e) {
          warning(sprintf(
            "Reconstruction failed for age %d Ma: %s\n  Keeping original coordinates.",
            age, conditionMessage(e)))
          NULL
        }
      )
      assign(cache_key, recon, envir = .map_cache)
    }

    recon <- get(cache_key, envir = .map_cache)
    if (is.null(recon)) next
    coords$longitude_present[idx] <- recon[, 1]
    coords$latitude_present[idx]  <- recon[, 2]
  }

  coords$.age_int <- NULL
  message("  Done.")
  coords
}

# Build the tree panel for plot_node()
.build_tree_panel <- function(phy, geo, node_nm, pt_col, tree_col = "grey30") {

  for (pkg in c("ggtree", "deeptime", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for the tree panel.\n",
           "  Install with: install.packages('", pkg, "')")
  }

  node_row  <- geo$nodes[geo$nodes$node_name == node_nm, ]
  node_taxa <- node_row$taxa[[1]]
  n_taxa    <- node_row$n_taxa

  # Use pre-computed ape_node if available (added by add_ape_nodes / phy arg)
  # otherwise compute on the fly
  if ("ape_node" %in% names(geo$nodes) && !is.na(node_row$ape_node)) {
    ape_node <- node_row$ape_node
  } else if (n_taxa == 1) {
    ape_node <- which(phy$tip.label == node_taxa[1])
    if (length(ape_node) == 0)
      stop("Tip '", node_taxa[1], "' not found in phy$tip.label.")
  } else {
    present <- node_taxa[node_taxa %in% phy$tip.label]
    if (length(present) < 2)
      stop("Fewer than 2 taxa from '", node_nm, "' found in phy$tip.label.")
    ape_node <- ape::getMRCA(phy, present)
  }

  tree_max <- max(ape::node.depth.edgelength(phy))
  x_min    <- -tree_max

  p_tree <- ggtree::ggtree(phy, colour = tree_col, linewidth = 0.4) |>
    ggtree::revts()

  p_tree <- p_tree +
    ggtree::geom_point2(
      mapping = ggplot2::aes(subset = (node == ape_node)),
      shape = 21, size = 4, fill = pt_col, colour = "white", stroke = 0.8
    ) +
    deeptime::coord_geo(
      xlim = c(x_min * 1.05, 0), ylim = c(0, ape::Ntip(phy) + 1),
      pos = "bottom", dat = "periods", abbrv = TRUE,
      size = 2.5, neg = TRUE, skip = NULL, expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = pretty(c(x_min, 0), n = 5),
      labels = function(x) abs(x)
    ) +
    ggplot2::labs(x = "Age (Ma)", y = NULL) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.line.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(size = 7),
      plot.margin  = ggplot2::margin(5, 10, 5, 5)
    )

  p_tree
}

# =============================================================================
# Utilities — exported helper functions for both workflows
# =============================================================================

#' List all unique node ages (rounded to integer Ma)
#' @param geo A geo_ancstates object.
get_node_ages <- function(geo) sort(unique(round(geo$nodes$height)))

#' Find tip taxon names matching a pattern
#' @param geo     A geo_ancstates object.
#' @param pattern Character. Regex or fixed string.
#' @param fixed   Logical. Fixed string match? (default FALSE).
find_tips <- function(geo, pattern, fixed = FALSE) {
  tip_taxa <- unlist(lapply(seq_len(nrow(geo$nodes)), function(i) {
    if (geo$nodes$n_taxa[i] == 1) geo$nodes$taxa[[i]] else NULL
  }))
  tip_taxa[grep(pattern, tip_taxa, fixed = fixed)]
}

#' Clear the in-session palaeomap cache
clear_map_cache <- function() {
  rm(list = ls(envir = .map_cache), envir = .map_cache)
  message("Palaeomap cache cleared.")
}

#' List palaeomaps currently held in the session cache
list_cached_maps <- function() ls(.map_cache)
