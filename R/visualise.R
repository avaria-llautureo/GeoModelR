# =============================================================================
# GeoModelR — Visualisation Functions
# Functions to explore BayesTraits Geo model output (*.AncStates.txt files).
# Run these AFTER your BayesTraits analysis.
#
# plot_node()    Posterior distribution for a single node
# plot_timeinterval() All nodes within a time window
# plot_clade()   All nodes within a clade
#
# Dependencies: rgplates, sf, ggplot2
#               ggtree, deeptime, patchwork (tree panel in plot_node)
#               ggdensity (optional, for HPD contours)
# =============================================================================

#' Plot the posterior distribution of coordinates for a single node
#'
#' @param geo          A geo_ancstates object.
#' @param node         Character. Node name or unique partial match.
#' @param burnin       Numeric 0-1. Burn-in fraction (default 0.1).
#' @param model        Character. rgplates model (default "PALEOMAP").
#' @param bbox         Numeric c(lon_min, lon_max, lat_min, lat_max).
#'                     NULL (default) = full world map.
#' @param show_points  Logical. Plot MCMC sample as points (default TRUE).
#' @param pt_col       Character. Point colour (default "#c0392b").
#' @param pt_size      Numeric. Point size (default 0.8).
#' @param pt_alpha     Numeric 0-1. Point transparency (default 0.15).
#' @param show_density Logical. Overlay 50%/95% HPD contours (default FALSE).
#'                     Requires ggdensity.
#' @param show_median  Logical. Mark posterior median with a diamond (default FALSE).
#' @param phy          A phylo object (ape). Adds a tree panel beside the map.
#' @param tree_width   Numeric. Relative width of the tree panel (default 0.5).
#' @param maps         A geo_paleomaps object. Avoids GWS calls.
#' @param title        Character or NULL. Custom plot title.
#' @return A ggplot2 or patchwork object.
#'
#' @examples
#' geo <- read_geo_ancstates("Primate.AncStates.txt")
#' plot_node(geo, "Node-00008")
#' plot_node(geo, "Node-00008", bbox = c(-82, -33, -56, -15))
#' plot_node(geo, "Node-00008", phy = phy, maps = maps)
plot_node <- function(geo,
                              node,
                              burnin       = 0.1,
                              model        = "PALEOMAP",
                              bbox         = NULL,
                              show_points  = TRUE,
                              pt_col       = "#c0392b",
                              pt_size      = 0.8,
                              pt_alpha     = 0.15,
                              show_density = FALSE,
                              show_median  = FALSE,
                              phy          = NULL,
                              tree_width   = 0.5,
                              maps         = NULL,
                              title        = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  coords  <- .get_node_coords(geo, node, burnin)
  node_nm <- unique(coords$node_name)
  height  <- unique(coords$height)
  is_tip  <- unique(coords$is_tip)
  age_int <- .snap_age(height)

  # Print tip names to console
  node_row  <- geo$nodes[geo$nodes$node_name == node_nm, ]
  node_taxa <- node_row$taxa[[1]]
  cat("\n--- Tips associated with", node_nm, "---\n")
  cat(paste0("  ", seq_along(node_taxa), ". ", node_taxa, collapse = "\n"), "\n\n")

  ext <- if (!is.null(bbox)) .unpack_bbox(bbox) else
           list(xlim = c(-180, 180), ylim = c(-90, 90))

  med_df <- data.frame(
    longitude = median(coords$longitude, na.rm = TRUE),
    latitude  = median(coords$latitude,  na.rm = TRUE)
  )

  p_map <- ggplot2::ggplot() +
    .palaeomap_base(age_int, model, ext$xlim, ext$ylim, maps = maps)

  if (show_points) {
    p_map <- p_map + ggplot2::geom_point(
      data = coords, ggplot2::aes(x = longitude, y = latitude),
      colour = pt_col, alpha = pt_alpha, size = pt_size, shape = 16
    )
  }

  if (show_density) {
    if (!requireNamespace("ggdensity", quietly = TRUE)) {
      warning("Package 'ggdensity' not installed — density layer skipped.")
    } else {
      # geom_hdr default behaviour: varying alpha per density level,
      # fill set to pt_col. Simple and reliable.
      p_map <- p_map +
        ggdensity::geom_hdr(
          data        = coords,
          mapping     = ggplot2::aes(x = longitude, y = latitude),
          probs       = c(0.99, 0.95, 0.80, 0.50),
          fill        = pt_col,
          show.legend = FALSE
        )
    }
  }

  if (show_median) {
    p_map <- p_map + ggplot2::geom_point(
      data = med_df, ggplot2::aes(x = longitude, y = latitude),
      colour = "white", fill = pt_col, size = 3.5, shape = 23
    )
  }

  node_label <- if (is_tip) paste0(node_nm, " (tip)") else node_nm
  age_label  <- if (age_int == 0) "Present" else paste0(age_int, " Ma")

  p_map <- p_map + ggplot2::labs(
    title    = if (!is.null(title)) title else
                 paste0("Posterior location - ", node_label,
                        "  [", age_label, "]  n = ", nrow(coords), " iterations"),
    subtitle = paste0("Palaeomap: ", model, "  |  burn-in: ", burnin * 100, "%"),
    x = "Palaeolongitude (deg)", y = "Palaeolatitude (deg)"
  )

  if (is.null(phy)) return(p_map)

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for the tree panel.")

  p_tree <- withCallingHandlers(
    .build_tree_panel(phy, geo, node_nm, pt_col),
    warning = function(w) {
      if (grepl("subgroup", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    }
  )

  result <- patchwork::wrap_plots(p_map, p_tree, ncol = 2,
                                   widths = c(1, tree_width), guides = "collect")

  withCallingHandlers(
    print(result),
    warning = function(w) {
      if (grepl("subgroup", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    }
  )
  invisible(result)
}

#' Plot posterior node locations for all nodes within a time window
#'
#' @param geo           A geo_ancstates object.
#' @param age_min       Numeric. Lower bound of the time window in Ma (default 0).
#' @param age_max       Numeric. Upper bound in Ma (default Inf = all nodes).
#' @param burnin        Numeric 0-1. Burn-in fraction (default 0.1).
#' @param model         Character. rgplates model (default "PALEOMAP").
#' @param modern_coords Logical. If TRUE (default), reconstruct coordinates to
#'                      present day to avoid land/sea mismatch across ages.
#' @param map_age       Numeric or NULL. Used when modern_coords = FALSE.
#'                      Palaeomap age in Ma; NULL = midpoint of the window.
#' @param show_tips     Logical. Include tip nodes? Default TRUE.
#' @param bbox          Numeric c(lon_min, lon_max, lat_min, lat_max).
#'                      NULL (default) = full world map.
#' @param show_points   Logical. Plot full MCMC sample as points (default FALSE).
#' @param pt_col        Character. Point colour (default "#555555").
#' @param pt_size       Numeric. Point size (default 0.4).
#' @param pt_alpha      Numeric 0-1. Point transparency (default 0.15).
#' @param show_tips     Logical. Include tip nodes in the plot (default TRUE).
#'                      Set FALSE to show only internal nodes.
#' @param show_median   Logical. Plot posterior median per node (default TRUE).
#' @param med_size      Numeric. Median symbol size (default 3).
#' @param maps          A geo_paleomaps object. Avoids GWS calls.
#' @return A ggplot2 object.
#'
#' @examples
#' geo <- read_geo_ancstates("VR.AncStates.txt")
#' plot_timeinterval(geo, age_min = 0, age_max = 8)
#' plot_timeinterval(geo, age_min = 0, age_max = 8, modern_coords = FALSE, maps = maps)
plot_timeinterval <- function(geo,
                                 age_min       = 0,
                                 age_max       = Inf,
                                 burnin        = 0.1,
                                 model         = "PALEOMAP",
                                 modern_coords = TRUE,
                                 map_age       = NULL,
                                 show_tips     = TRUE,
                                 bbox          = NULL,
                                 show_points   = FALSE,
                                 pt_col        = "#555555",
                                 pt_size       = 0.4,
                                 pt_alpha      = 0.15,
                                 show_median   = TRUE,
                                 med_size      = 3,
                                 maps          = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  in_bin    <- geo$nodes$height >= age_min &
               geo$nodes$height <= age_max &
               (show_tips | geo$nodes$n_taxa > 1)
  bin_nodes <- geo$nodes$node_name[in_bin]

  if (length(bin_nodes) == 0)
    stop("No nodes found between ", age_min, " and ", age_max, " Ma.")

  message(sprintf("[plot_timeinterval] %d nodes in [%.2f, %.2f] Ma",
                  length(bin_nodes), age_min, age_max))

  iters      <- sort(unique(geo$coords$iteration))
  keep       <- iters[(floor(length(iters) * burnin) + 1):length(iters)]
  coords_bin <- geo$coords[geo$coords$node_name %in% bin_nodes &
                             geo$coords$iteration %in% keep, ]

  if (modern_coords) {
    if (!"longitude_present" %in% names(coords_bin)) {
      coords_bin <- .reconstruct_to_present(coords_bin, model)
    } else {
      message("[plot_timeinterval] Using pre-computed present-day coordinates.")
    }
    lon_col <- "longitude_present"; lat_col <- "latitude_present"
    map_age_used <- 0L
    coord_label  <- "Modern longitude (deg)"; lat_label <- "Modern latitude (deg)"
    subtitle_map <- "Modern map (coordinates reconstructed to present)"
  } else {
    lon_col <- "longitude"; lat_col <- "latitude"
    if (is.null(map_age))
      map_age <- .snap_age((age_min + min(age_max, max(geo$nodes$height))) / 2)
    map_age_used <- .snap_age(map_age)
    coord_label  <- "Palaeolongitude (deg)"; lat_label <- "Palaeolatitude (deg)"
    subtitle_map <- paste0("Palaeomap at ", map_age_used, " Ma")
  }

  medians <- do.call(rbind, lapply(bin_nodes, function(nd) {
    s <- coords_bin[coords_bin$node_name == nd, ]
    data.frame(node_name = nd,
               longitude = median(s[[lon_col]], na.rm = TRUE),
               latitude  = median(s[[lat_col]], na.rm = TRUE),
               height    = unique(s$height),
               is_tip    = unique(s$is_tip),
               stringsAsFactors = FALSE)
  }))

  ext <- if (!is.null(bbox)) .unpack_bbox(bbox) else
           list(xlim = c(-180, 180), ylim = c(-90, 90))

  p <- ggplot2::ggplot() +
    .palaeomap_base(map_age_used, model, ext$xlim, ext$ylim, maps = maps)

  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = coords_bin,
      ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]]),
      colour = pt_col, alpha = pt_alpha, size = pt_size, shape = 16
    )
  }

  if (show_median) {
    p <- p +
      ggplot2::geom_point(
        data    = medians,
        mapping = ggplot2::aes(x = longitude, y = latitude,
                               colour = is_tip, shape = is_tip),
        size = med_size
      ) +
      ggplot2::scale_colour_manual(
        name = "Node type",
        values = c("TRUE" = "#2980b9", "FALSE" = "#c0392b"),
        labels = c("TRUE" = "Tip", "FALSE" = "Internal")
      ) +
      ggplot2::scale_shape_manual(
        name = "Node type",
        values = c("TRUE" = 16, "FALSE" = 18),
        labels = c("TRUE" = "Tip", "FALSE" = "Internal")
      )
  }

  p + ggplot2::labs(
    title    = paste0("Posterior node locations - ", age_min, "-",
                       ifelse(is.infinite(age_max), "root", age_max), " Ma"),
    subtitle = paste0(subtitle_map, "  |  ", length(bin_nodes),
                       " nodes  |  Model: ", model),
    x = coord_label, y = lat_label
  )
}

#' Plot posterior node locations for all nodes within a clade
#'
#' @param geo           A geo_ancstates object.
#' @param tips          Character vector. Tip names defining the clade.
#' @param burnin        Numeric 0-1. Burn-in fraction (default 0.1).
#' @param model         Character. rgplates model (default "PALEOMAP").
#' @param modern_coords Logical. If TRUE (default), reconstruct coordinates to
#'                      present day to avoid land/sea mismatch.
#' @param map_age       Numeric or NULL. Used when modern_coords = FALSE.
#'                      NULL = clade root age.
#' @param bbox          Numeric c(lon_min, lon_max, lat_min, lat_max).
#'                      NULL (default) = full world map.
#' @param show_points   Logical. Plot full MCMC sample as points (default FALSE).
#' @param pt_col        Character. Point colour (default "#555555").
#' @param pt_size       Numeric. Point size (default 0.4).
#' @param pt_alpha      Numeric 0-1. Point transparency (default 0.15).
#' @param show_median   Logical. Plot posterior median per node (default TRUE).
#' @param med_size      Numeric. Median symbol size (default 3).
#' @param label_nodes   Logical. Add node name labels (default FALSE).
#' @param maps          A geo_paleomaps object. Avoids GWS calls.
#' @return A ggplot2 object.
#'
#' @examples
#' geo <- read_geo_ancstates("VR.AncStates.txt")
#' plot_clade(geo, tips = c("A_ligtu_ligtu_PCM32", "A_ligtu_simsii_OTN173"))
#' plot_clade(geo, tips = c("A_ligtu_ligtu_PCM32", "A_ligtu_simsii_OTN173"),
#'                   modern_coords = FALSE, maps = maps)
plot_clade <- function(geo,
                               tips,
                               burnin        = 0.1,
                               model         = "PALEOMAP",
                               modern_coords = TRUE,
                               map_age       = NULL,
                               bbox          = NULL,
                               show_points   = FALSE,
                               pt_col        = "#555555",
                               pt_size       = 0.4,
                               pt_alpha      = 0.15,
                               show_tips     = TRUE,
                               show_median   = TRUE,
                               med_size      = 3,
                               label_nodes   = FALSE,
                               maps          = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  # A node belongs to the clade if ALL of its taxa are within the set of
  # clade-defining tips (subset condition). The previous condition
  # all(tips %in% tl) was inverted — it found nodes that CONTAIN all tips,
  # which includes every ancestor up to the root of the whole tree.
  clade_mask  <- sapply(geo$nodes$taxa, function(tl) {
    length(tl) > 0 && all(tl %in% tips)
  })
  clade_nodes <- geo$nodes$node_name[clade_mask]

  # Optionally exclude tip nodes
  if (!show_tips)
    clade_nodes <- clade_nodes[geo$nodes$n_taxa[clade_mask] > 1]

  if (length(clade_nodes) == 0)
    stop("No nodes found containing all of: ", paste(tips, collapse = ", "),
         "\nCheck tip names with: geo$nodes$taxa[[1]]")

  message(sprintf("[plot_clade] Clade defined by %d tip(s): %d nodes found (%s)",
                  length(tips), length(clade_nodes),
                  if (show_tips) "tips + internal" else "internal only"))

  iters        <- sort(unique(geo$coords$iteration))
  keep         <- iters[(floor(length(iters) * burnin) + 1):length(iters)]
  coords_clade <- geo$coords[geo$coords$node_name %in% clade_nodes &
                               geo$coords$iteration %in% keep, ]

  if (modern_coords) {
    if (!"longitude_present" %in% names(coords_clade)) {
      coords_clade <- .reconstruct_to_present(coords_clade, model)
    } else {
      message("[plot_clade] Using pre-computed present-day coordinates.")
    }
    lon_col <- "longitude_present"; lat_col <- "latitude_present"
    map_age_used <- 0L
    coord_label  <- "Modern longitude (deg)"; lat_label <- "Modern latitude (deg)"
    subtitle_map <- "Modern map (coordinates reconstructed to present)"
  } else {
    lon_col <- "longitude"; lat_col <- "latitude"
    if (is.null(map_age)) map_age <- max(geo$nodes$height[clade_mask])
    map_age_used <- .snap_age(map_age)
    coord_label  <- "Palaeolongitude (deg)"; lat_label <- "Palaeolatitude (deg)"
    subtitle_map <- paste0("Palaeomap at ", map_age_used, " Ma")
  }

  medians <- do.call(rbind, lapply(clade_nodes, function(nd) {
    s <- coords_clade[coords_clade$node_name == nd, ]
    data.frame(node_name = nd,
               longitude = median(s[[lon_col]], na.rm = TRUE),
               latitude  = median(s[[lat_col]], na.rm = TRUE),
               height    = unique(s$height),
               is_tip    = unique(s$is_tip),
               stringsAsFactors = FALSE)
  }))

  ext <- if (!is.null(bbox)) .unpack_bbox(bbox) else
           list(xlim = c(-180, 180), ylim = c(-90, 90))

  p <- ggplot2::ggplot() +
    .palaeomap_base(map_age_used, model, ext$xlim, ext$ylim, maps = maps)

  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = coords_clade,
      ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]]),
      colour = pt_col, alpha = pt_alpha, size = pt_size, shape = 16
    )
  }

  if (show_median) {
    p <- p +
      ggplot2::geom_point(
        data    = medians,
        mapping = ggplot2::aes(x = longitude, y = latitude,
                               colour = height, shape = is_tip),
        size = med_size
      ) +
      ggplot2::scale_colour_viridis_c(
        name = "Node age (Ma)", option = "magma", direction = -1
      ) +
      ggplot2::scale_shape_manual(
        name   = "Node type",
        values = c("TRUE" = 16, "FALSE" = 18),
        labels = c("TRUE" = "Tip", "FALSE" = "Internal")
      )

    if (label_nodes) {
      p <- p + ggplot2::geom_text(
        data    = medians,
        mapping = ggplot2::aes(x = longitude, y = latitude, label = node_name),
        size = 2.5, colour = "grey20", nudge_y = 0.8
      )
    }
  }

  p + ggplot2::labs(
    title    = paste0("Clade: ", .infer_clade_label(tips)),
    subtitle = paste0(subtitle_map, "  |  ", length(clade_nodes),
                       " nodes  |  Model: ", model),
    x = coord_label, y = lat_label
  )
}
