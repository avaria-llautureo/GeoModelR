# =============================================================================
# GeoModelR — Palaeomap Preparation Functions
# Functions to build geographic input files for BayesTraits Geo model analysis.
# Run these BEFORE your BayesTraits analysis.
#
# get_paleomaps()            Download palaeomaps from GPlates Web Service
# fill_missing_maps()        Re-download any ages that failed
# sf_to_csv()                Convert palaeomaps to binary land/ocean CSV grids
# export_restriction_maps()  Write raw polygon CSVs
# write_buildmaps_file()     Write the BayesTraits BuildMaps input file
#
# Dependencies: rgplates, sf, terra
# =============================================================================

# =============================================================================
# SECTION 1 — PALAEOMAP PREPARATION
# Functions to build geographic input files for the BayesTraits Geo model.
# Run these before your BayesTraits analysis.
# =============================================================================

#' Download and store palaeomaps for a Geo model analysis
#'
#' Fetches palaeocoastlines from the GPlates Web Service for every integer age
#' from 0 Ma to \code{max_age} Ma (1 Ma resolution). All maps are downloaded
#' once and stored in a \code{geo_paleomaps} object. Pass this object to
#' \code{sf_to_csv()} to generate BayesTraits input files, or to any plotting
#' function via the \code{maps} argument to avoid repeated GWS calls.
#'
#' @param geo     A \code{geo_ancstates} object. Used to determine
#'                \code{max_age} automatically. Supply either \code{geo}
#'                or \code{max_age}.
#' @param max_age Numeric. Maximum age in Ma. Overrides \code{geo}.
#' @param model   Character. rgplates plate model (default \code{"PALEOMAP"}).
#'                For alternative paleogeographic models look at the
#'                rgplates::reconstruction() function.
#' @param ages    Integer vector or NULL. Explicit ages to download.
#'                NULL (default) downloads every integer from 0 to max_age.
#' @return A \code{geo_paleomaps} object. Pass to plotting functions via
#'         the \code{maps} argument.
#'
#' @examples
#' # From a tree (no BayesTraits output needed for data preparation)
#' library(ape)
#' phy  <- read.nexus("Median.trees")
#' maps <- get_paleomaps(max_age = ceiling(max(node.depth.edgelength(phy))))
#'
#' # Or from a geo_ancstates object
#' geo  <- read_geo_ancstates("Geo_Primates_AncStates.txt")
#' maps <- get_paleomaps(geo)
#'
#' # Use in plotting
#' plot_node_palaeo(geo, "Node-00000", maps = maps)
get_paleomaps <- function(geo     = NULL,
                           max_age = NULL,
                           model   = "PALEOMAP",
                           ages    = NULL) {

  if (!requireNamespace("rgplates", quietly = TRUE))
    stop("Package 'rgplates' is required. Install with: install.packages('rgplates')")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.")

  if (is.null(ages)) {
    if (is.null(max_age)) {
      if (is.null(geo) || !inherits(geo, "geo_ancstates"))
        stop("Supply either a geo_ancstates object or max_age.")
      max_age <- ceiling(max(geo$nodes$height))
    }
    ages <- seq(0L, as.integer(max_age), by = 1L)
  } else {
    ages    <- as.integer(sort(unique(ages)))
    max_age <- max(ages)
  }

  n <- length(ages)
  message(sprintf(
    "[get_paleomaps] Downloading %d palaeomaps (0-%d Ma), model '%s'...",
    n, max_age, model
  ))
  message("  This may take a few minutes. Maps are stored locally for reuse.")

  coastlines <- vector("list", n)
  names(coastlines) <- as.character(ages)
  failed <- integer(0)

  for (i in seq_along(ages)) {
    age <- ages[i]
    if (i %% 10 == 0 || i == n)
      message(sprintf("  ... %d / %d ages downloaded", i, n))

    coastlines[[i]] <- tryCatch({
      raw <- rgplates::reconstruct("coastlines", age = age, model = model)
      cs <- sf::st_make_valid(raw)
      if (any(sf::st_geometry_type(cs) == "GEOMETRYCOLLECTION")) {
        cs <- tryCatch(
          sf::st_collection_extract(cs, "POLYGON"),
          error = function(e) cs   # already clean — skip extraction
        )
      }
      cs
    }, error = function(e) {
      message(sprintf("  WARNING: age %d Ma failed (%s) — skipping.",
                      age, conditionMessage(e)))
      failed <<- c(failed, age)
      NULL
    })

    if (age > 0 && i < n) Sys.sleep(0.3)   # avoid GWS rate-limiting
  }

  if (length(failed) > 0)
    warning(sprintf(
      "%d age(s) failed: %s\n  Re-run get_paleomaps() when the GWS is available.",
      length(failed), paste(failed, collapse = ", ")
    ))

  result <- structure(
    list(coastlines = coastlines, model = model,
         max_age = as.integer(max_age), ages = ages,
         n_ages = n - length(failed), failed = failed),
    class = "geo_paleomaps"
  )

  message(sprintf("[get_paleomaps] Done. %d/%d maps downloaded.", n - length(failed), n))
  result
}

#' Print method for geo_paleomaps
print.geo_paleomaps <- function(x, ...) {
  cat("Palaeomap collection (geo_paleomaps)\n")
  cat("  Model  :", x$model, "\n")
  cat("  Range  : 0 -", x$max_age, "Ma\n")
  cat("  Maps   :", x$n_ages, "ages downloaded\n")
  if (length(x$failed) > 0)
    cat("  Failed :", paste(x$failed, collapse = ", "), "Ma\n")
  cat("\nPass to plotting functions via maps = or to sf_to_csv().\n")
  invisible(x)
}

#' Convert palaeomaps to binary land/ocean CSV files for BayesTraits
#'
#' Rasterises palaeocoastlines into a binary land (1) / ocean (0) grid and
#' writes CSV files in the format required by the BayesTraits \code{BuildMaps}
#' command. For a single age, returns a data frame for inspection. For multiple
#' ages, writes one file per age to \code{dir}.
#'
#' @param age        Numeric vector. Age(s) in Ma. Default \code{0}.
#' @param maps       A \code{geo_paleomaps} object from \code{get_paleomaps()}.
#'                   Recommended for batch use — avoids re-downloading.
#' @param model      Character. rgplates model (default \code{"PALEOMAP"}).
#'                   Ignored when \code{maps} is supplied.
#' @param extent     Numeric \code{c(xmin, xmax, ymin, ymax)}.
#'                   Default: global \code{c(-180, 180, -90, 90)}.
#' @param resolution Numeric. Grid cell size in degrees (default \code{1}).
#'                   Use smaller values for analysis; 1-2 is suitable for testing.
#' @param dir        Character. Output directory. Default \code{"maps_csv"}.
#' @param write      Logical. Force file output for a single age (default FALSE).
#' @param prefix     Character. File name prefix (default \code{"mask_"}).
#' @param verbose    Logical. Print progress (default TRUE).
#' @return Single age + write=FALSE: a data frame (Lon, Lat, mask).
#'         Otherwise: invisibly returns the output directory path.
#'
#' @examples
#' # Inspect one map before committing to full resolution
#' df <- sf_to_csv(age = 50)
#' library(ggplot2)
#' ggplot(df, aes(x = Lon, y = Lat, fill = factor(mask))) +
#'   geom_raster() +
#'   scale_fill_manual(values = c("0" = "#b8d4e8", "1" = "#d2c8a0"),
#'                     labels = c("Ocean", "Land"), name = NULL) +
#'   coord_fixed(expand = FALSE) + theme_bw()
#'
#' # Batch: write all ages from a maps object
#' maps <- get_paleomaps(max_age = 74)
#' sf_to_csv(age = maps$ages, maps = maps, resolution = 1, dir = "maps_csv")
sf_to_csv <- function(age        = 0,
                       maps       = NULL,
                       model      = "PALEOMAP",
                       extent     = c(-180, 180, -90, 90),
                       resolution = 1,
                       dir        = "maps_csv",
                       write      = FALSE,
                       prefix     = "mask_",
                       verbose    = TRUE) {

  if (!is.numeric(age) || any(age < 0))
    stop("'age' must be a non-negative numeric vector.")
  if (length(extent) != 4 || !is.numeric(extent))
    stop("'extent' must be c(xmin, xmax, ymin, ymax).")
  if (extent[1] < -180 || extent[2] > 180 || extent[3] < -90 || extent[4] > 90)
    stop("'extent' values outside valid range: lon [-180, 180], lat [-90, 90].")
  if (!is.numeric(resolution) || resolution <= 0)
    stop("'resolution' must be a positive number.")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install with: install.packages('terra')")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.")
  if (!is.null(maps) && !inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps().")

  age       <- sort(unique(as.integer(round(age))))
  multi_age <- length(age) > 1 || write
  if (multi_age && !dir.exists(dir)) dir.create(dir, recursive = TRUE)

  # Build raster template once — reused for all ages
  grid_template <- terra::rast(
    xmin = extent[1], xmax = extent[2],
    ymin = extent[3], ymax = extent[4],
    res  = resolution, crs = "EPSG:4326"
  )

  if (verbose && multi_age)
    message(sprintf("[sf_to_csv] Writing %d map(s) at %.2f deg resolution to: %s/",
                    length(age), resolution, dir))

  results <- vector("list", length(age))
  names(results) <- as.character(age)

  for (i in seq_along(age)) {
    a <- age[i]
    if (verbose && multi_age && (i %% 10 == 0 || i == length(age)))
      message(sprintf("  ... %d / %d", i, length(age)))

    coast_sf <- tryCatch({
      if (!is.null(maps)) {
        key <- as.character(a)
        if (!key %in% names(maps$coastlines))
          stop(sprintf("Age %d Ma not in maps object.", a))
        maps$coastlines[[key]]
      } else {
        if (!requireNamespace("rgplates", quietly = TRUE))
          stop("Package 'rgplates' is required.")
        raw <- rgplates::reconstruct("coastlines", age = a, model = model)
        cs  <- sf::st_make_valid(raw)
        if (any(sf::st_geometry_type(cs) == "GEOMETRYCOLLECTION")) {
          cs <- tryCatch(
            sf::st_collection_extract(cs, "POLYGON"),
            error = function(e) cs
          )
        }
        cs
      }
    }, error = function(e) {
      warning(sprintf("Age %d Ma: could not get coastlines (%s) — skipping.",
                      a, conditionMessage(e)))
      NULL
    })
    if (is.null(coast_sf)) next

    # terra::rasterize requires a plain sf — strip any extra classes
    class(coast_sf) <- c("sf", "data.frame")

    map_rast <- tryCatch(
      terra::rasterize(coast_sf, grid_template),
      error = function(e) {
        warning(sprintf("Age %d Ma: rasterization failed (%s) — skipping.",
                        a, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(map_rast)) next

    # Binary: land = 1, ocean = 0
    vals <- terra::values(map_rast)
    terra::values(map_rast) <- ifelse(is.na(vals), 0L, 1L)

    df <- as.data.frame(map_rast, xy = TRUE, na.rm = FALSE)
    colnames(df) <- c("Lon", "Lat", "mask")

    # Enforce BayesTraits coordinate limits
    df <- df[df$Lon >= -180 & df$Lon < 180 &
             df$Lat >  -90  & df$Lat <  90, ]

    if (multi_age || write) {
      write.csv(df, file = file.path(dir, sprintf("%s%d.csv", prefix, a)),
                row.names = FALSE, quote = FALSE)
    } else {
      results[[1]] <- df
    }
  }

  if (multi_age || write) {
    if (verbose)
      message(sprintf("[sf_to_csv] Done. Files written to: %s/",
                      normalizePath(dir)))
    invisible(normalizePath(dir))
  } else {
    results[[1]]
  }
}

#' Export palaeomaps as restriction map CSVs for BayesTraits
#'
#' Writes one CSV file per age containing coastline polygon vertex coordinates
#' in the BayesTraits restriction map format. Each file is named
#' \code{mask_N.csv} where N is the age in Ma.
#'
#' Note: \code{sf_to_csv()} produces rasterised binary grids (Lon/Lat/mask),
#' which is the format required by BayesTraits \code{BuildMaps}. This function
#' exports raw polygon vertices instead, which may be useful for custom
#' processing or visualisation.
#'
#' @param maps    A \code{geo_paleomaps} object from \code{get_paleomaps()}.
#' @param dir     Output directory. Default \code{"restriction_maps"}.
#' @param ages    Integer vector or NULL. Subset of ages. NULL = all ages.
#' @param verbose Logical. Print progress (default TRUE).
#' @return Invisibly returns the output directory path.
#'
#' @examples
#' maps <- get_paleomaps(geo)
#' export_restriction_maps(maps)
#' export_restriction_maps(maps, ages = get_node_ages(geo))
export_restriction_maps <- function(maps,
                                     dir     = "restriction_maps",
                                     ages    = NULL,
                                     verbose = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps().")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.")

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  export_ages <- if (is.null(ages)) maps$ages else as.integer(ages)
  export_ages <- export_ages[export_ages %in% maps$ages]
  if (length(export_ages) == 0) stop("No matching ages found in the maps object.")

  n_written <- 0L

  for (age in export_ages) {
    coast_sf <- maps$coastlines[[as.character(age)]]
    if (is.null(coast_sf)) {
      if (verbose) message(sprintf("  Skipping age %d Ma (download failed).", age))
      next
    }

    coords_list <- lapply(seq_len(nrow(coast_sf)), function(i) {
      coords <- tryCatch(
        as.data.frame(sf::st_coordinates(sf::st_geometry(coast_sf[i, ]))[, 1:2]),
        error = function(e) NULL
      )
      if (is.null(coords) || nrow(coords) == 0) return(NULL)
      names(coords) <- c("longitude", "latitude")
      coords
    })
    coords_list <- Filter(Negate(is.null), coords_list)
    if (length(coords_list) == 0) next

    all_coords <- do.call(rbind, coords_list)
    write.table(all_coords,
                file      = file.path(dir, sprintf("mask_%d.csv", age)),
                sep       = ",", row.names = FALSE, col.names = FALSE,
                quote     = FALSE)
    n_written <- n_written + 1L
  }

  if (verbose)
    message(sprintf("[export_restriction_maps] %d files written to: %s/",
                    n_written, normalizePath(dir)))
  invisible(normalizePath(dir))
}


#' Write a BayesTraits BuildMaps input file
#'
#' Creates the text file required by the BayesTraits \code{BuildMaps} command
#' to compile individual CSV map files into a binary map file. Each CSV map
#' is assigned a time interval following the convention used in BayesTraits:
#'
#' \itemize{
#'   \item \code{mask_0.csv}   covers \code{[0, 0.5)} Ma
#'   \item \code{mask_1.csv}   covers \code{[0.5, 1.5)} Ma
#'   \item \code{mask_N.csv}   covers \code{[N-0.5, N+0.5)} Ma  (middle maps)
#'   \item \code{mask_MAX.csv} covers \code{[MAX-0.5, MAX+0.5)} Ma
#' }
#'
#' The output file can be passed directly to BayesTraits with the
#' \code{BuildMaps} command. After \code{BuildMaps} has run, load the resulting
#' binary file in BayesTraits with \code{LoadMaps BinaryMap.bin}.
#'
#' @param maps       A \code{geo_paleomaps} object from \code{get_paleomaps()}.
#'                   Used to determine the set of available ages.
#' @param bin_file   Character. Name of the binary map file that BayesTraits
#'                   will create (default \code{"BinaryMap.bin"}).
#' @param prefix     Character. CSV file name prefix (default \code{"mask_"}).
#'                   Must match the prefix used in \code{sf_to_csv()}.
#' @param lon_boxes  Integer. Number of longitude grid boxes for BayesTraits
#'                   spatial indexing (default 250, as recommended in the manual).
#' @param lat_boxes  Integer. Number of latitude grid boxes (default 150).
#' @param out_file   Character. Output file path (default \code{"BuildMaps.txt"}).
#' @param ages       Integer vector or NULL. Subset of ages to include.
#'                   NULL (default) uses all ages in the maps object.
#' @param verbose    Logical. Print a summary on completion (default TRUE).
#' @return Invisibly returns the output file path.
#'
#' @examples
#' maps <- get_paleomaps(geo)
#' sf_to_csv(age = maps$ages, maps = maps, dir = "maps_csv")
#'
#' # Write the BuildMaps file
#' write_buildmaps_file(maps)
#'
#' # Custom output name and binary file name
#' write_buildmaps_file(maps,
#'                      bin_file = "PrimatesMaps.bin",
#'                      out_file = "BuildPrimatesMaps.txt")
#'
#' # Then in BayesTraits (command line):
#' # BuildMaps BuildPrimatesMaps.txt
#' # LoadMaps  PrimatesMaps.bin
write_buildmaps_file <- function(maps,
                                  bin_file  = "BinaryMap.bin",
                                  prefix    = "mask_",
                                  lon_boxes = 250L,
                                  lat_boxes = 150L,
                                  out_file  = "BuildMaps.txt",
                                  ages      = NULL,
                                  verbose   = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps().")
  if (!is.character(bin_file) || nchar(trimws(bin_file)) == 0)
    stop("'bin_file' must be a non-empty character string.")
  if (grepl(" ", bin_file))
    stop("'bin_file' cannot contain spaces (BayesTraits requirement).")

  # Determine ages to include — skip any that failed to download
  use_ages <- if (is.null(ages)) maps$ages else as.integer(sort(unique(ages)))
  available <- use_ages[!use_ages %in% maps$failed &
                        !vapply(maps$coastlines[as.character(use_ages)],
                                is.null, logical(1))]

  if (length(available) == 0)
    stop("No available ages found. Run get_paleomaps() or fill_missing_maps() first.")

  if (length(available) < length(use_ages)) {
    missing <- setdiff(use_ages, available)
    warning(sprintf(
      "%d age(s) skipped (download failed): %s
  Run fill_missing_maps() to fix.",
      length(missing), paste(missing, collapse = ", ")
    ))
  }

  available <- sort(available)
  n         <- length(available)

  # Build time intervals using midpoints between adjacent ages so that every
  # point in time is covered by exactly one map — works for both consecutive
  # integers (e.g. 0,1,2,3) and irregular spacing (e.g. 0,3,5,10):
  #   First map  : [0, midpoint(age[1], age[2]))
  #   Middle maps: [midpoint(age[i-1], age[i]), midpoint(age[i], age[i+1]))
  #   Last map   : [midpoint(age[n-1], age[n]), age[n] + 0.5)
  midpoints  <- (available[-n] + available[-1]) / 2.0
  start_ages <- c(0.0,       midpoints)
  end_ages   <- c(midpoints, available[n] + 0.5)

  # Build lines
  header <- sprintf("%s\t%d\t%d", bin_file, as.integer(lon_boxes), as.integer(lat_boxes))
  map_lines <- sprintf("%s%d.csv\t%.1f\t%.1f",
                       prefix, available, start_ages, end_ages)

  out <- c(header, map_lines)

  writeLines(out, con = out_file)

  if (verbose) {
    message(sprintf("[write_buildmaps_file] Written: %s", normalizePath(out_file)))
    message(sprintf("  Binary output : %s", bin_file))
    message(sprintf("  Ages included : %d  (%d - %d Ma)",
                    n, min(available), max(available)))
    message(sprintf("  Grid          : %d lon x %d lat boxes",
                    as.integer(lon_boxes), as.integer(lat_boxes)))
    message(sprintf("
To build the binary map file, run in BayesTraits:"))
    message(sprintf("  BuildMaps %s", out_file))
    message(sprintf("Then load it in your analysis with:"))
    message(sprintf("  LoadMaps %s", bin_file))
  }

  invisible(normalizePath(out_file))
}




fill_missing_maps <- function(maps, ages = NULL, verbose = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps().")
  if (!requireNamespace("rgplates", quietly = TRUE))
    stop("Package 'rgplates' is required. Install with: install.packages('rgplates')")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required.")

  # Determine which ages to attempt
  if (is.null(ages)) {
    # NULL entries in the coastlines list
    null_ages <- as.integer(names(Filter(is.null, maps$coastlines)))
    # Ages listed as failed
    fill_ages <- sort(unique(c(null_ages, maps$failed)))
  } else {
    fill_ages <- as.integer(sort(unique(ages)))
    # For user-specified ages not yet in the object, add them
    new_ages <- fill_ages[!as.character(fill_ages) %in% names(maps$coastlines)]
    if (length(new_ages) > 0) {
      maps$coastlines[as.character(new_ages)] <- vector("list", length(new_ages))
      maps$ages   <- sort(unique(c(maps$ages, new_ages)))
      maps$failed <- sort(unique(c(maps$failed, new_ages)))
    }
  }

  if (length(fill_ages) == 0) {
    message("[fill_missing_maps] No missing maps found — nothing to do.")
    return(invisible(maps))
  }

  message(sprintf(
    "[fill_missing_maps] Attempting to download %d missing age(s): %s",
    length(fill_ages), paste(fill_ages, collapse = ", ")
  ))

  newly_fixed  <- integer(0)
  still_failed <- integer(0)

  for (age in fill_ages) {
    key <- as.character(age)
    result <- tryCatch({
      raw <- rgplates::reconstruct("coastlines", age = age, model = maps$model)
      cs <- sf::st_make_valid(raw)
      if (any(sf::st_geometry_type(cs) == "GEOMETRYCOLLECTION")) {
        cs <- tryCatch(
          sf::st_collection_extract(cs, "POLYGON"),
          error = function(e) cs   # already clean — skip extraction
        )
      }
      cs
    }, error = function(e) {
      if (verbose)
        message(sprintf("  Age %d Ma: still failing (%s)", age, conditionMessage(e)))
      NULL
    })

    if (!is.null(result)) {
      maps$coastlines[[key]] <- result
      newly_fixed  <- c(newly_fixed, age)
      if (verbose) message(sprintf("  Age %d Ma: downloaded successfully.", age))
    } else {
      still_failed <- c(still_failed, age)
    }

    Sys.sleep(0.3)   # avoid rate-limiting
  }

  # Update metadata
  maps$failed <- sort(still_failed)
  maps$n_ages <- sum(!vapply(maps$coastlines, is.null, logical(1)))

  if (verbose) {
    if (length(newly_fixed) > 0)
      message(sprintf("[fill_missing_maps] Fixed: %s",
                      paste(newly_fixed, collapse = ", ")))
    if (length(still_failed) > 0)
      message(sprintf("[fill_missing_maps] Still missing: %s — try again later.",
                      paste(still_failed, collapse = ", ")))
    else
      message("[fill_missing_maps] All missing ages resolved.")
  }

  maps
}
