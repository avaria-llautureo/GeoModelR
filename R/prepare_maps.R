# =============================================================================
# GeoModelR - Palaeomap Preparation Functions
# Functions to build geographic input files for BayesTraits Geo model analysis.
# Run these BEFORE your BayesTraits analysis.
#
# get_paleomaps()            Download palaeocoastlines from the GPlates Web Service
# fill_missing_maps()        Re-download any GPlates ages that failed
# load_raster_maps()         Register local raster maps (a time series, or one
#                            present-day file such as ETOPO 2022)
# fetch_etopo2022()          Download and cache the ETOPO 2022 present-day DEM
# paleomap_to_csv()          Convert palaeomaps (GPlates or raster) to binary
#                            1/0 CSV grids for BayesTraits BuildMaps
# sf_to_csv()                Deprecated alias for paleomap_to_csv(source = "gplates")
# export_restriction_maps()  Write restriction boundaries as CSVs (polygon
#                            vertices for GPlates, traced contours for rasters)
# write_buildmaps_file()     Write the BayesTraits BuildMaps input file
#
# Dependencies: rgplates, sf, terra
# =============================================================================

# =============================================================================
# SECTION 1 - PALAEOMAP PREPARATION
# Functions to build geographic input files for the BayesTraits Geo model.
# Run these before your BayesTraits analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers (not exported)
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Repair GPlates coastline polygons and keep only polygon parts.
.clean_coastlines <- function(raw) {
  cs <- sf::st_make_valid(raw)
  if (any(sf::st_geometry_type(cs) == "GEOMETRYCOLLECTION")) {
    cs <- tryCatch(
      sf::st_collection_extract(cs, "POLYGON"),
      error = function(e) cs   # already clean - skip extraction
    )
  }
  cs
}

# Source type of a geo_paleomaps object. Objects made before 0.3.1 have no
# 'type' field; they are GPlates coastlines.
.maps_type <- function(maps) maps$type %||% "coastlines"

# Check mask_range and return it as a list of c(min, max) vectors.
.normalize_mask_range <- function(mask_range) {
  if (is.null(mask_range))
    stop("'mask_range' must be supplied for raster maps, e.g. c(0, Inf) for ",
         "land, c(-200, 0) for a continental-shelf band, or ",
         "list(c(-200, 0), c(1000, Inf)) to combine several bands into one ",
         "1/0 mask.")
  if (is.numeric(mask_range)) mask_range <- list(mask_range)
  ok <- is.list(mask_range) && length(mask_range) > 0 &&
        all(vapply(mask_range, function(r)
          is.numeric(r) && length(r) == 2 && !anyNA(r), logical(1)))
  if (!ok)
    stop("'mask_range' must be a numeric vector c(min, max), or a list of ",
         "such vectors.")
  for (r in mask_range) if (r[1] > r[2])
    stop("Each 'mask_range' entry must be c(min, max) with min <= max.")
  mask_range
}

# Turn a continuous raster into a binary 1/0 mask.
# A cell is 1 when its value falls inside ANY of the ranges (bounds
# inclusive), otherwise 0. 'invert' swaps 1 and 0. Cells with no value (NA)
# always become 0, because BayesTraits reads only 1 and 0.
.apply_mask_range <- function(rast_layer, mask_range, invert = FALSE) {
  ranges <- .normalize_mask_range(mask_range)
  hit <- NULL
  for (r in ranges) {
    inside <- rast_layer >= r[1] & rast_layer <= r[2]
    hit <- if (is.null(hit)) inside else (hit | inside)
  }
  masked <- terra::ifel(hit, 1L, 0L)
  if (invert) masked <- 1L - masked
  masked <- terra::ifel(is.na(masked), 0L, masked)
  names(masked) <- "mask"
  masked
}

# Read one raster file, optionally selecting a variable (NetCDF sub-dataset or
# layer name).
.read_raster <- function(path, values_var = NULL) {
  if (is.null(values_var) || is.na(values_var)) return(terra::rast(path))
  # 'subds' selects a NetCDF variable; other formats (e.g. GeoTIFF) ignore it
  # and return every layer, so then the layer is selected by name.
  r <- tryCatch(terra::rast(path, subds = values_var), error = function(e) NULL)
  if (is.null(r) || terra::nlyr(r) > 1) {
    r <- terra::rast(path)
    if (!values_var %in% names(r))
      stop(sprintf("Variable '%s' not found in '%s'. Available: %s.",
                   values_var, basename(path), paste(names(r), collapse = ", ")))
    r <- r[[values_var]]
  }
  r
}

# Find raster files in 'dir' and take each file's age from its name with
# 'age_regex'. Returns data.frame(age, path) ordered present -> past.
.discover_raster_files <- function(dir, age_regex, pattern = "\\.nc$",
                                   age_group = 1) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0)
    stop(sprintf("No files found in '%s'%s.", dir,
                 if (!is.null(pattern))
                   sprintf(" matching pattern '%s'", pattern) else ""))

  groups  <- regmatches(basename(files),
                        regexec(age_regex, basename(files), perl = TRUE))
  has_hit <- vapply(groups, function(g) length(g) >= age_group + 1, logical(1))

  if (!any(has_hit))
    stop(sprintf("None of the %d file(s) in '%s' matched age_regex '%s'.",
                 length(files), dir, age_regex))

  age_chr <- vapply(groups[has_hit], function(g) g[age_group + 1], character(1))
  ages    <- suppressWarnings(as.numeric(age_chr))
  if (any(is.na(ages)))
    stop("'age_regex' matched non-numeric text for one or more files. ",
         "Check the capture group and 'age_group'.")

  out <- data.frame(age = ages, path = files[has_hit], stringsAsFactors = FALSE)
  out <- out[order(out$age), ]   # present (0 Ma) -> past

  if (any(duplicated(out$age))) {
    dupes <- unique(out$age[duplicated(out$age)])
    stop(sprintf(
      "Several files matched the same age(s): %s. Refine 'pattern' or 'age_regex'.",
      paste(dupes, collapse = ", ")))
  }
  rownames(out) <- NULL
  out
}

# Ages that have usable data in a geo_paleomaps object (either source type).
.paleomaps_available_ages <- function(maps, ages = NULL) {
  if (identical(.maps_type(maps), "raster")) {
    use_ages <- if (is.null(ages)) maps$ages else sort(unique(as.numeric(ages)))
    use_ages[use_ages %in% maps$files$age]
  } else {
    use_ages <- if (is.null(ages)) maps$ages else as.integer(sort(unique(ages)))
    use_ages <- use_ages[as.character(use_ages) %in% names(maps$coastlines)]
    use_ages[!use_ages %in% maps$failed &
             !vapply(maps$coastlines[as.character(use_ages)], is.null, logical(1))]
  }
}

# Warn when any raster that will be written at native resolution is finer
# than 0.5 degrees. Checks each file, so a fine present-day layer (e.g.
# ETOPO 2022 via age0_override) is flagged even when the palaeo series is
# coarse. Only file headers are read.
.warn_fine_resolution <- function(maps, ages) {
  rows <- maps$files[maps$files$age %in% ages, , drop = FALSE]
  fine <- character(0)
  for (i in seq_len(nrow(rows))) {
    r <- tryCatch(.read_raster(rows$path[i], rows$values_var[i]),
                  error = function(e) NULL)
    if (is.null(r)) next
    res <- terra::res(r)
    if (any(res < 0.5))
      fine <- c(fine, sprintf("%g Ma (%s deg, about %s million cells)",
                              rows$age[i],
                              paste(signif(res, 3), collapse = " x "),
                              format(round(terra::ncell(r) / 1e6, 1))))
  }
  if (length(fine) > 0)
    warning(sprintf(paste0(
      "Native resolution is finer than 0.5 degrees for %d map(s): %s.\n",
      "  Output CSVs will be very large and may exhaust memory. Set ",
      "'resolution' (e.g. resolution = 1) to resample."),
      length(fine), paste(fine, collapse = "; ")), call. = FALSE)
  invisible(length(fine) > 0)
}


#' Download and store palaeomaps for a Geo model analysis
#'
#' Fetches palaeocoastlines from the GPlates Web Service for every integer age
#' from 0 Ma to \code{max_age} Ma (1 Ma resolution). All maps are downloaded
#' once and stored in a \code{geo_paleomaps} object. Pass this object to
#' \code{paleomap_to_csv()} to generate BayesTraits input files, or to any
#' plotting function via the \code{maps} argument to avoid repeated GWS calls.
#'
#' For a single present-day map (for example a time-invariant restriction for
#' a young clade), use \code{ages = 0}.
#'
#' For palaeomaps from local raster files (palaeo-elevation and bathymetry
#' such as the Scotese PaleoDEMs, or any other gridded variable), use \code{load_raster_maps()} instead.
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
#' \dontrun{
#' # Requires an internet connection (GPlates Web Service)
#' tree_file <- system.file("extdata", "Median.trees", package = "GeoModelR")
#' anc_file  <- system.file("extdata", "Primate.AncStates.txt",
#'                          package = "GeoModelR")
#'
#' # From a tree (no BayesTraits output needed for data preparation)
#' phy  <- ape::read.nexus(tree_file)
#' maps <- get_paleomaps(max_age = ceiling(max(ape::node.depth.edgelength(phy))))
#'
#' # Or from a geo_ancstates object
#' geo  <- read_geo_ancstates(anc_file)
#' maps <- get_paleomaps(geo)
#'
#' # Use in plotting
#' plot_node(geo, "Node-00000", maps = maps)
#' }
#' @export
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

    # Single-bracket assignment keeps a NULL placeholder for failed ages;
    # coastlines[[i]] <- NULL would delete the element and shift later ages.
    coastlines[i] <- list(tryCatch({
      raw <- rgplates::reconstruct("coastlines", age = age, model = model)
      .clean_coastlines(raw)
    }, error = function(e) {
      message(sprintf("  WARNING: age %d Ma failed (%s) - skipping.",
                      age, conditionMessage(e)))
      failed <<- c(failed, age)
      NULL
    }))

    if (age > 0 && i < n) Sys.sleep(0.3)   # avoid GWS rate-limiting
  }

  if (length(failed) > 0)
    warning(sprintf(
      "%d age(s) failed: %s\n  Re-run get_paleomaps() when the GWS is available.",
      length(failed), paste(failed, collapse = ", ")
    ))

  result <- structure(
    list(type = "coastlines",
         coastlines = coastlines, model = model,
         max_age = as.integer(max_age), ages = ages,
         n_ages = n - length(failed), failed = failed),
    class = "geo_paleomaps"
  )

  message(sprintf("[get_paleomaps] Done. %d/%d maps downloaded.", n - length(failed), n))
  result
}

#' Register local raster maps for a Geo model analysis
#'
#' Registers local raster files as a \code{geo_paleomaps} object, for use with
#' \code{paleomap_to_csv()}, \code{export_restriction_maps()} and
#' \code{write_buildmaps_file()}. Unlike \code{get_paleomaps()}, which
#' downloads coastline polygons from the GPlates Web Service, this function
#' works with any continuous gridded variable in any format \code{terra} can
#' read (NetCDF, GeoTIFF, etc.): palaeo-elevation and bathymetry (for example
#' the Scotese PaleoDEMs), palaeotemperature, or any
#' other raster that you want to turn into a 1/0 restriction with
#' \code{mask_range} in \code{paleomap_to_csv()}.
#'
#' There are two ways to register maps:
#' \itemize{
#'   \item \strong{A time series} (\code{dir} and \code{age_regex}): one file
#'         per age in a folder. Each file's age is read from its file name
#'         with \code{age_regex}, so no naming convention is assumed. Files are
#'         always ordered from the present to the past.
#'   \item \strong{One present-day map} (\code{file}): a single raster
#'         registered as the 0 Ma map, for example ETOPO 2022 from
#'         \code{fetch_etopo2022()}. Use this for a time-invariant restriction,
#'         e.g. a young clade or a within-species tree that evolved on a
#'         practically unchanged Earth. Give \code{write_buildmaps_file()} a
#'         \code{max_age} so that this one map covers every node.
#' }
#'
#' Rasters must be in geographic coordinates (EPSG:4326, longitude and
#' latitude in degrees); reproject with \code{terra::project()} first if they
#' are not. Rasters are not loaded into memory here: only the file paths are
#' recorded (\code{$files}), and each file is read when it is needed.
#'
#' @param dir        Folder containing the raster files, one file per age.
#'                   Use with \code{age_regex}. Leave \code{NULL} when using
#'                   \code{file}.
#' @param age_regex  Perl-compatible regular expression with one capture group
#'                   around the numeric age in the file name. Examples:
#'                   \code{"_([0-9]+)Ma\\\\.nc$"} for Scotese-style names such
#'                   as \code{"Map07_PALEOMAP_6min_Early_Miocene_20Ma.nc"}, and
#'                   \code{"^([0-9]+)_"} for names that start with the age,
#'                   such as \code{"100_elevation.tif"}.
#' @param file       Path to one raster file to register as the present-day
#'                   (0 Ma) map. Use instead of \code{dir} and
#'                   \code{age_regex}.
#' @param age_group  Which capture group in \code{age_regex} holds the age
#'                   (default 1).
#' @param pattern    File name filter passed to \code{list.files()} before
#'                   \code{age_regex} is applied (default \code{"\\\\.nc$"}).
#'                   Use e.g. \code{"\\\\.tif$"} for GeoTIFF files, or
#'                   \code{NULL} for every file in \code{dir}.
#' @param values_var Character or \code{NULL}. Name of the variable (NetCDF
#'                   sub-dataset or layer) that holds the values, e.g.
#'                   \code{"z"} (Scotese PaleoDEMs). Required when a file contains more than one
#'                   variable. If \code{NULL} and a file has exactly one
#'                   layer, that layer is used.
#' @param ages       Numeric vector or \code{NULL}. Subset of ages to register
#'                   from \code{dir} (e.g. a time span). \code{NULL} (default)
#'                   registers every age found.
#' @param label      Character or \code{NULL}. Short name for this map source
#'                   (e.g. \code{"Scotese PaleoDEM"}), used only when printing.
#' @param check_crs  Logical. Stop if a raster is not in geographic lon/lat
#'                   coordinates (default \code{TRUE}). Set \code{FALSE} only
#'                   if you are certain the coordinates are lon/lat degrees
#'                   despite a missing or unrecognised CRS.
#' @param age0_override Character path, \code{list(path = , values_var = )},
#'                   or \code{NULL} (default). Only with \code{dir}.
#'                   Palaeo-elevation models are reconstructions even at
#'                   0 Ma; an observational DEM such as ETOPO 2022 (see
#'                   \code{fetch_etopo2022()}) is the real present day. When
#'                   supplied, this raster replaces the series' own 0 Ma map
#'                   (or is added as 0 Ma if the series has none). It must use
#'                   the same sign convention as the series (positive = above
#'                   present sea level), because the same \code{mask_range}
#'                   is applied to it. It is usually much finer than the
#'                   palaeo series, so always set \code{resolution} in
#'                   \code{paleomap_to_csv()} when using it.
#' @param verbose    Logical. Print a summary (default \code{TRUE}).
#' @return A \code{geo_paleomaps} object with \code{type = "raster"}.
#'
#' @seealso \code{\link{paleomap_to_csv}()}, \code{\link{fetch_etopo2022}()},
#'   \code{\link{write_buildmaps_file}()}
#'
#' @examples
#' \dontrun{
#' # The raster files are not part of the package: download them first.
#'
#' # Scotese PaleoDEMs, one file every 5 Ma, only 0-100 Ma
#' scotese <- load_raster_maps(dir = "scotese_nc/",
#'                             age_regex = "_([0-9]+)Ma\\.nc$",
#'                             values_var = "z", ages = seq(0, 100, 5),
#'                             label = "Scotese PaleoDEM")
#'
#' # The same, with the real present-day DEM (ETOPO 2022) replacing the
#' # series' own 0 Ma map
#' etopo    <- fetch_etopo2022()
#' scotese2 <- load_raster_maps(dir = "scotese_nc/",
#'                              age_regex = "_([0-9]+)Ma\\.nc$",
#'                              values_var = "z", age0_override = etopo)
#' scotese2$files     # check the age-0 row
#'
#' # One present-day map only (time-invariant restriction)
#' present <- load_raster_maps(file = etopo, label = "ETOPO 2022")
#' }
#' @export
load_raster_maps <- function(dir           = NULL,
                             age_regex     = NULL,
                             file          = NULL,
                             age_group     = 1,
                             pattern       = "\\.nc$",
                             values_var    = NULL,
                             ages          = NULL,
                             label         = NULL,
                             check_crs     = TRUE,
                             age0_override = NULL,
                             verbose       = TRUE) {

  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install with: install.packages('terra')")

  use_file <- !is.null(file)
  if (use_file && (!is.null(dir) || !is.null(age_regex)))
    stop("Supply either 'file' (one present-day map) or 'dir' and ",
         "'age_regex' (a time series), not both.")
  if (!use_file && (is.null(dir) || is.null(age_regex)))
    stop("Supply 'dir' and 'age_regex' for a time series, or 'file' for one ",
         "present-day map.")

  if (use_file) {
    if (!is.character(file) || length(file) != 1)
      stop("'file' must be a single file path.")
    if (!file.exists(file))
      stop(sprintf("File not found: %s", file))
    if (!is.null(ages))
      stop("'ages' cannot be used with 'file': the file is registered as 0 Ma.")
    if (!is.null(age0_override))
      stop("'age0_override' cannot be used with 'file': 'file' is already ",
           "the present-day map.")
    file_index <- data.frame(age = 0, path = file, stringsAsFactors = FALSE)
  } else {
    if (!dir.exists(dir))
      stop(sprintf("Folder not found: %s", dir))
    file_index <- .discover_raster_files(dir, age_regex = age_regex,
                                         pattern = pattern,
                                         age_group = age_group)
    if (!is.null(ages)) {
      missing <- setdiff(ages, file_index$age)
      if (length(missing) > 0)
        warning(sprintf("No file found for requested age(s): %s",
                        paste(missing, collapse = ", ")))
      file_index <- file_index[file_index$age %in% ages, , drop = FALSE]
    }
    if (nrow(file_index) == 0)
      stop("No matching raster files found for the requested ages.")
  }

  .check_lonlat <- function(r, what) {
    if (check_crs && !isTRUE(suppressWarnings(terra::is.lonlat(r))))
      stop(sprintf(paste0(
        "%s is not in geographic lon/lat coordinates (EPSG:4326). ",
        "Reproject with terra::project(r, 'EPSG:4326') first, or set ",
        "check_crs = FALSE if you are certain the coordinates are already ",
        "lon/lat degrees."), what))
  }

  # Read one file's header, check it has one layer and lon/lat coordinates,
  # and return its variable name and resolution.
  .probe <- function(path, var, what, hint) {
    r <- tryCatch(.read_raster(path, var), error = function(e)
      stop(sprintf("Could not read %s: %s", what, conditionMessage(e)),
           call. = FALSE))
    if (terra::nlyr(r) > 1)
      stop(sprintf("%s contains %d layers (%s). %s", what, terra::nlyr(r),
                   paste(names(r), collapse = ", "), hint), call. = FALSE)
    .check_lonlat(r, what)
    list(var = if (is.null(var)) names(r)[1] else var, res = terra::res(r))
  }

  first <- .probe(file_index$path[1], values_var,
                  sprintf("'%s'", basename(file_index$path[1])),
                  "Set 'values_var' to the variable to use.")
  file_index$values_var <- first$var

  if (!is.null(age0_override)) {
    ov_path <- if (is.list(age0_override)) age0_override$path else age0_override
    ov_var  <- if (is.list(age0_override)) age0_override$values_var else NULL
    if (!is.character(ov_path) || length(ov_path) != 1 || !file.exists(ov_path))
      stop("'age0_override' must be an existing file path, or ",
           "list(path = , values_var = ).")
    ov <- .probe(ov_path, ov_var,
                 sprintf("age0_override '%s'", basename(ov_path)),
                 "Supply age0_override = list(path = , values_var = ).")

    if (0 %in% file_index$age) {
      file_index$path[file_index$age == 0]       <- ov_path
      file_index$values_var[file_index$age == 0] <- ov$var
    } else {
      file_index <- rbind(
        data.frame(age = 0, path = ov_path, values_var = ov$var,
                   stringsAsFactors = FALSE),
        file_index)
    }
    rownames(file_index) <- NULL
    if (verbose)
      message(sprintf(paste0(
        "[load_raster_maps] Age 0 Ma uses '%s' (variable '%s') instead of ",
        "the series' own 0 Ma map. Set 'resolution' in paleomap_to_csv()."),
        basename(ov_path), ov$var))
  }

  if (verbose)
    message(sprintf(paste0(
      "[load_raster_maps] Registered %d map(s), %g-%g Ma%s.\n",
      "  Variable: '%s'  |  Native resolution: %s deg."),
      nrow(file_index), min(file_index$age), max(file_index$age),
      if (!is.null(label)) sprintf(" (%s)", label) else "",
      first$var, paste(signif(first$res, 4), collapse = " x ")))

  structure(
    list(type       = "raster",
         files      = file_index,
         values_var = first$var,
         native_res = first$res,
         model      = if (is.null(label)) "local raster" else label,
         max_age    = max(file_index$age),
         ages       = file_index$age,
         n_ages     = nrow(file_index),
         failed     = numeric(0)),
    class = "geo_paleomaps"
  )
}

#' Download and cache the ETOPO 2022 present-day global relief model
#'
#' Downloads the ETOPO 2022 60 arc-second bedrock global relief model from
#' NOAA's National Centers for Environmental Information (NCEI). It is one
#' global NetCDF file of land elevation and ocean depth (positive values above
#' present sea level, negative below), at about 1.85 km resolution at the
#' equator. Use it as a present-day map in \code{load_raster_maps()}, either
#' alone (\code{file = }) or in place of a palaeo series' own 0 Ma map
#' (\code{age0_override = }).
#'
#' The file is large. It is cached at \code{dest} and reused on later calls
#' unless \code{force = TRUE}. If a download fails, the partial file is
#' deleted so it is not reused.
#'
#' ETOPO 2022 is released by NOAA NCEI under a CC0 public domain dedication.
#' Please cite it when you use it (see References).
#'
#' @param dest    File path for the cached NetCDF. By default a folder for
#'                this package in the user's cache directory
#'                (\code{tools::R_user_dir()}).
#' @param force   Logical. Download again even if \code{dest} exists
#'                (default \code{FALSE}).
#' @param timeout Numeric, seconds. Download time limit (default 1800, i.e.
#'                30 minutes). R's own default of 60 seconds is too short for
#'                this file. The setting is restored after the call.
#' @return The local file path, invisibly. It can be passed directly to
#'         \code{file} or \code{age0_override} in \code{load_raster_maps()}.
#'
#' @references
#' MacFerrin, M., Amante, C., Carignan, K., Love, M., & Lim, E. (2025). The
#' Earth Topography 2022 (ETOPO 2022) global DEM dataset. \emph{Earth System
#' Science Data}, 17, 1835-1849. \doi{10.5194/essd-17-1835-2025}
#'
#' NOAA National Centers for Environmental Information (2022). ETOPO 2022
#' 15 Arc-Second Global Relief Model. \doi{10.25921/fd45-gt74}. (NCEI's
#' recommended dataset citation. This DOI record is titled for the 15
#' arc-second version; the 30 and 60 arc-second files are listed in the same
#' record.)
#'
#' @seealso \code{\link{load_raster_maps}()}
#'
#' @examples
#' \dontrun{
#' # Requires an internet connection (large download)
#' etopo   <- fetch_etopo2022()
#' present <- load_raster_maps(file = etopo, label = "ETOPO 2022")
#' }
#' @export
fetch_etopo2022 <- function(dest    = file.path(
                                        tools::R_user_dir(utils::packageName(),
                                                          "cache"),
                                        "ETOPO_2022_v1_60s_N90W180_bed.nc"),
                            force   = FALSE,
                            timeout = 1800) {
  url <- paste0(
    "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/",
    "60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc")

  if (!force && file.exists(dest)) {
    message(sprintf("[fetch_etopo2022] Using cached file: %s", dest))
    return(invisible(dest))
  }
  if (!dir.exists(dirname(dest))) dir.create(dirname(dest), recursive = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  message(sprintf(paste0(
    "[fetch_etopo2022] Downloading ETOPO 2022 (60 arc-second bedrock, ",
    "large file, timeout = %d s)..."), as.integer(timeout)))

  status <- tryCatch(
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
    error = function(e) {
      if (file.exists(dest)) file.remove(dest)   # never cache a partial file
      stop(sprintf(paste0(
        "ETOPO 2022 download failed (%s). If it timed out, try a larger ",
        "'timeout'."), conditionMessage(e)), call. = FALSE)
    })
  if (!identical(as.integer(status), 0L)) {
    if (file.exists(dest)) file.remove(dest)
    stop(sprintf("ETOPO 2022 download failed (status %s).", status),
         call. = FALSE)
  }

  message(sprintf("[fetch_etopo2022] Saved to: %s", dest))
  invisible(dest)
}

#' Print method for geo_paleomaps
#' @param x   A \code{geo_paleomaps} object from \code{get_paleomaps()} or
#'            \code{load_raster_maps()}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.geo_paleomaps <- function(x, ...) {
  raster <- identical(.maps_type(x), "raster")
  cat("Palaeomap collection (geo_paleomaps)\n")
  cat("  Source :", if (raster) "local raster" else "GPlates coastlines", "\n")
  cat("  Model  :", x$model, "\n")
  cat("  Range  : 0 -", x$max_age, "Ma\n")
  cat("  Maps   :", x$n_ages, if (raster) "ages registered\n" else "ages downloaded\n")
  if (raster) {
    cat("  Variable   :", x$values_var, "\n")
    cat("  Native res :", paste(signif(x$native_res, 4), collapse = " x "),
        "deg\n")
  }
  if (length(x$failed) > 0)
    cat("  Failed :", paste(x$failed, collapse = ", "), "Ma\n")
  if (raster)
    cat("\nPass to paleomap_to_csv(), export_restriction_maps() or",
        "write_buildmaps_file().\n")
  else
    cat("\nPass to plotting functions via maps = or to paleomap_to_csv().\n")
  invisible(x)
}

#' Convert palaeomaps to binary 1/0 CSV files for BayesTraits
#'
#' Turns palaeomaps into the binary grid CSV format required by the
#' BayesTraits \code{BuildMaps} command. Two sources are supported:
#'
#' \itemize{
#'   \item \code{source = "gplates"} (default): palaeocoastline polygons, from
#'         \code{get_paleomaps()} or downloaded on the fly, are rasterised.
#'         Cells inside a polygon become 1 (land), all others 0 (ocean).
#'   \item \code{source = "raster"}: a continuous raster from
#'         \code{load_raster_maps()} (palaeo-elevation and bathymetry, or any
#'         other gridded variable) is thresholded with \code{mask_range}.
#'         Cells whose value falls inside \code{mask_range} become 1, all
#'         others 0. Because the source is a continuous surface,
#'         \code{mask_range} can define land, ocean, a shoreline or shelf
#'         band, a depth or elevation zone, etc. \code{invert} swaps 1 and 0.
#'         Cells with no value (\code{NA}) always become 0, because
#'         BayesTraits reads only 1 and 0.
#' }
#'
#' For a single age, returns a data frame for inspection. For several ages
#' (or \code{write = TRUE}), writes one file per age to \code{dir}, named
#' \code{<prefix><age>.csv}.
#'
#' @section Resolution:
#' \code{resolution = 1} gives 1 x 1 degree cells: about 111 km at the
#' equator, and narrower in longitude towards the poles. For raster sources
#' the native resolution is used by default; a warning is given when any map
#' to be written is finer than 0.5 degrees (e.g. the Scotese PaleoDEMs at
#' 0.1 degrees, or ETOPO 2022 at 1 arc-minute, about 233 million cells),
#' because the CSVs become very large. \code{resolution} applies to every
#' age, including a present-day map from \code{age0_override}.
#'
#' Setting \code{resolution} finer than a coarse source (e.g.
#' \code{resolution = 0.25} for a 1 degree source) interpolates the
#' values before thresholding, giving a smoother 1/0 boundary instead of
#' large blocky plateaus, which may help the MCMC move. This is interpolation
#' of the existing data, not new information, and should be described that
#' way in any write-up.
#'
#' @param age        Numeric vector. Age(s) in Ma. Default \code{0}.
#' @param maps       A \code{geo_paleomaps} object: from
#'                   \code{get_paleomaps()} for \code{source = "gplates"}
#'                   (optional, recommended for batch use), or from
#'                   \code{load_raster_maps()} for \code{source = "raster"}
#'                   (required).
#' @param source     Character. \code{"gplates"} (default) or \code{"raster"}.
#' @param model      Character. rgplates model (default \code{"PALEOMAP"}).
#'                   Only used for \code{source = "gplates"} when \code{maps}
#'                   is not supplied.
#' @param mask_range Raster sources only. A numeric \code{c(min, max)}, or a
#'                   list of such vectors to combine several bands into one
#'                   mask. Bounds are inclusive. There is no default: you
#'                   decide what 1 means for your analysis.
#' @param invert     Raster sources only. Logical; swaps 1 and 0 (default
#'                   \code{FALSE}). \code{NA} cells stay 0.
#' @param extent     Numeric \code{c(xmin, xmax, ymin, ymax)}.
#'                   Default: global \code{c(-180, 180, -90, 90)}.
#' @param resolution Numeric or \code{NULL}. Grid cell size in degrees.
#'                   \code{NULL} (default) means 1 degree for GPlates sources
#'                   and the native resolution for raster sources. See the
#'                   Resolution section.
#' @param resample_method Raster sources only. \code{"bilinear"} (default) or
#'                   \code{"cubicspline"}, passed to \code{terra::resample()}
#'                   when \code{resolution} is set. \code{"cubicspline"} gives
#'                   a smoother surface at extra compute cost.
#' @param dir        Character. Output folder. Default \code{"maps_csv"}.
#' @param write      Logical. Write a file even for a single age (default
#'                   \code{FALSE}).
#' @param prefix     Character. File name prefix (default \code{"mask_"}).
#' @param verbose    Logical. Print progress (default \code{TRUE}).
#' @return Single age with \code{write = FALSE}: a data frame
#'         (\code{Lon}, \code{Lat}, \code{mask}). Otherwise: the output
#'         folder path, invisibly.
#'
#' @seealso \code{\link{get_paleomaps}()}, \code{\link{load_raster_maps}()},
#'   \code{\link{write_buildmaps_file}()}
#'
#' @examples
#' \dontrun{
#' # --- GPlates coastlines (requires the GPlates Web Service) ---
#' df <- paleomap_to_csv(age = 50)
#' library(ggplot2)
#' ggplot(df, aes(x = Lon, y = Lat, fill = factor(mask))) +
#'   geom_raster() +
#'   scale_fill_manual(values = c("0" = "#b8d4e8", "1" = "#d2c8a0"),
#'                     labels = c("Ocean", "Land"), name = NULL) +
#'   coord_fixed(expand = FALSE) + theme_bw()
#'
#' maps <- get_paleomaps(max_age = 74)
#' paleomap_to_csv(age = maps$ages, maps = maps, dir = "maps_csv")
#'
#' # --- Local raster maps (files must be downloaded first) ---
#' scotese <- load_raster_maps(dir = "scotese_nc/",
#'                             age_regex = "_([0-9]+)Ma\\.nc$",
#'                             values_var = "z",
#'                             age0_override = fetch_etopo2022())
#'
#' # Check one age: land (elevation >= 0 m)
#' df <- paleomap_to_csv(age = 100, maps = scotese, source = "raster",
#'                       mask_range = c(0, Inf), resolution = 1)
#'
#' # All ages: shelf band (-200 to 0 m), smoother boundary
#' paleomap_to_csv(age = scotese$ages, maps = scotese, source = "raster",
#'                 mask_range = c(-200, 0), resolution = 0.5,
#'                 resample_method = "cubicspline", dir = "maps_csv_shelf")
#' }
#' @export
paleomap_to_csv <- function(age             = 0,
                            maps            = NULL,
                            source          = c("gplates", "raster"),
                            model           = "PALEOMAP",
                            mask_range      = NULL,
                            invert          = FALSE,
                            extent          = c(-180, 180, -90, 90),
                            resolution      = NULL,
                            resample_method = c("bilinear", "cubicspline"),
                            dir             = "maps_csv",
                            write           = FALSE,
                            prefix          = "mask_",
                            verbose         = TRUE) {

  source          <- match.arg(source)
  resample_method <- match.arg(resample_method)

  if (!is.numeric(age) || length(age) == 0 || anyNA(age) || any(age < 0))
    stop("'age' must be a non-negative numeric vector.")
  if (length(extent) != 4 || !is.numeric(extent))
    stop("'extent' must be c(xmin, xmax, ymin, ymax).")
  if (extent[1] < -180 || extent[2] > 180 || extent[3] < -90 || extent[4] > 90)
    stop("'extent' values outside valid range: lon [-180, 180], lat [-90, 90].")
  if (!is.null(resolution) &&
      (!is.numeric(resolution) || length(resolution) != 1 || resolution <= 0))
    stop("'resolution' must be a positive number.")
  if (!is.logical(invert) || length(invert) != 1 || is.na(invert))
    stop("'invert' must be TRUE or FALSE.")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install with: install.packages('terra')")
  if (!is.null(maps) && !inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps() or ",
         "load_raster_maps().")
  if (!is.null(maps)) {
    maps_source <- if (identical(.maps_type(maps), "raster")) "raster" else "gplates"
    if (!identical(maps_source, source))
      stop(sprintf(paste0(
        "source = '%s' was requested, but 'maps' comes from %s. ",
        "Use source = '%s'."), source,
        if (maps_source == "raster") "load_raster_maps()" else "get_paleomaps()",
        maps_source))
  }

  age       <- sort(unique(age))
  multi_age <- length(age) > 1 || write
  if (multi_age && !dir.exists(dir)) dir.create(dir, recursive = TRUE)

  results <- vector("list", length(age))
  names(results) <- as.character(age)

  # Shared final step: data frame in BayesTraits coordinate limits, then
  # write it or keep it.
  .finish <- function(map_rast, a, fmt) {
    df <- as.data.frame(map_rast, xy = TRUE, na.rm = FALSE)
    colnames(df) <- c("Lon", "Lat", "mask")
    df$mask <- as.integer(df$mask)
    df <- df[df$Lon >= -180 & df$Lon < 180 &
             df$Lat >  -90  & df$Lat <  90, ]
    rownames(df) <- NULL
    if (multi_age) {
      utils::write.csv(df, file = file.path(dir, sprintf(fmt, prefix, a)),
                       row.names = FALSE, quote = FALSE)
      NULL
    } else df
  }

  if (source == "raster") {
    # -------------------------------------------------------------------
    # RASTER SOURCE (load_raster_maps)
    # -------------------------------------------------------------------
    if (is.null(maps))
      stop("For source = 'raster', supply a 'maps' object from load_raster_maps().")
    mask_range <- .normalize_mask_range(mask_range)
    req_ext    <- terra::ext(extent[1], extent[2], extent[3], extent[4])

    if (is.null(resolution)) .warn_fine_resolution(maps, age)

    if (verbose && multi_age)
      message(sprintf("[paleomap_to_csv] Writing %d map(s) to: %s/",
                      length(age), dir))

    for (i in seq_along(age)) {
      a <- age[i]
      if (verbose && multi_age && (i %% 10 == 0 || i == length(age)))
        message(sprintf("  ... %d / %d", i, length(age)))

      row <- maps$files[maps$files$age == a, , drop = FALSE]
      if (nrow(row) == 0) {
        warning(sprintf("Age %g Ma not found in the maps object - skipping.", a))
        next
      }

      r <- tryCatch(.read_raster(row$path[1], row$values_var[1]),
                    error = function(e) {
                      warning(sprintf("Age %g Ma: could not read raster (%s) - skipping.",
                                      a, conditionMessage(e)))
                      NULL
                    })
      if (is.null(r)) next
      if (terra::nlyr(r) > 1) {
        warning(sprintf("Age %g Ma: raster has %d layers; expected 1 - skipping.",
                        a, terra::nlyr(r)))
        next
      }

      r <- terra::crop(r, req_ext)
      if (!is.null(resolution)) {
        tmpl <- terra::rast(extent = req_ext, resolution = resolution,
                            crs = "EPSG:4326")
        r <- terra::resample(r, tmpl, method = resample_method)
      }

      masked <- .apply_mask_range(r, mask_range, invert = invert)
      res_i  <- .finish(masked, a, "%s%g.csv")
      if (!multi_age) results[[1]] <- res_i
    }

  } else {
    # -------------------------------------------------------------------
    # GPLATES SOURCE (get_paleomaps, or on-the-fly download)
    # -------------------------------------------------------------------
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required.")
    if (!is.null(mask_range) || isTRUE(invert))
      warning("'mask_range' and 'invert' are only used for raster maps ",
              "(source = 'raster'); ignored.")

    age <- sort(unique(as.integer(round(age))))
    if (is.null(resolution)) resolution <- 1

    # Raster template built once and reused for all ages
    grid_template <- terra::rast(
      xmin = extent[1], xmax = extent[2],
      ymin = extent[3], ymax = extent[4],
      res  = resolution, crs = "EPSG:4326"
    )

    if (verbose && multi_age)
      message(sprintf("[paleomap_to_csv] Writing %d map(s) at %.2f deg resolution to: %s/",
                      length(age), resolution, dir))

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
          .clean_coastlines(
            rgplates::reconstruct("coastlines", age = a, model = model))
        }
      }, error = function(e) {
        warning(sprintf("Age %d Ma: could not get coastlines (%s) - skipping.",
                        a, conditionMessage(e)))
        NULL
      })
      if (is.null(coast_sf)) next

      # terra::rasterize requires a plain sf - strip any extra classes
      class(coast_sf) <- c("sf", "data.frame")

      map_rast <- tryCatch(
        terra::rasterize(coast_sf, grid_template),
        error = function(e) {
          warning(sprintf("Age %d Ma: rasterization failed (%s) - skipping.",
                          a, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(map_rast)) next

      # Binary: land = 1, ocean = 0
      vals <- terra::values(map_rast)
      terra::values(map_rast) <- ifelse(is.na(vals), 0L, 1L)

      res_i <- .finish(map_rast, a, "%s%d.csv")
      if (!multi_age) results[[1]] <- res_i
    }
  }

  if (multi_age) {
    if (verbose)
      message(sprintf("[paleomap_to_csv] Done. Files written to: %s/",
                      normalizePath(dir)))
    invisible(normalizePath(dir))
  } else {
    results[[1]]
  }
}

#' Convert GPlates palaeomaps to binary CSV files (deprecated)
#'
#' \code{sf_to_csv()} is deprecated. Use \code{paleomap_to_csv()}, which works
#' with both GPlates coastlines and local raster maps. \code{sf_to_csv()}
#' keeps its arguments and defaults so existing code still runs; it calls
#' \code{paleomap_to_csv(..., source = "gplates")}.
#'
#' @inheritParams paleomap_to_csv
#' @param resolution Numeric. Grid cell size in degrees (default \code{1}).
#' @return As \code{paleomap_to_csv()}.
#' @seealso \code{\link{paleomap_to_csv}()}
#' @keywords internal
#' @export
sf_to_csv <- function(age        = 0,
                      maps       = NULL,
                      model      = "PALEOMAP",
                      extent     = c(-180, 180, -90, 90),
                      resolution = 1,
                      dir        = "maps_csv",
                      write      = FALSE,
                      prefix     = "mask_",
                      verbose    = TRUE) {
  .Deprecated("paleomap_to_csv", package = utils::packageName())
  if (!is.numeric(resolution) || length(resolution) != 1 || resolution <= 0)
    stop("'resolution' must be a positive number.")
  paleomap_to_csv(age = age, maps = maps, source = "gplates", model = model,
                  extent = extent, resolution = resolution, dir = dir,
                  write = write, prefix = prefix, verbose = verbose)
}

#' Export palaeomap restriction boundaries as CSVs
#'
#' Writes one CSV file per age with the coordinates of the restriction
#' boundary, rather than the rasterised grid. Each file is named
#' \code{mask_N.csv}, where N is the age in Ma. The output depends on the
#' source of \code{maps}:
#' \itemize{
#'   \item GPlates maps (\code{get_paleomaps()}): the coastline polygon
#'         vertex coordinates.
#'   \item Raster maps (\code{load_raster_maps()}): the \code{mask_range}
#'         boundary traced as a contour line (e.g. the shoreline, or the
#'         edge of an elevation or depth band), at the raster's native
#'         resolution, and its vertex coordinates.
#' }
#'
#' Note: \code{paleomap_to_csv()} produces rasterised binary grids
#' (Lon/Lat/mask), which is the format required by BayesTraits
#' \code{BuildMaps}. This function exports the boundary geometry instead,
#' which may be useful for custom processing or visualisation.
#'
#' @param maps       A \code{geo_paleomaps} object from \code{get_paleomaps()}
#'                   or \code{load_raster_maps()}.
#' @param dir        Output folder. Default \code{"restriction_maps"}.
#' @param ages       Numeric vector or \code{NULL}. Subset of ages.
#'                   \code{NULL} = all ages.
#' @param mask_range Raster maps only, and required for them. See
#'                   \code{paleomap_to_csv()}.
#' @param invert     Raster maps only. See \code{paleomap_to_csv()}.
#' @param verbose    Logical. Print progress (default \code{TRUE}).
#' @return The output folder path, invisibly.
#'
#' @examples
#' \dontrun{
#' # Requires an internet connection (GPlates Web Service)
#' geo  <- read_geo_ancstates(system.file("extdata", "Primate.AncStates.txt",
#'                                        package = "GeoModelR"))
#' maps <- get_paleomaps(geo)
#' export_restriction_maps(maps, dir = tempdir())
#' export_restriction_maps(maps, dir = tempdir(), ages = get_node_ages(geo))
#'
#' # Raster maps: traced shoreline (raster files must be downloaded first)
#' scotese <- load_raster_maps(dir = "scotese_nc/",
#'                             age_regex = "_([0-9]+)Ma\\.nc$",
#'                             values_var = "z")
#' export_restriction_maps(scotese, dir = tempdir(), mask_range = c(0, Inf))
#' }
#' @export
export_restriction_maps <- function(maps,
                                    dir        = "restriction_maps",
                                    ages       = NULL,
                                    mask_range = NULL,
                                    invert     = FALSE,
                                    verbose    = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps() or ",
         "load_raster_maps().")
  raster <- identical(.maps_type(maps), "raster")
  if (raster) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("Package 'terra' is required.")
    mask_range <- .normalize_mask_range(mask_range)
  } else {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required.")
    if (!is.null(mask_range) || isTRUE(invert))
      warning("'mask_range' and 'invert' are only used for raster maps; ignored.")
  }

  export_ages <- if (is.null(ages)) maps$ages
                 else if (raster) as.numeric(ages) else as.integer(ages)
  export_ages <- export_ages[export_ages %in% maps$ages]
  if (length(export_ages) == 0) stop("No matching ages found in the maps object.")

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  n_written <- 0L

  for (age in export_ages) {
    if (raster) {
      row <- maps$files[maps$files$age == age, , drop = FALSE]
      masked <- .apply_mask_range(.read_raster(row$path[1], row$values_var[1]),
                                  mask_range, invert = invert)
      contour <- tryCatch(terra::as.contour(masked, levels = 0.5),
                          error = function(e) NULL)
      if (is.null(contour) || nrow(contour) == 0) {
        if (verbose)
          message(sprintf("  Age %g Ma: no boundary at this mask_range - skipping.", age))
        next
      }
      all_coords <- as.data.frame(terra::crds(contour))
      names(all_coords) <- c("longitude", "latitude")
      fname <- sprintf("mask_%g.csv", age)
    } else {
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
      fname <- sprintf("mask_%d.csv", age)
    }

    utils::write.table(all_coords, file = file.path(dir, fname),
                       sep = ",", row.names = FALSE, col.names = FALSE,
                       quote = FALSE)
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
#' is assigned the time interval closest to its age, using midpoints between
#' neighbouring maps, so every point in time is covered by exactly one map.
#' With maps every 1 Ma:
#'
#' \itemize{
#'   \item \code{mask_0.csv}   covers \code{[0, 0.5)} Ma
#'   \item \code{mask_1.csv}   covers \code{[0.5, 1.5)} Ma
#'   \item \code{mask_N.csv}   covers \code{[N-0.5, N+0.5)} Ma  (middle maps)
#'   \item \code{mask_MAX.csv} covers \code{[MAX-0.5, MAX+0.5)} Ma
#' }
#'
#' Irregular spacing works the same way (e.g. Scotese maps every 5 Ma:
#' \code{mask_5.csv} covers \code{[2.5, 7.5)} Ma).
#'
#' By default the last map ends 0.5 Ma after its age. Set \code{max_age}
#' (e.g. the root age of the tree) to extend the last map so every node gets
#' a restriction map. This is needed when the oldest map is younger than the
#' root, and in particular for a single present-day map used as a
#' time-invariant restriction (\code{get_paleomaps(ages = 0)} or
#' \code{load_raster_maps(file = )}).
#'
#' The output file can be passed directly to BayesTraits with the
#' \code{BuildMaps} command. After \code{BuildMaps} has run, load the resulting
#' binary file in BayesTraits with \code{LoadMaps BinaryMap.bin}.
#'
#' @param maps       A \code{geo_paleomaps} object from \code{get_paleomaps()}
#'                   or \code{load_raster_maps()}. Used to determine the set of
#'                   available ages.
#' @param bin_file   Character. Name of the binary map file that BayesTraits
#'                   will create (default \code{"BinaryMap.bin"}).
#' @param prefix     Character. CSV file name prefix (default \code{"mask_"}).
#'                   Must match the prefix used in \code{paleomap_to_csv()}.
#' @param lon_boxes  Integer. Number of longitude grid boxes for BayesTraits
#'                   spatial indexing (default 250, as recommended in the manual).
#' @param lat_boxes  Integer. Number of latitude grid boxes (default 150).
#' @param out_file   Character. Output file path (default \code{"BuildMaps.txt"}).
#' @param ages       Numeric vector or \code{NULL}. Subset of ages to include.
#'                   \code{NULL} (default) uses all ages in the maps object.
#' @param max_age    Numeric or \code{NULL}. Oldest age (Ma) that must be
#'                   covered, e.g. the root age of the tree. When set, the
#'                   last map covers up to \code{ceiling(max_age) + 0.5} Ma.
#'                   Must not be younger than the oldest map. \code{NULL}
#'                   (default) keeps the last map's end at its age + 0.5 Ma.
#' @param verbose    Logical. Print a summary on completion (default TRUE).
#' @return Invisibly returns the output file path.
#'
#' @examples
#' \dontrun{
#' # Requires an internet connection (GPlates Web Service)
#' maps <- get_paleomaps(max_age = 74)
#' paleomap_to_csv(age = maps$ages, maps = maps, dir = "maps_csv")
#'
#' # Write the BuildMaps file
#' write_buildmaps_file(maps)
#'
#' # Custom output name and binary file name
#' write_buildmaps_file(maps,
#'                      bin_file = "PrimatesMaps.bin",
#'                      out_file = "BuildPrimatesMaps.txt")
#'
#' # One present-day map for the whole tree (time-invariant restriction)
#' present <- get_paleomaps(ages = 0)
#' paleomap_to_csv(age = 0, maps = present, write = TRUE, dir = "maps_csv")
#' write_buildmaps_file(present, max_age = 1.2)   # root age of the tree
#'
#' # Then in BayesTraits (command line):
#' # BuildMaps BuildPrimatesMaps.txt
#' # LoadMaps  PrimatesMaps.bin
#' }
#' @export
write_buildmaps_file <- function(maps,
                                 bin_file  = "BinaryMap.bin",
                                 prefix    = "mask_",
                                 lon_boxes = 250L,
                                 lat_boxes = 150L,
                                 out_file  = "BuildMaps.txt",
                                 ages      = NULL,
                                 max_age   = NULL,
                                 verbose   = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps() or ",
         "load_raster_maps().")
  if (!is.character(bin_file) || nchar(trimws(bin_file)) == 0)
    stop("'bin_file' must be a non-empty character string.")
  if (grepl(" ", bin_file))
    stop("'bin_file' cannot contain spaces (BayesTraits requirement).")
  if (!is.null(max_age) &&
      (!is.numeric(max_age) || length(max_age) != 1 || is.na(max_age) ||
       max_age < 0))
    stop("'max_age' must be a single non-negative number (Ma).")

  raster   <- identical(.maps_type(maps), "raster")
  use_ages <- if (is.null(ages)) maps$ages
              else if (raster) sort(unique(as.numeric(ages)))
              else as.integer(sort(unique(ages)))
  available <- .paleomaps_available_ages(maps, use_ages)

  if (length(available) == 0)
    stop("No available ages found. Run get_paleomaps() (and fill_missing_maps()) ",
         "or load_raster_maps() first.")

  if (length(available) < length(use_ages)) {
    missing <- setdiff(use_ages, available)
    warning(sprintf(
      "%d age(s) skipped (%s): %s",
      length(missing),
      if (raster) "not registered" else "download failed; run fill_missing_maps() to fix",
      paste(missing, collapse = ", ")
    ))
  }

  available <- sort(available)
  n         <- length(available)

  if (!is.null(max_age) && max_age < available[n])
    stop(sprintf("'max_age' (%g Ma) is younger than the oldest map (%g Ma).",
                 max_age, available[n]))
  if (n == 1 && is.null(max_age))
    warning(sprintf(paste0(
      "Only one map (%g Ma): it covers [0, %g) Ma, so older nodes get no ",
      "restriction map. Set 'max_age' (e.g. the root age) to cover the whole tree."),
      available[1], available[1] + 0.5), call. = FALSE)

  # Midpoints between neighbouring ages, so every point in time is covered
  # by exactly one map (consecutive or irregular spacing):
  #   First map  : [0, midpoint(age[1], age[2]))
  #   Middle maps: [midpoint(age[i-1], age[i]), midpoint(age[i], age[i+1]))
  #   Last map   : [midpoint(age[n-1], age[n]), end)
  #   end = age[n] + 0.5, or ceiling(max_age) + 0.5 when max_age is set
  last_end   <- if (is.null(max_age)) available[n] + 0.5 else ceiling(max_age) + 0.5
  midpoints  <- (available[-n] + available[-1]) / 2.0
  start_ages <- c(0.0,       midpoints)
  end_ages   <- c(midpoints, last_end)

  header    <- sprintf("%s\t%d\t%d", bin_file, as.integer(lon_boxes), as.integer(lat_boxes))
  map_lines <- sprintf("%s%g.csv\t%.1f\t%.1f",
                       prefix, available, start_ages, end_ages)

  writeLines(c(header, map_lines), con = out_file)

  if (verbose) {
    message(sprintf("[write_buildmaps_file] Written: %s", normalizePath(out_file)))
    message(sprintf("  Binary output : %s", bin_file))
    message(sprintf("  Ages included : %d  (%g - %g Ma), covering 0 - %g Ma",
                    n, min(available), max(available), last_end))
    message(sprintf("  Grid          : %d lon x %d lat boxes",
                    as.integer(lon_boxes), as.integer(lat_boxes)))
    message("\nTo build the binary map file, run in BayesTraits:")
    message(sprintf("  BuildMaps %s", out_file))
    message("Then load it in your analysis with:")
    message(sprintf("  LoadMaps %s", bin_file))
  }

  invisible(normalizePath(out_file))
}


#' Re-download palaeomaps that failed in get_paleomaps()
#'
#' The GPlates Web Service occasionally fails for individual ages (time-outs,
#' rate limiting). \code{get_paleomaps()} records those ages in
#' \code{maps$failed} and keeps an empty (\code{NULL}) slot for them. This
#' function retries the download for those ages, or for any ages you supply,
#' and returns the updated \code{geo_paleomaps} object.
#'
#' @param maps    A \code{geo_paleomaps} object from \code{get_paleomaps()}.
#' @param ages    Integer vector or \code{NULL}. Ages (Ma) to (re-)download.
#'                \code{NULL} (default) retries every age that is missing or
#'                listed in \code{maps$failed}. Ages not yet in \code{maps}
#'                are added.
#' @param verbose Logical. Report progress for each age (default \code{TRUE}).
#' @return The updated \code{geo_paleomaps} object. \code{maps$failed}
#'         lists any ages that still could not be downloaded.
#'
#' @examples
#' \dontrun{
#' # Requires an internet connection (GPlates Web Service)
#' maps <- get_paleomaps(max_age = 74)
#' maps <- fill_missing_maps(maps)          # retry failed ages
#' maps <- fill_missing_maps(maps, ages = 75:80)  # add new ages
#' }
#' @export
fill_missing_maps <- function(maps, ages = NULL, verbose = TRUE) {

  if (!inherits(maps, "geo_paleomaps"))
    stop("'maps' must be a geo_paleomaps object from get_paleomaps().")
  if (identical(maps$type, "raster"))
    stop("fill_missing_maps() only applies to maps from get_paleomaps(). ",
         "For local raster maps, add or replace files in the source folder ",
         "and re-run load_raster_maps().")
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
    message("[fill_missing_maps] No missing maps found - nothing to do.")
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
      .clean_coastlines(raw)
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
      message(sprintf("[fill_missing_maps] Still missing: %s - try again later.",
                      paste(still_failed, collapse = ", ")))
    else
      message("[fill_missing_maps] All missing ages resolved.")
  }

  maps
}
