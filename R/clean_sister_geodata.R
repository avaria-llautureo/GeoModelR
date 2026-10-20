# =============================================================================
# GeoModelR - Geographic input data cleaning
#
# clean_sister_geodata()
#   Check and clean BayesTraits v5+ Geo model input coordinates so that
#   sister taxa do not share (near-)identical locations, which can compress
#   branches to zero length under the variable-rates model.
#
# Authors: Jorge Avaria-Llautureo & Andrew Meade
# =============================================================================

#' Clean overlapping sister-taxon coordinates for the BayesTraits Geo model
#'
#' Avaria-Llautureo et al. (2025) found that the variable-rates Geo model can
#' compress the branches between sister taxa to zero length when their input
#' coordinates are less than about 50 m apart. This function finds sister
#' taxa whose coordinates overlap (closer than \code{distance_threshold}) and
#' removes the overlapping distribution points before the data are given to
#' BayesTraits. Overlapping single centroids are never changed automatically:
#' they are reported, and you decide which taxon to remove.
#'
#' @section Input data rules (checked first):
#' \enumerate{
#'   \item Taxa in \code{singledata} or \code{distdata} that are not in the
#'         tree are dropped and reported in \code{$not_in_tree}.
#'   \item Every tree tip must have a centroid (one row) in
#'         \code{singledata}; otherwise the function stops.
#'   \item Every taxon in \code{distdata} must also have a centroid in
#'         \code{singledata}; otherwise the function stops.
#' }
#'
#' @section Effective coordinates:
#' For a taxon with two or more distribution points, BayesTraits uses only
#' those points and ignores its centroid (the centroid row must still be in
#' the file). A taxon with fewer than two distribution points uses its
#' centroid. Overlap is always checked between these \emph{effective}
#' coordinate sets.
#'
#' @section Which pairs are compared:
#' \itemize{
#'   \item Sister tips: every pair of tips that share the same parent node
#'         (cherries, and all tip pairs within a polytomy).
#'   \item Triples (\code{include_triples = TRUE}): a tip that is sister to a
#'         two-tip clade is compared with each tip of that clade.
#' }
#'
#' @section How overlaps are resolved:
#' \describe{
#'   \item{Both taxa use distribution points}{The overlapping points are
#'     deleted from the taxon with \strong{more} distribution points; the
#'     other taxon is left unchanged. If both have the same number of points,
#'     the taxon is chosen at random (set \code{seed} for reproducible
#'     results). Every deleted point is recorded in \code{$deletion_log}.}
#'   \item{One taxon uses distribution points, the other its centroid}{The
#'     distribution points that overlap the centroid are deleted.}
#'   \item{Both taxa use their centroid}{Nothing is deleted or moved. The
#'     function stops and lists the conflicting pairs. Resolve them by
#'     correcting the coordinates, or by naming the taxa to remove in
#'     \code{delete} and running the function again. With
#'     \code{ask = TRUE} (interactive sessions only) a menu asks which taxon
#'     of each pair to remove and fills in \code{delete} for you.}
#' }
#' When deletions leave a taxon with fewer than two distribution points it
#' falls back to its centroid (reported in \code{$fallback_to_centroid}).
#' Because a fallback can create a new overlap, all pairs are compared again
#' until nothing changes or \code{max_iter} passes have run.
#'
#' Taxa named in \code{delete} are removed from both data sets and from the
#' tree before the comparison (their former sisters are then compared with
#' their new sisters). Use the returned \code{$clean_tree} in BayesTraits
#' together with the cleaned data.
#'
#' The conflict error has class \code{"geo_centroid_conflict"} and carries
#' a \code{conflicts} data frame, so it can also be handled in code with
#' \code{tryCatch()}.
#'
#' @param tree               A \code{phylo} or \code{multiPhylo} object. Sister
#'                           relationships are taken from the first tree.
#' @param singledata         Data frame with columns \code{V1} (taxon),
#'                           \code{V2} (longitude), \code{V3} (latitude): the
#'                           single-coordinate BayesTraits input file, e.g.
#'                           read with \code{read.table(file, sep = "\\t")}.
#' @param distdata           Optional data frame with columns \code{V1}
#'                           (taxon), \code{V2} (e.g. \code{"Linked"}),
#'                           \code{V3} (comma-separated longitudes) and
#'                           \code{V4} (comma-separated latitudes): the
#'                           distribution-data BayesTraits input file.
#' @param distance_threshold Numeric. Great-circle (Haversine) distance in
#'                           kilometres below which two points overlap
#'                           (default \code{0.05}, i.e. 50 m).
#' @param include_distdata   Logical. Use \code{distdata} in the check
#'                           (default \code{TRUE} when \code{distdata} is
#'                           supplied).
#' @param include_triples    Logical. Also compare triples (default
#'                           \code{TRUE}).
#' @param delete             Character vector or \code{NULL}. Taxa to remove
#'                           (from the data and the tree) to resolve centroid
#'                           conflicts.
#' @param ask                Logical. In an interactive session, show a menu
#'                           for each centroid conflict instead of stopping
#'                           straight away (default \code{FALSE}).
#' @param seed               Integer or \code{NULL}. Seed for the random
#'                           choice between two taxa with the same number of
#'                           distribution points. The global random number
#'                           stream is restored afterwards.
#' @param max_iter           Integer. Maximum number of comparison passes
#'                           (default 10).
#' @param verbose            Logical. Print a summary (default \code{TRUE}).
#'
#' @return An object of class \code{geo_cleaned} (a list):
#' \describe{
#'   \item{\code{clean_singledata}}{The input \code{singledata} rows for the
#'         retained taxa, unchanged.}
#'   \item{\code{clean_distdata}}{Distribution data (\code{V1}--\code{V4})
#'         for taxa that still have two or more points. Retained points keep
#'         their original text formatting.}
#'   \item{\code{clean_tree}}{\code{tree} with the \code{delete} taxa
#'         removed (identical to \code{tree} if \code{delete} is
#'         \code{NULL}).}
#'   \item{\code{deleted_taxa}}{Taxa removed through \code{delete}.}
#'   \item{\code{not_in_tree}}{List with elements \code{singledata} and
#'         \code{distdata}: taxa dropped because they are not in the tree.}
#'   \item{\code{fallback_to_centroid}}{Taxa left with fewer than two
#'         distribution points, which BayesTraits will place at their
#'         centroid. Worth reviewing.}
#'   \item{\code{deletion_log}}{Data frame with one row per deleted
#'         distribution point: \code{iteration}, \code{sister_pair},
#'         \code{deleted_from}, \code{lon}, \code{lat}, \code{distance_km}
#'         (to the nearest point of the sister) and \code{reason}.}
#'   \item{\code{taxa_list}}{Per-taxon list with the \code{centroid} and the
#'         remaining distribution points (\code{dist}).}
#'   \item{\code{n_pairs}, \code{n_iterations}}{Number of taxon pairs
#'         compared and comparison passes run.}
#'   \item{\code{settings}}{The arguments used.}
#' }
#'
#' @references Avaria-Llautureo, J., et al. (2025). Proceedings of the
#'   National Academy of Sciences of the USA. \doi{10.1073/pnas.2423833122}
#'
#' @author Jorge Avaria-Llautureo and Andrew Meade
#'
#' @examples
#' tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
#'
#' # A and B have centroids about 20 m apart
#' single <- data.frame(V1 = c("A", "B", "C", "D", "E"),
#'                      V2 = c(-70, -70.0002, -60, 10, 10.2),
#'                      V3 = c(-30, -30, -20, 45, 45.2))
#'
#' # D and E share one distribution point; D has more points
#' dist <- data.frame(V1 = c("D", "E"), V2 = "Linked",
#'                    V3 = c("10,10.5,11", "10.00001,12"),
#'                    V4 = c("45,45.5,46", "45,47"))
#'
#' # Stops: the A and B centroids overlap
#' res <- try(clean_sister_geodata(tree, single, dist), silent = TRUE)
#' cat(res)
#'
#' # Decide which taxon to remove and run again
#' res <- clean_sister_geodata(tree, single, dist, delete = "B")
#' res$deletion_log
#' res$clean_distdata
#'
#' \dontrun{
#' # Write the cleaned files and tree for BayesTraits
#' write.table(res$clean_singledata, "Single_Coord_Data_clean.txt", sep = "\t",
#'             row.names = FALSE, col.names = FALSE, quote = FALSE)
#' write.table(res$clean_distdata, "Distribution_Data_clean.txt", sep = "\t",
#'             row.names = FALSE, col.names = FALSE, quote = FALSE)
#' ape::write.nexus(res$clean_tree, file = "tree_clean.nex")
#' }
#' @export
clean_sister_geodata <- function(tree,
                                 singledata,
                                 distdata           = NULL,
                                 distance_threshold = 0.05,
                                 include_distdata   = !is.null(distdata),
                                 include_triples    = TRUE,
                                 delete             = NULL,
                                 ask                = FALSE,
                                 seed               = NULL,
                                 max_iter           = 10,
                                 verbose            = TRUE) {

  # Inputs as supplied, for the re-run after an interactive choice
  input <- list(tree = tree, singledata = singledata, distdata = distdata)

  # ---------------------------------------------------------------------------
  # 0. Argument checks
  # ---------------------------------------------------------------------------
  if (!inherits(tree, "phylo") && !inherits(tree, "multiPhylo"))
    stop("'tree' must be a phylo or multiPhylo object.")
  if (!is.data.frame(singledata) || !all(c("V1", "V2", "V3") %in% names(singledata)))
    stop("'singledata' must be a data frame with columns V1 (taxon), V2 (lon), V3 (lat).")
  if (include_distdata) {
    if (is.null(distdata))
      stop("include_distdata = TRUE but 'distdata' is NULL.")
    if (!is.data.frame(distdata) || !all(c("V1", "V3", "V4") %in% names(distdata)))
      stop("'distdata' must be a data frame with columns V1 (taxon), V3 (lons), V4 (lats).")
  } else {
    distdata <- NULL
  }
  if (!is.numeric(distance_threshold) || length(distance_threshold) != 1 ||
      is.na(distance_threshold) || distance_threshold < 0)
    stop("'distance_threshold' must be a single non-negative number (km).")
  if (!is.null(delete) && !is.character(delete))
    stop("'delete' must be a character vector of taxon names, or NULL.")

  if (!is.null(seed)) {
    old_seed <- .get_rng_state()
    on.exit(.set_rng_state(old_seed), add = TRUE)
    set.seed(seed)
  }

  settings <- list(distance_threshold = distance_threshold,
                   include_distdata = include_distdata,
                   include_triples = include_triples,
                   delete = delete, seed = seed, max_iter = max_iter)

  # ---------------------------------------------------------------------------
  # 1. Input data rules
  # ---------------------------------------------------------------------------
  phy       <- if (inherits(tree, "multiPhylo")) tree[[1]] else tree
  tree_tips <- phy$tip.label

  singledata$V1 <- as.character(singledata$V1)
  if (!is.null(distdata)) distdata$V1 <- as.character(distdata$V1)

  dup_s <- unique(singledata$V1[duplicated(singledata$V1)])
  if (length(dup_s) > 0)
    stop("Taxa with more than one row in 'singledata': ",
         paste(dup_s, collapse = ", "))
  if (!is.null(distdata)) {
    dup_d <- unique(distdata$V1[duplicated(distdata$V1)])
    if (length(dup_d) > 0)
      stop("Taxa with more than one row in 'distdata': ",
           paste(dup_d, collapse = ", "))
  }

  # Rule 1 - drop taxa that are not in the tree
  not_in_tree <- list(
    singledata = setdiff(singledata$V1, tree_tips),
    distdata   = if (is.null(distdata)) character(0) else
                   setdiff(distdata$V1, tree_tips)
  )
  if (verbose) {
    for (src in names(not_in_tree)) {
      if (length(not_in_tree[[src]]) > 0)
        message(sprintf("%d taxon/taxa in %s are not in the tree and are ignored: %s",
                        length(not_in_tree[[src]]), src,
                        paste(not_in_tree[[src]], collapse = ", ")))
    }
  }
  singledata <- singledata[singledata$V1 %in% tree_tips, , drop = FALSE]
  if (!is.null(distdata))
    distdata <- distdata[distdata$V1 %in% tree_tips, , drop = FALSE]

  # Rule 2 - every tip needs a centroid
  missing_centroid <- setdiff(tree_tips, singledata$V1)
  if (length(missing_centroid) > 0)
    stop(sprintf(paste0(
      "%d tree tip(s) have no centroid in singledata: %s\n",
      "  Every tip in the tree must have a row in singledata."),
      length(missing_centroid), paste(missing_centroid, collapse = ", ")))

  # Rule 3 - every distdata taxon needs a centroid (implied by rules 1-2,
  # kept explicit for clarity)
  if (!is.null(distdata)) {
    dist_no_centroid <- setdiff(distdata$V1, singledata$V1)
    if (length(dist_no_centroid) > 0)
      stop(sprintf(paste0(
        "%d taxon/taxa have distribution data but no centroid in singledata: %s"),
        length(dist_no_centroid), paste(dist_no_centroid, collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # 2. User deletions (resolution of centroid conflicts)
  # ---------------------------------------------------------------------------
  delete <- unique(delete)
  if (length(delete) > 0) {
    unknown <- setdiff(delete, singledata$V1)
    if (length(unknown) > 0)
      stop("Taxa in 'delete' not found in the data/tree: ",
           paste(unknown, collapse = ", "))
    if (length(setdiff(tree_tips, delete)) < 2)
      stop("'delete' would leave fewer than two tips in the tree.")
    singledata <- singledata[!singledata$V1 %in% delete, , drop = FALSE]
    if (!is.null(distdata))
      distdata <- distdata[!distdata$V1 %in% delete, , drop = FALSE]
    phy  <- ape::drop.tip(phy, delete)
    tree <- if (inherits(tree, "multiPhylo")) {
      structure(lapply(seq_along(tree), function(i) ape::drop.tip(tree[[i]], delete)),
                class = "multiPhylo")
    } else phy
  }

  # ---------------------------------------------------------------------------
  # 3. Per-taxon coordinates and the pairs to compare
  # ---------------------------------------------------------------------------
  taxa_list <- .build_taxa_list(singledata, distdata)
  pairs     <- .get_sister_pairs(phy, include_triples)

  # ---------------------------------------------------------------------------
  # 4. Compare pairs until nothing changes
  # ---------------------------------------------------------------------------
  log_rows  <- list()
  conflicts <- list()
  fallback  <- character(0)
  iter      <- 0L
  changed   <- TRUE

  while (changed && iter < max_iter) {
    iter      <- iter + 1L
    changed   <- FALSE
    conflicts <- list()          # only conflicts of the final pass matter

    for (pr in pairs) {
      res <- .resolve_pair(taxa_list, pr[1], pr[2], distance_threshold, iter)
      taxa_list <- res$taxa_list
      if (length(res$log) > 0) {
        log_rows <- c(log_rows, res$log)
        changed  <- TRUE
      }
      fallback <- c(fallback, res$fallback)
      if (!is.null(res$conflict))
        conflicts[[length(conflicts) + 1L]] <- res$conflict
    }
  }
  if (changed)
    warning(sprintf(paste0(
      "Reached max_iter = %d with deletions still happening in the last pass. ",
      "Increase max_iter."), max_iter))

  # ---------------------------------------------------------------------------
  # 5. Centroid conflicts: ask (optional) or stop
  # ---------------------------------------------------------------------------
  if (length(conflicts) > 0) {
    conflicts_df <- do.call(rbind, conflicts)
    rownames(conflicts_df) <- NULL

    if (isTRUE(ask) && .is_interactive()) {
      chosen <- .ask_conflicts(conflicts_df)
      if (length(chosen) > 0) {
        new_delete <- unique(c(delete, chosen))
        message("Re-running with delete = ", .deparse_chr(new_delete),
                "\n  (use this argument to reproduce the result)")
        return(clean_sister_geodata(
          tree = input$tree, singledata = input$singledata,
          distdata = input$distdata, distance_threshold = distance_threshold,
          include_distdata = include_distdata, include_triples = include_triples,
          delete = new_delete, ask = ask, seed = seed, max_iter = max_iter,
          verbose = verbose))
      }
    }

    stop(.conflict_condition(conflicts_df, delete, distance_threshold))
  }

  # ---------------------------------------------------------------------------
  # 6. Assemble outputs
  # ---------------------------------------------------------------------------
  deletion_log <- if (length(log_rows) > 0) do.call(rbind, log_rows) else
    data.frame(iteration = integer(), sister_pair = character(),
               deleted_from = character(), lon = numeric(), lat = numeric(),
               distance_km = numeric(), reason = character(),
               stringsAsFactors = FALSE)
  rownames(deletion_log) <- NULL

  fallback <- unique(fallback)
  fallback <- fallback[vapply(fallback, function(nm) nrow(taxa_list[[nm]]$dist) < 2,
                              logical(1))]

  clean_singledata <- singledata
  rownames(clean_singledata) <- NULL
  clean_distdata <- .rebuild_distdata(taxa_list)

  result <- structure(
    list(
      clean_singledata     = clean_singledata,
      clean_distdata       = clean_distdata,
      clean_tree           = tree,
      deleted_taxa         = delete,
      not_in_tree          = not_in_tree,
      fallback_to_centroid = fallback,
      deletion_log         = deletion_log,
      taxa_list            = lapply(taxa_list, function(t)
                               list(centroid = t$centroid,
                                    dist = t$dist[, c("lon", "lat"), drop = FALSE])),
      n_pairs              = length(pairs),
      n_iterations         = iter,
      settings             = settings
    ),
    class = "geo_cleaned"
  )

  if (verbose) print(result)
  invisible(result)
}


# =============================================================================
# Internal helpers
# =============================================================================

# All taxon pairs to compare: tips sharing a parent (cherries and polytomies)
# and, optionally, a tip against each member of its sister cherry (triples).
# Returns a list of length-2 character vectors.
.get_sister_pairs <- function(phy, include_triples = TRUE) {
  n_tips   <- ape::Ntip(phy)
  children <- split(phy$edge[, 2], phy$edge[, 1])
  lab      <- phy$tip.label
  pairs    <- list()

  for (nd in names(children)) {
    ch   <- children[[nd]]
    tips <- ch[ch <= n_tips]

    if (length(tips) >= 2) {
      cmb <- utils::combn(tips, 2)
      for (j in seq_len(ncol(cmb)))
        pairs[[length(pairs) + 1L]] <- lab[cmb[, j]]
    }

    if (include_triples && length(ch) == 2 && length(tips) == 1) {
      gch <- children[[as.character(ch[ch > n_tips])]]
      if (length(gch) == 2 && all(gch <= n_tips)) {
        pairs[[length(pairs) + 1L]] <- lab[c(tips, gch[1])]
        pairs[[length(pairs) + 1L]] <- lab[c(tips, gch[2])]
      }
    }
  }
  pairs
}


# Split a comma-separated coordinate string, keeping the original text.
.split_coords <- function(x) {
  txt <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  txt <- txt[nzchar(txt)]
  list(txt = txt, num = suppressWarnings(as.numeric(txt)))
}


# Per-taxon list: $centroid (1-row data.frame lon, lat) and $dist
# (data.frame lon, lat, lon_txt, lat_txt; possibly 0 rows), plus the
# original V2 value of the distribution row ($dist_v2).
.build_taxa_list <- function(singledata, distdata) {

  lon <- suppressWarnings(as.numeric(singledata$V2))
  lat <- suppressWarnings(as.numeric(singledata$V3))
  bad <- singledata$V1[is.na(lon) | is.na(lat) | abs(lon) > 180 | abs(lat) > 90]
  if (length(bad) > 0)
    stop("Missing or invalid centroid coordinates in singledata for: ",
         paste(bad, collapse = ", "))

  empty_dist <- data.frame(lon = numeric(0), lat = numeric(0),
                           lon_txt = character(0), lat_txt = character(0),
                           stringsAsFactors = FALSE)

  taxa_list <- lapply(seq_len(nrow(singledata)), function(i)
    list(centroid = data.frame(lon = lon[i], lat = lat[i]),
         dist = empty_dist, dist_v2 = "Linked"))
  names(taxa_list) <- singledata$V1

  if (!is.null(distdata) && nrow(distdata) > 0) {
    bad <- character(0)
    for (i in seq_len(nrow(distdata))) {
      nm  <- distdata$V1[i]
      lo  <- .split_coords(distdata$V3[i])
      la  <- .split_coords(distdata$V4[i])
      if (length(lo$num) != length(la$num) || anyNA(lo$num) || anyNA(la$num) ||
          any(abs(lo$num) > 180) || any(abs(la$num) > 90)) {
        bad <- c(bad, nm)
        next
      }
      taxa_list[[nm]]$dist <- data.frame(lon = lo$num, lat = la$num,
                                         lon_txt = lo$txt, lat_txt = la$txt,
                                         stringsAsFactors = FALSE)
      if ("V2" %in% names(distdata) && !is.na(distdata$V2[i]))
        taxa_list[[nm]]$dist_v2 <- as.character(distdata$V2[i])
    }
    if (length(bad) > 0)
      stop(paste0("Invalid distribution coordinates in distdata for: ",
                  paste(bad, collapse = ", "),
                  "\n  Check that V3 and V4 hold the same number of ",
                  "comma-separated numeric values within valid ranges."))
  }
  taxa_list
}


# Effective coordinate set: distribution points if >= 2, else the centroid.
.effective_set <- function(taxon) {
  if (nrow(taxon$dist) >= 2)
    return(list(coords = taxon$dist[, c("lon", "lat")], type = "dist"))
  list(coords = taxon$centroid, type = "centroid")
}


# Great-circle distance matrix (km) between two lon/lat sets.
.dist_km <- function(a, b) {
  geosphere::distm(as.matrix(a[, c("lon", "lat")]),
                   as.matrix(b[, c("lon", "lat")]),
                   fun = geosphere::distHaversine) / 1000
}


# Compare one pair and resolve any overlap.
# Returns list(taxa_list, log = list of 1-row data.frames, fallback = chr,
#              conflict = NULL or 1-row data.frame)
.resolve_pair <- function(taxa_list, s1, s2, threshold_km, iter) {

  out <- list(taxa_list = taxa_list, log = list(), fallback = character(0),
              conflict = NULL)
  t1 <- taxa_list[[s1]]
  t2 <- taxa_list[[s2]]
  if (is.null(t1) || is.null(t2)) return(out)

  e1 <- .effective_set(t1)
  e2 <- .effective_set(t2)
  dmat  <- .dist_km(e1$coords, e2$coords)
  close <- dmat < threshold_km
  if (!any(close)) return(out)

  pair_lab <- paste(s1, "vs", s2)

  # --- Both centroids: never changed automatically --------------------------
  if (e1$type == "centroid" && e2$type == "centroid") {
    out$conflict <- data.frame(taxon_1 = s1, taxon_2 = s2,
                               distance_km = round(dmat[1, 1], 5),
                               stringsAsFactors = FALSE)
    return(out)
  }

  # --- Choose the taxon that loses points ------------------------------------
  if (e1$type == "dist" && e2$type == "dist") {
    n1 <- nrow(e1$coords); n2 <- nrow(e2$coords)
    if (n1 == n2) {
      loser  <- if (stats::runif(1) < 0.5) 1L else 2L
      reason <- "overlapping distribution point; equal number of points, taxon chosen at random"
    } else {
      loser  <- if (n1 > n2) 1L else 2L
      reason <- "overlapping distribution point; deleted from the taxon with more points"
    }
  } else {
    loser  <- if (e1$type == "dist") 1L else 2L
    reason <- "distribution point overlapping the sister's centroid"
  }

  if (loser == 1L) {
    nm <- s1; rows <- which(rowSums(close) > 0)
    dmin <- apply(dmat[rows, , drop = FALSE], 1, min)
  } else {
    nm <- s2; rows <- which(colSums(close) > 0)
    dmin <- apply(dmat[, rows, drop = FALSE], 2, min)
  }

  pts <- taxa_list[[nm]]$dist[rows, , drop = FALSE]
  out$log <- list(data.frame(
    iteration    = iter,
    sister_pair  = pair_lab,
    deleted_from = nm,
    lon          = pts$lon,
    lat          = pts$lat,
    distance_km  = round(unname(dmin), 5),
    reason       = reason,
    stringsAsFactors = FALSE))

  remaining <- taxa_list[[nm]]$dist[-rows, , drop = FALSE]
  rownames(remaining) <- NULL
  taxa_list[[nm]]$dist <- remaining
  if (nrow(remaining) < 2) out$fallback <- nm

  out$taxa_list <- taxa_list
  out
}


# Rebuild the distribution data frame (V1-V4) from the remaining points,
# re-using the original text of each coordinate.
.rebuild_distdata <- function(taxa_list) {
  keep <- names(taxa_list)[vapply(taxa_list, function(t) nrow(t$dist) >= 2,
                                  logical(1))]
  if (length(keep) == 0)
    return(data.frame(V1 = character(), V2 = character(), V3 = character(),
                      V4 = character(), stringsAsFactors = FALSE))
  data.frame(
    V1 = keep,
    V2 = vapply(taxa_list[keep], `[[`, character(1), "dist_v2"),
    V3 = vapply(taxa_list[keep], function(t) paste(t$dist$lon_txt, collapse = ","),
                character(1)),
    V4 = vapply(taxa_list[keep], function(t) paste(t$dist$lat_txt, collapse = ","),
                character(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}


# Build the classed error for centroid conflicts.
.conflict_condition <- function(conflicts_df, delete, threshold_km) {
  lines <- sprintf("  %-4s %s  vs  %s   (%.1f m)",
                   paste0(seq_len(nrow(conflicts_df)), "."),
                   conflicts_df$taxon_1, conflicts_df$taxon_2,
                   conflicts_df$distance_km * 1000)
  ex_del  <- conflicts_df$taxon_1[1]
  ex_keep <- conflicts_df$taxon_2[1]
  example <- .deparse_chr(unique(c(delete, ex_del)))
  msg <- paste0(
    nrow(conflicts_df), " sister pair(s) have centroids closer than ",
    threshold_km * 1000, " m. Centroids are not changed automatically:\n",
    paste(lines, collapse = "\n"), "\n",
    "Resolve each pair by correcting the coordinates in singledata, or by ",
    "removing one taxon per pair and running again. For example, to keep ",
    ex_keep, " and remove ", ex_del, ":\n",
    "  clean_sister_geodata(..., delete = ", example, ")\n",
    "(add ask = TRUE in an interactive session to choose from a menu).")
  structure(class = c("geo_centroid_conflict", "error", "condition"),
            list(message = msg, call = NULL, conflicts = conflicts_df))
}


# Interactive menu: one question per conflicting pair. Returns the chosen
# taxa to delete (possibly empty).
.ask_conflicts <- function(conflicts_df) {
  chosen <- character(0)
  for (i in seq_len(nrow(conflicts_df))) {
    a <- conflicts_df$taxon_1[i]; b <- conflicts_df$taxon_2[i]
    if (a %in% chosen || b %in% chosen) next
    pick <- menu(
      c(paste("Delete", a), paste("Delete", b), "Keep both (fix the coordinates myself)"),
      title = sprintf("Conflict %d of %d: %s vs %s are %.1f m apart. Which taxon should be removed?",
                      i, nrow(conflicts_df), a, b, conflicts_df$distance_km[i] * 1000))
    if (pick == 1L) chosen <- c(chosen, a)
    if (pick == 2L) chosen <- c(chosen, b)
  }
  chosen
}

.is_interactive <- function() interactive()

.deparse_chr <- function(x) paste0("c(", paste0('"', x, '"', collapse = ", "), ")")

.get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
}
.set_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
}


# =============================================================================
# S3 print method
# =============================================================================

#' Print a geo_cleaned object
#'
#' @param x   A \code{geo_cleaned} object from \code{clean_sister_geodata()}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.geo_cleaned <- function(x, ...) {
  s <- x$settings
  cat("Sister-taxon coordinate cleaning (clean_sister_geodata)\n")
  cat("  Threshold          :", s$distance_threshold * 1000, "m\n")
  cat("  Distribution data  :", s$include_distdata, "\n")
  cat("  Triples checked    :", s$include_triples, "\n")
  cat("  Pairs compared     :", x$n_pairs, "\n")
  cat("  Passes run         :", x$n_iterations, "\n")
  if (!is.null(s$seed)) cat("  Random seed        :", s$seed, "\n")

  n_out <- length(x$not_in_tree$singledata) + length(x$not_in_tree$distdata)
  if (n_out > 0) {
    cat("\n  Dropped (not in tree):\n")
    if (length(x$not_in_tree$singledata) > 0)
      cat("    singledata:", paste(x$not_in_tree$singledata, collapse = ", "), "\n")
    if (length(x$not_in_tree$distdata) > 0)
      cat("    distdata  :", paste(x$not_in_tree$distdata, collapse = ", "), "\n")
  }

  if (length(x$deleted_taxa) > 0) {
    cat("\n  Removed via 'delete' (data and tree):\n")
    cat("   ", paste(x$deleted_taxa, collapse = ", "), "\n")
  }

  cat("\n  Distribution points deleted:", nrow(x$deletion_log), "\n")
  if (nrow(x$deletion_log) > 0) {
    tab <- table(x$deletion_log$deleted_from)
    cat("    by taxon:", paste0(names(tab), " (", tab, ")", collapse = ", "), "\n")
  }
  cat("  Taxa in clean_singledata  :", nrow(x$clean_singledata), "\n")
  cat("  Taxa in clean_distdata    :", nrow(x$clean_distdata), "\n")

  if (length(x$fallback_to_centroid) > 0) {
    cat("\n  NOTE - now placed at their centroid (< 2 distribution points left):\n")
    cat("   ", paste(x$fallback_to_centroid, collapse = ", "), "\n")
  }

  cat("\nAccess: $clean_singledata | $clean_distdata | $clean_tree |",
      "$deletion_log | $fallback_to_centroid\n")
  invisible(x)
}
