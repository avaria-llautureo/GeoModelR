# ---------------------------------------------------------------------------
# Internal helper: per-edge, per-iteration great-circle distances, with full
# spatial/temporal context attached (coordinates, ape node IDs, node ages,
# fixed tree branch time, and VR-scaled branch time).
# ---------------------------------------------------------------------------
.edge_distances <- function(geo, phy) {

  edges <- as.data.frame(phy$edge)
  names(edges) <- c("parent", "child")
  edges$edge_id        <- seq_len(nrow(edges))
  edges$branch_time_myr <- phy$edge.length   # fixed chronological branch time

  coords <- geo$coords[, c("node_name", "ape_node", "iteration", "longitude", "latitude",
                            "height", "branch_length_scaled")]

  pc <- merge(edges, coords, by.x = "parent", by.y = "ape_node", all.x = TRUE)
  names(pc)[names(pc) %in% c("node_name", "longitude", "latitude", "height")] <-
    c("parent_node_name", "parent_lon", "parent_lat", "parent_height")
  pc$branch_length_scaled <- NULL  # parent's own scaled length is irrelevant here

  full <- merge(pc, coords, by.x = c("child", "iteration"), by.y = c("ape_node", "iteration"),
                all.x = TRUE)
  names(full)[names(full) %in% c("node_name", "longitude", "latitude", "height",
                                  "branch_length_scaled")] <-
    c("child_node_name", "child_lon", "child_lat", "child_height", "branch_length_scaled")

  if (anyNA(full[, c("parent_node_name", "child_node_name")])) {
    warning("[get_dispersal_distance] Some edges could not be matched to nodes in 'geo'. ",
            "Check that 'geo' was parsed with the same tree used here (see add_ape_nodes()).")
  }

  full$distance_km <- geosphere::distHaversine(
    cbind(full$parent_lon, full$parent_lat),
    cbind(full$child_lon,  full$child_lat)
  ) / 1000

  full[order(full$iteration, full$edge_id), ]
}

# ---------------------------------------------------------------------------
# Shared input validation
# ---------------------------------------------------------------------------
.check_geo_phy <- function(geo, phy) {
  if (!inherits(geo, "geo_ancstates")) {
    stop("'geo' must be a geo_ancstates object from read_geo_ancstates().")
  }
  if (!inherits(phy, "phylo")) {
    stop("'phy' must be a phylo object (ape).")
  }
  if (is.null(geo$coords$ape_node) || anyNA(geo$coords$ape_node)) {
    stop("'geo' has no ape_node mapping. Re-run read_geo_ancstates(file, phy = phy) ",
         "or call add_ape_nodes(geo, phy) first.")
  }
}

# ---------------------------------------------------------------------------
# Internal: branch distances, long format, full context attached
# ---------------------------------------------------------------------------
# nodes = "all"  -> every branch (default; matches original behaviour)
# nodes = "tips" -> terminal branches only (parent -> tip), with a tip_label
#                   column giving the actual species name for the child tip.
.branch_distances <- function(ed, phy, nodes = c("all", "tips")) {

  nodes <- match.arg(nodes)
  n_tip <- length(phy$tip.label)

  child_is_tip <- ed$child <= n_tip

  if (nodes == "tips") {
    ed <- ed[child_is_tip, ]
    child_is_tip <- child_is_tip[child_is_tip]
  }

  tip_label <- rep(NA_character_, nrow(ed))
  tip_label[child_is_tip] <- phy$tip.label[ed$child[child_is_tip]]

  out <- data.frame(
    from_node            = ed$parent_node_name,
    to_node               = ed$child_node_name,
    tip_label             = tip_label,
    parent_ape_node        = ed$parent,
    child_ape_node         = ed$child,
    iteration              = ed$iteration,
    from_lon               = ed$parent_lon,
    from_lat               = ed$parent_lat,
    to_lon                 = ed$child_lon,
    to_lat                 = ed$child_lat,
    parent_height          = ed$parent_height,
    child_height           = ed$child_height,
    branch_time_myr        = ed$branch_time_myr,
    branch_length_scaled   = ed$branch_length_scaled,
    distance_km            = ed$distance_km,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Internal: root-to-tip cumulative distances, long format, full context
# attached (current coordinates, node age, and total lineage time elapsed).
# ---------------------------------------------------------------------------
.root_to_tip_distances <- function(ed, geo, phy, nodes = c("tips", "all")) {

  nodes <- match.arg(nodes)

  n_tip   <- length(phy$tip.label)
  n_nodes <- n_tip + phy$Nnode
  root    <- n_tip + 1L

  # Topology is identical across iterations: build the parent-of lookup once.
  topo <- ed[ed$iteration == ed$iteration[1], c("parent", "child")]
  parent_of <- integer(n_nodes)
  parent_of[topo$child] <- topo$parent

  # Node metadata (constant across iterations): name, tip status, age
  node_lookup <- unique(geo$coords[, c("ape_node", "node_name", "is_tip", "height")])
  node_lookup$tip_label <- NA_character_
  node_lookup$tip_label[node_lookup$is_tip] <- phy$tip.label[node_lookup$ape_node[node_lookup$is_tip]]

  root_height <- node_lookup$height[node_lookup$ape_node == root]
  node_lookup$lineage_time_myr <- root_height - node_lookup$height

  # Per-iteration coordinates, for attaching the node's own current position
  coords_it <- geo$coords[, c("ape_node", "iteration", "longitude", "latitude")]

  iterations <- sort(unique(ed$iteration))
  out_list <- vector("list", length(iterations))

  for (i in seq_along(iterations)) {
    it  <- iterations[i]
    sub <- ed[ed$iteration == it, ]

    # dist_to_parent[node] = branch length (km) from parent_of[node] -> node
    dist_to_parent <- numeric(n_nodes)
    dist_to_parent[sub$child] <- sub$distance_km

    cum <- rep(NA_real_, n_nodes)
    cum[root] <- 0

    # Recursive traversal over ape's edge matrix, memoised so each node's
    # cumulative distance is computed once regardless of how many tips
    # share that ancestral path.
    get_cum_dist <- function(node) {
      if (!is.na(cum[node])) return(cum[node])
      val <- get_cum_dist(parent_of[node]) + dist_to_parent[node]
      cum[node] <<- val
      val
    }

    query_nodes <- if (nodes == "tips") seq_len(n_tip) else seq_len(n_nodes)
    for (nd in query_nodes) get_cum_dist(nd)

    out_list[[i]] <- data.frame(
      ape_node          = query_nodes,
      iteration         = it,
      path_distance_km  = cum[query_nodes]
    )
  }

  out <- do.call(rbind, out_list)
  out <- merge(out, node_lookup, by = "ape_node", all.x = TRUE)
  out <- merge(out, coords_it, by = c("ape_node", "iteration"), all.x = TRUE)
  out <- out[, c("node_name", "ape_node", "tip_label", "iteration",
                 "longitude", "latitude", "height", "lineage_time_myr",
                 "path_distance_km", "is_tip")]
  out <- out[order(out$iteration, out$ape_node), ]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Plate frame: reconstruct posterior coordinates to present-day positions
# ---------------------------------------------------------------------------

# Reverse-reconstruct a lon/lat matrix from 'age' Ma to the present.
# Never errors: points that cannot be reconstructed are returned as NA.
# Sent to the GPlates Web Service in chunks so one bad request only affects
# its own chunk. Returns list(xy = matrix, error = character or NULL).
.reconstruct_points <- function(xy, age, model, chunk = 2000L) {
  out    <- matrix(NA_real_, nrow(xy), 2)
  errors <- character(0)
  starts <- seq(1L, nrow(xy), by = chunk)
  for (st in starts) {
    rows <- st:min(st + chunk - 1L, nrow(xy))
    res <- tryCatch(
      withCallingHandlers(
        rgplates::reconstruct(x = xy[rows, , drop = FALSE], age = age,
                              model = model, reverse = TRUE),
        warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) {
        errors <<- c(errors, conditionMessage(e))
        NULL
      })
    if (is.null(res)) next
    res <- as.matrix(res)
    if (nrow(res) != length(rows) || ncol(res) < 2) {
      errors <- c(errors, "unexpected result size from rgplates::reconstruct()")
      next
    }
    out[rows, ] <- suppressWarnings(matrix(as.numeric(res[, 1:2]), ncol = 2))
  }
  list(xy = out, error = if (length(errors) > 0) unique(errors) else NULL)
}

# Replace geo$coords longitude/latitude by present-day positions. Each node
# is reconstructed from its age rounded to the nearest integer Ma (the age
# of the palaeomap it was sampled on). Nodes of age 0 are left unchanged.
# Returns list(geo = geo, failures = data.frame per node).
.coords_to_plate_frame <- function(geo, model, verbose = TRUE) {

  if (!requireNamespace("rgplates", quietly = TRUE))
    stop("Package 'rgplates' is required for frame = \"plate\". ",
         "Install with: install.packages('rgplates')")

  co  <- geo$coords
  age <- vapply(co$height, .snap_age, integer(1))
  ages <- sort(unique(age[age > 0]))
  failed <- rep(FALSE, nrow(co))
  age_errors <- list()

  if (verbose && length(ages) > 0)
    message(sprintf(paste0(
      "[get_dispersal_distance] frame = \"plate\": reconstructing %d ",
      "posterior coordinates at %d age(s) to the present (model '%s')..."),
      sum(age > 0), length(ages), model))

  for (a in ages) {
    idx <- which(age == a)
    rec <- .reconstruct_points(cbind(co$longitude[idx], co$latitude[idx]),
                               age = a, model = model)
    bad <- !stats::complete.cases(rec$xy)
    co$longitude[idx] <- rec$xy[, 1]
    co$latitude[idx]  <- rec$xy[, 2]
    failed[idx] <- bad
    if (!is.null(rec$error)) age_errors[[as.character(a)]] <- rec$error
  }
  geo$coords <- co

  per_node <- stats::aggregate(
    data.frame(n_samples = 1L, n_failed = as.integer(failed)),
    by = list(node_name = co$node_name), FUN = sum)
  meta <- unique(co[, c("node_name", "ape_node", "height")])
  per_node <- merge(meta, per_node, by = "node_name")
  per_node$reconstruction_age <- vapply(per_node$height, .snap_age, integer(1))
  per_node <- per_node[order(per_node$ape_node),
                       c("node_name", "ape_node", "height", "reconstruction_age",
                         "n_samples", "n_failed")]
  rownames(per_node) <- NULL
  attr(per_node, "errors") <- age_errors

  if (any(per_node$n_failed > 0)) {
    nb <- per_node[per_node$n_failed > 0, ]
    warning(sprintf(paste0(
      "[get_dispersal_distance] %d of %d nodes have posterior samples that ",
      "could not be reconstructed to present-day coordinates ",
      "(%d of %d samples in total). These samples are kept with NA ",
      "distances, not dropped.\n",
      "  Per-node counts: attr(result, \"reconstruction\")"),
      nrow(nb), nrow(per_node), sum(nb$n_failed), sum(per_node$n_samples)),
      call. = FALSE)
  } else if (verbose) {
    message("  All samples reconstructed.")
  }

  list(geo = geo, failures = per_node)
}

# ---------------------------------------------------------------------------
# get_dispersal_distance()
# ---------------------------------------------------------------------------
#' Great-circle dispersal distances across the posterior
#'
#' Computes great-circle (Haversine) dispersal distances from a
#' \code{geo_ancstates} object and its associated tree, for every MCMC
#' iteration. Two distance types are available, selected via \code{type}:
#'
#' \itemize{
#'   \item \code{"branch"} - the distance between parent and child node
#'         coordinates for every branch in the tree (one row per branch x
#'         iteration), together with the coordinates, node ages, and both
#'         a fixed-time and a variable-rates (VR) scaled branch time, so
#'         dispersal rates can be computed directly from the output.
#'   \item \code{"path"} - the cumulative distance from the root to each
#'         tip (or, optionally, to every node), obtained by recursively
#'         summing branch distances over \code{phy$edge}. This is the
#'         total posterior dispersal distance per lineage per iteration,
#'         together with each node's age and total elapsed lineage time.
#'   \item \code{"both"} - both of the above, returned together.
#' }
#'
#' The \code{nodes} argument restricts which rows are returned, and its
#' default differs by \code{type} so that existing code is unaffected:
#' \itemize{
#'   \item For \code{type = "branch"}: default \code{"all"} (every branch,
#'         internal and terminal). \code{nodes = "tips"} restricts to
#'         terminal branches only - i.e. the branch connecting each tip to
#'         its immediate ancestor.
#'   \item For \code{type = "path"}: default \code{"tips"} (cumulative
#'         root-to-tip distance per lineage). \code{nodes = "all"} also
#'         returns cumulative distance at every internal node.
#' }
#'
#' @section Reference frame:
#' \describe{
#'   \item{\code{frame = "palaeo"} (default)}{Distances are measured
#'     between the posterior coordinates exactly as sampled by the Geo
#'     model, i.e. palaeocoordinates on each node's palaeomap. A lineage
#'     carried by a moving plate therefore accumulates distance from plate
#'     motion as well as from dispersal. This is the approach of
#'     O'Donovan et al. (2018) and Avaria-Llautureo et al. (2025), and the
#'     behaviour of earlier GeoModelR versions.}
#'   \item{\code{frame = "plate"}}{Each node's posterior coordinates are
#'     first reconstructed to their present-day position on the plate they
#'     sit on (\code{rgplates::reconstruct(..., reverse = TRUE)}), and
#'     distances are measured between those positions. This removes plate
#'     motion, so a lineage that stays on the same spot of a drifting plate
#'     has zero distance. Each node is reconstructed from its age rounded to
#'     the nearest integer Ma (the palaeomap it was sampled on); nodes of
#'     age 0 are not reconstructed. Requires the GPlates Web Service
#'     (internet).
#'
#' Samples that cannot be reconstructed (for example points not assigned
#' to a plate, or a failed request) are \strong{kept} with \code{NA}
#' coordinates and distances; cumulative \code{"path"} distances through
#' them are also \code{NA}. A warning gives the totals, and a per-node
#' table (\code{node_name}, \code{ape_node}, \code{height},
#' \code{reconstruction_age}, \code{n_samples}, \code{n_failed}) is
#' attached as \code{attr(result, "reconstruction")}.}
#' }
#' The frame used is stored in \code{attr(result, "frame")}.
#'
#' @param geo A \code{geo_ancstates} object from \code{\link{read_geo_ancstates}()},
#'   parsed with (or subsequently linked to via \code{\link{add_ape_nodes}()}) the
#'   same tree supplied as \code{phy}.
#' @param phy A \code{phylo} object (\pkg{ape}). Must be the same tree used to
#'   generate \code{geo} (i.e. \code{geo$coords$ape_node} must index this tree).
#' @param type Character. \code{"branch"} (default), \code{"path"}, or
#'   \code{"both"}. See Details.
#' @param nodes Character, \code{"tips"} or \code{"all"}. Optional; if left
#'   as \code{NULL} (the default), the sensible default for \code{type} is
#'   used automatically (see Details).
#' @param frame Character. \code{"palaeo"} (default) or \code{"plate"}. See
#'   the Reference frame section.
#' @param model Character. Plate model used when \code{frame = "plate"}
#'   (default \code{"PALEOMAP"}); use the model your palaeomaps were built
#'   with.
#' @param verbose Logical. Report progress of the plate reconstruction
#'   (default \code{TRUE}).
#'
#' @return
#'   If \code{type = "branch"}, a long-format \code{data.frame} with one row
#'   per branch x iteration:
#'   \itemize{
#'     \item \code{from_node}, \code{to_node} - BayesTraits node names.
#'     \item \code{tip_label} - species name if \code{to_node} is a tip,
#'           else \code{NA}.
#'     \item \code{parent_ape_node}, \code{child_ape_node} - \pkg{ape} node
#'           numbers, for joining back onto \code{geo$coords} or \code{phy}.
#'     \item \code{iteration} - MCMC iteration.
#'     \item \code{from_lon}, \code{from_lat}, \code{to_lon}, \code{to_lat}
#'           - parent/child coordinates for this iteration, in the chosen
#'           frame (present-day positions when \code{frame = "plate"}).
#'     \item \code{parent_height}, \code{child_height} - node ages (fixed,
#'           same units as the input tree).
#'     \item \code{branch_time_myr} - fixed chronological branch length
#'           (from \code{phy$edge.length}; constant across iterations).
#'     \item \code{branch_length_scaled} - variable-rates-scaled branch
#'           length for this edge (from \code{geo$coords}; varies by
#'           iteration). Dividing \code{distance_km} by either time column
#'           gives a fixed-time or VR-scaled dispersal rate respectively.
#'     \item \code{distance_km} - great-circle distance in kilometres.
#'   }
#'
#'   If \code{type = "path"}, a long-format \code{data.frame} with one row
#'   per node x iteration:
#'   \itemize{
#'     \item \code{node_name}, \code{ape_node}, \code{tip_label} - as above.
#'     \item \code{iteration} - MCMC iteration.
#'     \item \code{longitude}, \code{latitude} - this node's own coordinates
#'           for this iteration, in the chosen frame.
#'     \item \code{height} - this node's age (fixed).
#'     \item \code{lineage_time_myr} - total chronological time elapsed
#'           from the root to this node (\code{root height - node height};
#'           fixed). Dividing \code{path_distance_km} by this gives a
#'           whole-lineage dispersal rate.
#'     \item \code{path_distance_km} - cumulative great-circle distance
#'           (km) from the root to this node.
#'     \item \code{is_tip} - logical, \code{TRUE} for tip nodes.
#'   }
#'
#'   If \code{type = "both"}, a list with elements \code{$branch} and
#'   \code{$path}, each as described above.
#'
#' @references
#' O'Donovan, C., Meade, A., & Venditti, C. (2018). Dinosaurs reveal the
#' geographical signature of an evolutionary radiation. Nature Ecology &
#' Evolution, 2, 452-458. \doi{10.1038/s41559-017-0454-6}
#'
#' Avaria-Llautureo, J., Puschel, T. A., Meade, A., Baker, J.,
#' Nicholson, S. L., & Venditti, C. (2025). The radiation and geographic
#' expansion of primates through diverse climates. Proceedings of the
#' National Academy of Sciences, 122(32), e2423833122.
#' \doi{10.1073/pnas.2423833122}
#'
#' @examples
#' phy <- ape::read.nexus(system.file("extdata", "Median.trees",
#'                                     package = "GeoModelR"))
#' geo <- read_geo_ancstates(system.file("extdata", "Primate.AncStates.txt",
#'                                       package = "GeoModelR"), phy = phy)
#'
#' # Terminal branches, with a dispersal rate (km / Myr) per iteration
#' terminal_bd <- get_dispersal_distance(geo, phy, type = "branch", nodes = "tips")
#' terminal_bd$rate_km_per_myr <- terminal_bd$distance_km / terminal_bd$branch_time_myr
#'
#' # Root-to-tip cumulative distance vs. lineage age
#' rt <- get_dispersal_distance(geo, phy, type = "path")
#' rt$rate_km_per_myr <- rt$path_distance_km / rt$lineage_time_myr
#' cor(rt$height, rt$path_distance_km)
#'
#' \dontrun{
#' # Remove plate motion (needs the GPlates Web Service)
#' bd_plate <- get_dispersal_distance(geo, phy, frame = "plate")
#' attr(bd_plate, "reconstruction")   # failed samples per node
#' }
#'
#' @export
get_dispersal_distance <- function(geo, phy,
                                   type    = c("branch", "path", "both"),
                                   nodes   = NULL,
                                   frame   = c("palaeo", "plate"),
                                   model   = "PALEOMAP",
                                   verbose = TRUE) {

  type  <- match.arg(type)
  frame <- match.arg(frame)

  if (!is.null(nodes)) {
    nodes <- match.arg(nodes, choices = c("tips", "all"))
  }

  .check_geo_phy(geo, phy)

  recon <- NULL
  if (frame == "plate") {
    pf    <- .coords_to_plate_frame(geo, model = model, verbose = verbose)
    geo   <- pf$geo
    recon <- pf$failures
  }

  ed <- .edge_distances(geo, phy)

  branch_nodes <- if (is.null(nodes)) "all"  else nodes
  path_nodes   <- if (is.null(nodes)) "tips" else nodes

  tag <- function(x) {
    attr(x, "frame") <- frame
    if (!is.null(recon)) attr(x, "reconstruction") <- recon
    x
  }

  if (type == "branch")
    return(tag(.branch_distances(ed, phy, nodes = branch_nodes)))

  if (type == "path")
    return(tag(.root_to_tip_distances(ed, geo, phy, nodes = path_nodes)))

  # type == "both"
  tag(list(
    branch = tag(.branch_distances(ed, phy, nodes = branch_nodes)),
    path   = tag(.root_to_tip_distances(ed, geo, phy, nodes = path_nodes))
  ))
}
