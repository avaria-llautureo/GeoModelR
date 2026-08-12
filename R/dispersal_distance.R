# ---------------------------------------------------------------------------
# Internal helper: per-edge, per-iteration great-circle distances
# ---------------------------------------------------------------------------
# Not exported. Joins ape's edge matrix (phy$edge, columns parent/child by
# ape node number) onto geo$coords (long: node_name / ape_node / iteration /
# longitude / latitude), then computes the Haversine distance for every
# edge x iteration combination in one vectorised pass.
.edge_distances <- function(geo, phy) {

  edges <- as.data.frame(phy$edge)
  names(edges) <- c("parent", "child")
  edges$edge_id <- seq_len(nrow(edges))

  coords <- geo$coords[, c("node_name", "ape_node", "iteration", "longitude", "latitude")]

  pc <- merge(edges, coords, by.x = "parent", by.y = "ape_node", all.x = TRUE)
  names(pc)[names(pc) %in% c("node_name", "longitude", "latitude")] <-
    c("parent_node_name", "parent_lon", "parent_lat")

  full <- merge(pc, coords, by.x = c("child", "iteration"), by.y = c("ape_node", "iteration"),
                all.x = TRUE)
  names(full)[names(full) %in% c("node_name", "longitude", "latitude")] <-
    c("child_node_name", "child_lon", "child_lat")

  if (anyNA(full[, c("parent_lon", "parent_lat", "child_lon", "child_lat")])) {
    warning("[get_dispersal_distance] Some edges could not be matched to coordinates. ",
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
# Internal: branch distances, long format
# ---------------------------------------------------------------------------
.branch_distances <- function(ed) {
  out <- data.frame(
    from_node   = ed$parent_node_name,
    to_node     = ed$child_node_name,
    iteration   = ed$iteration,
    distance_km = ed$distance_km,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Internal: root-to-tip cumulative distances, long format
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

  # Node-name / tip-status lookup (constant across iterations)
  node_lookup <- unique(geo$coords[, c("ape_node", "node_name", "is_tip")])

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
  out <- out[, c("node_name", "ape_node", "iteration", "path_distance_km", "is_tip")]
  out <- out[order(out$iteration, out$ape_node), ]
  rownames(out) <- NULL
  out
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
#'   \item \code{"branch"} — the distance between parent and child node
#'         coordinates for every branch in the tree (one row per branch x
#'         iteration).
#'   \item \code{"path"} — the cumulative distance from the root to each
#'         tip (or, optionally, to every node), obtained by recursively
#'         summing branch distances over \code{phy$edge}. This is the
#'         total posterior dispersal distance per lineage per iteration.
#'   \item \code{"both"} — both of the above, returned together.
#' }
#'
#' @param geo A \code{geo_ancstates} object from \code{\link{read_geo_ancstates}()},
#'   parsed with (or subsequently linked to via \code{\link{add_ape_nodes}()}) the
#'   same tree supplied as \code{phy}.
#' @param phy A \code{phylo} object (\pkg{ape}). Must be the same tree used to
#'   generate \code{geo} (i.e. \code{geo$coords$ape_node} must index this tree).
#' @param type Character. \code{"branch"} (default), \code{"path"}, or
#'   \code{"both"}. See Details.
#' @param nodes Character. Only used when \code{type} is \code{"path"} or
#'   \code{"both"}. \code{"tips"} (default) returns cumulative distance for
#'   tips only; \code{"all"} also returns cumulative distance at every
#'   internal node (root is always 0).
#'
#' @return
#'   If \code{type = "branch"}, a long-format \code{data.frame} with one row
#'   per branch x iteration: \code{from_node}, \code{to_node},
#'   \code{iteration}, \code{distance_km}.
#'
#'   If \code{type = "path"}, a long-format \code{data.frame} with one row
#'   per node x iteration: \code{node_name}, \code{ape_node},
#'   \code{iteration}, \code{path_distance_km}, \code{is_tip}.
#'
#'   If \code{type = "both"}, a list with elements \code{$branch} and
#'   \code{$path}, each as described above.
#'
#' @examples
#' \dontrun{
#' phy <- ape::read.nexus("Median.trees")
#' geo <- read_geo_ancstates("Primate.AncStates.txt", phy = phy)
#'
#' # Per-branch great-circle distances
#' bd <- get_dispersal_distance(geo, phy, type = "branch")
#'
#' # Root-to-tip cumulative dispersal distance
#' rt <- get_dispersal_distance(geo, phy, type = "path")
#' aggregate(path_distance_km ~ node_name, data = rt, FUN = mean)
#'
#' # Both at once
#' both <- get_dispersal_distance(geo, phy, type = "both")
#' head(both$branch)
#' head(both$path)
#' }
#'
#' @export
get_dispersal_distance <- function(geo, phy,
                                    type  = c("branch", "path", "both"),
                                    nodes = c("tips", "all")) {

  type  <- match.arg(type)
  nodes <- match.arg(nodes)

  .check_geo_phy(geo, phy)

  ed <- .edge_distances(geo, phy)

  if (type == "branch") {
    return(.branch_distances(ed))
  }

  if (type == "path") {
    return(.root_to_tip_distances(ed, geo, phy, nodes = nodes))
  }

  # type == "both"
  list(
    branch = .branch_distances(ed),
    path   = .root_to_tip_distances(ed, geo, phy, nodes = nodes)
  )
}
