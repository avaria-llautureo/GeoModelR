# =============================================================================
# read_geo_ancstates()
# Parser for BayesTraits v5+ Geo model output files (*.AncStates.txt)
#
# Returns a geo_ancstates object with:
#   $nodes  — data.frame of node metadata (one row per node)
#   $mcmc   — data.frame of MCMC posterior samples (wide format)
#   $coords — tidy long-format table (one row per node × iteration)
#
# If a phylo object is supplied via the phy argument, an ape_node column is
# added to $nodes and $coords mapping each BayesTraits node name to its
# corresponding R node number in the ape tree. This can also be done after
# parsing with add_ape_nodes(geo, phy). The ape tree object must be the tree used for analyses.
#
# Author:  Jorge Avaria-Llautureo
# Depends: base R only (ape required only when phy is supplied)
# =============================================================================

read_geo_ancstates <- function(file, phy = NULL, verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Read all lines
  # ---------------------------------------------------------------------------
  raw <- readLines(file, warn = FALSE)
  raw <- raw[nzchar(trimws(raw))]

  # ---------------------------------------------------------------------------
  # 2. Locate block boundaries
  #    Block 1 header: "Node Name"
  #    Block 2 header: "Itter"
  # ---------------------------------------------------------------------------
  block1_header_idx <- which(startsWith(raw, "Node Name"))[1]
  block2_header_idx <- which(startsWith(raw, "Itter"))[1]

  if (is.na(block1_header_idx) || is.na(block2_header_idx)) {
    stop(
      "Could not locate expected block headers.\n",
      "  Block 1 must start with 'Node Name'\n",
      "  Block 2 must start with 'Itter'\n",
      "  Check that the file is a valid BayesTraits AncStates.txt output."
    )
  }

  node_lines <- raw[block1_header_idx:(block2_header_idx - 1)]
  mcmc_lines <- raw[block2_header_idx:length(raw)]

  if (verbose) {
    message(sprintf(
      "[read_geo_ancstates] File: %s\n  Node block: %d lines | MCMC block: %d lines (1 header + %d iterations)",
      basename(file),
      length(node_lines),
      length(mcmc_lines),
      length(mcmc_lines) - 1L
    ))
  }

  # ---------------------------------------------------------------------------
  # 3. Parse Block 1 — Node metadata
  # ---------------------------------------------------------------------------
  nodes <- .parse_node_block(node_lines)

  # ---------------------------------------------------------------------------
  # 4. Parse Block 2 — MCMC posterior samples
  # ---------------------------------------------------------------------------
  mcmc <- .parse_mcmc_block(mcmc_lines)

  # ---------------------------------------------------------------------------
  # 5. Build tidy long-format coords table
  # ---------------------------------------------------------------------------
  coords <- .build_coords_long(mcmc, nodes)

  # ---------------------------------------------------------------------------
  # 6. Assemble geo_ancstates object
  # ---------------------------------------------------------------------------
  geo <- structure(
    list(
      nodes        = nodes,
      mcmc         = mcmc,
      coords       = coords,
      file         = normalizePath(file),
      n_nodes      = nrow(nodes),
      n_iterations = nrow(mcmc)
    ),
    class = "geo_ancstates"
  )

  # ---------------------------------------------------------------------------
  # 7. Optionally match ape node numbers
  # ---------------------------------------------------------------------------
  if (!is.null(phy)) {
    geo <- add_ape_nodes(geo, phy, verbose = verbose)
  }

  geo
}


# =============================================================================
# add_ape_nodes()
# Match BayesTraits node names to ape node numbers and add as a column.
# Can be called at parse time (via phy argument) or afterwards.
# =============================================================================

#' Match BayesTraits node names to ape node numbers
#'
#' Adds an \code{ape_node} column to \code{geo$nodes} and \code{geo$coords}
#' giving the R node number for each BayesTraits node as it appears in the
#' supplied \code{phylo} object. This lets users refer to nodes by either
#' their BayesTraits name (e.g. \code{"Node-00005"}) or their ape number
#' (e.g. \code{62}).
#'
#' Matching strategy:
#' \itemize{
#'   \item Tip nodes: matched directly by \code{phy$tip.label}.
#'         Tip numbers in ape are \code{1:Ntip(phy)}.
#'   \item Internal nodes: matched via \code{ape::getMRCA()} on the full
#'         set of taxa in each node's taxa list. The ape node number for
#'         internal nodes is \code{Ntip(phy) + 1 : Ntip(phy) + Nnode(phy)}.
#' }
#'
#' @param geo     A \code{geo_ancstates} object from \code{read_geo_ancstates()}.
#' @param phy     A \code{phylo} object (\pkg{ape}). Must be the same tree
#'                used in the Geo model analysis.
#' @param verbose Logical. Report matching statistics (default \code{TRUE}).
#' @return The geo_ancstates object with \code{ape_node} added to
#'         \code{geo$nodes} and \code{geo$coords}.
#'
#' @examples
#' geo <- read_geo_ancstates("Primate.AncStates.txt")
#' phy <- ape::read.nexus("Median.trees")
#'
#' # Option 1 — at parse time
#' geo <- read_geo_ancstates("Primate.AncStates.txt", phy = phy)
#'
#' # Option 2 — after parsing
#' geo <- add_ape_nodes(geo, phy)
#'
#' # Access by BayesTraits name
#' geo$nodes[geo$nodes$node_name == "Node-00008", "ape_node"] # Primate crown node
#'
#' # Access by ape node number
#' geo$nodes[geo$nodes$ape_node == 903, "node_name"]
#'
#' # Filter coords by ape node number
#' geo$coords[geo$coords$ape_node == 903, ]
add_ape_nodes <- function(geo, phy, verbose = TRUE) {

  if (!inherits(geo, "geo_ancstates"))
    stop("'geo' must be a geo_ancstates object from read_geo_ancstates().")
  if (!inherits(phy, "phylo"))
    stop("'phy' must be a phylo object (ape::read.tree / ape::read.nexus).")
  if (!requireNamespace("ape", quietly = TRUE))
    stop("Package 'ape' is required. Install with: install.packages('ape')")

  n_tips <- ape::Ntip(phy)
  nodes  <- geo$nodes

  ape_node <- vapply(seq_len(nrow(nodes)), function(i) {

    taxa   <- nodes$taxa[[i]]
    n_taxa <- nodes$n_taxa[i]

    if (n_taxa == 0) return(NA_integer_)

    if (n_taxa == 1) {
      # Tip node — look up directly in tip.label
      idx <- match(taxa[1], phy$tip.label)
      if (is.na(idx)) {
        warning(sprintf(
          "Tip '%s' (node %s) not found in phy$tip.label — ape_node set to NA.",
          taxa[1], nodes$node_name[i]
        ))
      }
      return(if (is.na(idx)) NA_integer_ else as.integer(idx))
    }

    # Internal node — find MRCA of all taxa present in phy
    present <- taxa[taxa %in% phy$tip.label]
    if (length(present) < 2) {
      warning(sprintf(
        "Node %s: fewer than 2 taxa found in phy$tip.label — ape_node set to NA.",
        nodes$node_name[i]
      ))
      return(NA_integer_)
    }
    as.integer(ape::getMRCA(phy, present))

  }, integer(1))

  # Add to $nodes
  geo$nodes$ape_node <- ape_node

  # Add to $coords via merge
  node_map <- data.frame(
    node_name = nodes$node_name,
    ape_node  = ape_node,
    stringsAsFactors = FALSE
  )
  geo$coords <- merge(geo$coords, node_map, by = "node_name",
                      all.x = TRUE, sort = FALSE)
  geo$coords <- geo$coords[order(geo$coords$iteration,
                                  geo$coords$node_name), ]
  rownames(geo$coords) <- NULL

  if (verbose) {
    n_matched <- sum(!is.na(ape_node))
    n_tips_m  <- sum(!is.na(ape_node) & nodes$n_taxa == 1)
    n_int_m   <- sum(!is.na(ape_node) & nodes$n_taxa  > 1)
    message(sprintf(
      "[add_ape_nodes] Matched %d / %d nodes  (%d tips, %d internal).",
      n_matched, nrow(nodes), n_tips_m, n_int_m
    ))
    if (any(is.na(ape_node)))
      message(sprintf(
        "  %d node(s) could not be matched — check that phy$tip.label matches the BayesTraits taxon names.",
        sum(is.na(ape_node))
      ))
  }

  geo
}


# =============================================================================
# Internal helpers (not exported)
# =============================================================================

.parse_node_block <- function(lines) {

  split_lines <- strsplit(lines, "\t")
  header      <- split_lines[[1]]
  data_rows   <- split_lines[-1]

  result <- lapply(data_rows, function(row) {
    row <- trimws(row)
    row <- row[nzchar(row)]
    n   <- length(row)
    if (n < 4) {
      warning("Skipping malformed node row: ", paste(row, collapse = "|"))
      return(NULL)
    }
    list(
      node_name       = row[1],
      branch_length   = as.numeric(row[2]),
      height          = as.numeric(row[3]),
      restriction_map = row[4],
      taxa            = if (n > 4) row[5:n] else character(0),
      n_taxa          = max(0L, n - 4L)
    )
  })

  result <- Filter(Negate(is.null), result)

  data.frame(
    node_name       = sapply(result, `[[`, "node_name"),
    branch_length   = sapply(result, `[[`, "branch_length"),
    height          = sapply(result, `[[`, "height"),
    restriction_map = sapply(result, `[[`, "restriction_map"),
    n_taxa          = sapply(result, `[[`, "n_taxa"),
    taxa            = I(lapply(result, `[[`, "taxa")),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}


.parse_mcmc_block <- function(lines) {

  header_raw <- strsplit(lines[1], "\t")[[1]]
  header_raw <- trimws(header_raw)
  header_raw <- header_raw[nzchar(header_raw)]

  col_names <- gsub(" - Branch Length$", "_BranchLength", header_raw)
  col_names <- gsub(" - Long$",          "_Long",         col_names)
  col_names <- gsub(" - Lat$",           "_Lat",          col_names)
  col_names <- gsub("-",                 "_",              col_names)
  col_names <- gsub(" ",                 "",               col_names)

  data_lines <- sub("\\t+$", "", lines[-1])

  con <- textConnection(data_lines)
  on.exit(close(con), add = TRUE)

  mcmc <- tryCatch(
    read.table(
      con,
      sep              = "\t",
      header           = FALSE,
      col.names        = col_names,
      colClasses       = "numeric",
      stringsAsFactors = FALSE,
      comment.char     = ""
    ),
    error = function(e) {
      stop(
        "Failed to parse MCMC block.\n",
        "  Original error: ", conditionMessage(e), "\n",
        "  Check that all MCMC rows have the same number of columns as the header."
      )
    }
  )

  mcmc
}


.build_coords_long <- function(mcmc, nodes) {

  long_cols <- grep("_Long$", names(mcmc), value = TRUE)
  node_ids  <- sub("_Long$", "", long_cols)

  n_iter      <- nrow(mcmc)
  n_nodes     <- length(node_ids)
  total_rows  <- n_iter * n_nodes

  iteration     <- integer(total_rows)
  node_name_out <- character(total_rows)
  long_out      <- numeric(total_rows)
  lat_out       <- numeric(total_rows)
  bl_out        <- numeric(total_rows)

  idx <- 1L
  for (nd in node_ids) {
    long_col <- paste0(nd, "_Long")
    lat_col  <- paste0(nd, "_Lat")
    bl_col   <- paste0(nd, "_BranchLength")
    rows     <- idx:(idx + n_iter - 1L)
    bt_name  <- sub("_", "-", nd)

    iteration[rows]     <- mcmc$Itter
    node_name_out[rows] <- bt_name
    long_out[rows]      <- mcmc[[long_col]]
    lat_out[rows]       <- mcmc[[lat_col]]
    bl_out[rows]        <- if (bl_col %in% names(mcmc)) mcmc[[bl_col]] else NA_real_

    idx <- idx + n_iter
  }

  coords <- data.frame(
    iteration            = iteration,
    node_name            = node_name_out,
    longitude            = long_out,
    latitude             = lat_out,
    branch_length_scaled = bl_out,
    stringsAsFactors     = FALSE
  )

  meta          <- nodes[, c("node_name", "height", "n_taxa", "restriction_map")]
  meta$is_tip   <- meta$n_taxa == 1L
  coords        <- merge(coords, meta, by = "node_name", all.x = TRUE, sort = FALSE)
  coords        <- coords[order(coords$iteration, coords$node_name), ]
  rownames(coords) <- NULL
  coords
}


# =============================================================================
# S3 methods
# =============================================================================

print.geo_ancstates <- function(x, ...) {
  cat("BayesTraits Geo Model — AncStates\n")
  cat("  File      :", x$file, "\n")
  cat("  Nodes     :", x$n_nodes,
      sprintf("(%d tips, %d internal)",
              sum(x$nodes$n_taxa == 1),
              sum(x$nodes$n_taxa != 1)), "\n")
  cat("  Iterations:", x$n_iterations, "\n")
  cat("  Heights   : min =", round(min(x$nodes$height), 4),
      "| max =", round(max(x$nodes$height), 4), "\n")
  has_ape <- "ape_node" %in% names(x$nodes)
  cat("  ape nodes :", if (has_ape) "yes (ape_node column present)" else
                         "no  (supply phy to add_ape_nodes())", "\n")
  cat("\nSlots: $nodes | $mcmc | $coords\n")
  invisible(x)
}

summary.geo_ancstates <- function(object, ...) {
  cat("=== geo_ancstates summary ===\n\n")

  node_cols <- c("node_name", "branch_length", "height", "restriction_map", "n_taxa")
  if ("ape_node" %in% names(object$nodes))
    node_cols <- c(node_cols, "ape_node")

  cat("-- Node table ($nodes) --\n")
  print(head(object$nodes[, node_cols], 6))
  cat("  ... [", nrow(object$nodes), "rows total ]\n\n")

  cat("-- MCMC iterations ($mcmc) --\n")
  cat("  Iterations :", nrow(object$mcmc), "\n")
  cat("  Lh range   :", round(range(object$mcmc$Lh), 2), "\n")
  cat("  Scale range:", round(range(object$mcmc$Scale), 4), "\n\n")

  cat("-- Tidy coords ($coords) --\n")
  cat("  Rows       :", nrow(object$coords), "\n")
  cat("  Long range :", round(range(object$coords$longitude, na.rm = TRUE), 3), "\n")
  cat("  Lat range  :", round(range(object$coords$latitude,  na.rm = TRUE), 3), "\n")
  invisible(object)
}
