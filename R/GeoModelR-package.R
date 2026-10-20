#' GeoModelR: tools for the BayesTraits Geo model
#'
#' Prepare palaeomaps and geographic input data for the Geo model in
#' BayesTraits v5+, and explore its output: ancestral locations on
#' age-matched palaeomaps and posterior dispersal distances.
#'
#' Main functions, by workflow:
#' \itemize{
#'   \item Input data: \code{\link{clean_sister_geodata}()}.
#'   \item Palaeomaps: \code{\link{get_paleomaps}()},
#'         \code{\link{fill_missing_maps}()}, \code{\link{sf_to_csv}()},
#'         \code{\link{write_buildmaps_file}()}.
#'   \item Output: \code{\link{read_geo_ancstates}()},
#'         \code{\link{add_ape_nodes}()}, \code{\link{plot_node}()},
#'         \code{\link{plot_timeinterval}()}, \code{\link{plot_clade}()},
#'         \code{\link{get_dispersal_distance}()}.
#' }
#'
#' Example data (a primate tree and a Geo model AncStates file with 50
#' posterior samples) are installed in \code{system.file("extdata",
#' package = "GeoModelR")}; see \code{vignette("GeoModelR")} for a
#' walkthrough.
#'
#' @importFrom stats median runif setNames
#' @importFrom utils head menu read.table write.csv write.table
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"

# Column names used inside ggplot2::aes() calls (non-standard evaluation).
utils::globalVariables(c(
  "longitude", "latitude", "is_tip", "height", "node_name", "node"
))
