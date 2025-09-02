# pixels_commented.R
# -----------------------------------------------------------------------------
# Title: Create a regular grid of square pixels over an INLA mesh domain
# Author: <Your name>
# Date: <YYYY-MM-DD>
#
# Description
# -----------
# This function builds a regular grid of square polygons ("pixels") spanning
# the spatial extent of a 2D INLA mesh. It optionally masks the grid to a
# supplied polygonal region (e.g., study area) or to a domain derived from the
# mesh itself. The output is an sf object of pixel polygons and the (constant)
# pixel area, which can be used for projection, integration, or mapping.
#
# Provenance
# ----------
# Adapted from the original *pixels* helpers used in SPDE workflows (INLA);
# this version adds publication-grade comments, minor robustness tweaks, and a
# safer masking step for multi-feature domains.
#
# Dependencies: sf (>= 1.0)
# -----------------------------------------------------------------------------

#' Create square pixel polygons from an INLA mesh
#'
#' @param mesh An `inla.mesh` object (manifold 'R2'). The function uses
#'   `mesh$loc` to determine x/y ranges. If the mesh carries a CRS in
#'   `mesh$crs$input`, it is propagated to the output if possible.
#' @param nx,ny Either a single integer (number of pixels along x/y; covering
#'   the min-max mesh extent) or a numeric vector of explicit x/y breakpoints.
#'   If vectors are supplied, they are sorted and the minimum positive spacing
#'   is used as the grid cell size in each direction.
#' @param mask Logical or `sf` object. If `TRUE`, the grid is masked by a
#'   polygon derived from the mesh (attempts `st_as_sf(mesh)`, then falls back
#'   to the convex hull of `mesh$loc`). If an `sf`/`sfc` object is supplied,
#'   the grid is masked by intersection with that geometry. If `FALSE`, no
#'   masking is applied.
#'
#' @return A list with two elements:
#'   - `pixels_sf`: an `sf` object of square polygons
#'   - `area_pix` : numeric, the area of one pixel (constant across the grid)
#'
#' @examples
#' # px <- pixels(mesh, nx = 100, ny = 100, mask = TRUE)
#' # plot(px$pixels_sf[1:100, ], col = 'grey80', border = NA)
#'
#' @note For irregular domains, masking with `st_intersects` (used below) keeps
#'   any pixel that intersects the mask polygon(s). If you want only pixels
#'   fully contained within the mask, switch to `st_within`.
#'
#' @import sf
#' @export
pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE) {
  # Ensure sf is available
  if (!requireNamespace("sf", quietly = TRUE))
    stop("The 'sf' package is required but not installed.")

  # -------------------------------
  # Build x/y coordinate sequences
  # -------------------------------
  # If nx/ny are scalars, create equally spaced sequences spanning mesh extent.
  # If vectors, sort and use them directly.
  if (length(nx) == 1) {
    x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length.out = nx)
  } else {
    x <- sort(as.numeric(nx))
  }
  if (length(ny) == 1) {
    y <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length.out = ny)
  } else {
    y <- sort(as.numeric(ny))
  }

  # -------------------------------
  # Compute cell size (handles irregular spacing by using min positive diff)
  # -------------------------------
  dx <- if (length(x) > 1) min(diff(x)) else 0
  dy <- if (length(y) > 1) min(diff(y)) else 0
  if (dx <= 0 || dy <= 0) stop("Cannot determine positive grid cell size.")

  # -------------------------------
  # Create grid of squares with st_make_grid
  # -------------------------------
  # Build a bbox and attach CRS if available on the mesh
  crs_out <- tryCatch(mesh$crs$input, error = function(e) NULL)
  bb <- sf::st_bbox(c(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)))
  if (!is.null(crs_out)) attr(bb, "crs") <- crs_out

  grid_sf <- sf::st_make_grid(sf::st_as_sfc(bb), cellsize = c(dx, dy), what = "polygons")

  # -------------------------------
  # Optional masking
  # -------------------------------
  if (is.logical(mask)) {
    if (isTRUE(mask)) {
      # Try to convert mesh to sf; if it fails, use convex hull of mesh nodes
      dom_sf <- tryCatch(sf::st_as_sf(mesh), error = function(e) NULL)
      if (is.null(dom_sf)) {
        pts <- mesh$loc[, 1:2, drop = FALSE]
        hull <- pts[chull(pts), , drop = FALSE]
        dom_sf <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(hull))), crs = crs_out)
      }
      keep <- sf::st_intersects(grid_sf, dom_sf, sparse = FALSE)
      # keep is a logical matrix if dom_sf has multiple features
      keep_rows <- apply(keep, 1, any)
      grid_sf <- grid_sf[keep_rows]
    }
  } else if (inherits(mask, "sf") || inherits(mask, "sfc")) {
    keep <- sf::st_intersects(grid_sf, mask, sparse = FALSE)
    keep_rows <- apply(keep, 1, any)
    grid_sf <- grid_sf[keep_rows]
  } # else: unrecognized mask type -> no masking

  # Pixel area (constant for a regular grid)
  area_pix <- dx * dy

  # Return as requested
  return(list(pixels_sf = sf::st_sf(geometry = grid_sf), area_pix = area_pix))
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
