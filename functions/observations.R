# ---------------------------------------------------------------------------
# Helper function: chargeDON()
# Purpose: Load collision and live-observation data for selected species and a
#          time window; harmonize CRSs; clip to the study domain; and compute
#          per-pixel counts on predefined prediction grids.
# ---------------------------------------------------------------------------

#' Load and prepare species-specific data within a time window
#'
#' This helper reads objects from a .rda file, filters by species and year, and
#' returns (i) sf point layers (live observations and collisions) in a
#' kilometer-based LCC CRS, and (ii) per-pixel counts for two prediction grids.
#'
#' @param espece character vector of species names to retain (must match the
#'   content of roadkill_sp$spNEW and datVIV$nom_vern).
#' @param Tmin numeric/integer lower year bound (exclusive in this implementation).
#' @param Tmax numeric/integer upper year bound (exclusive in this implementation).
#'   # NOTE: If inclusive bounds are desired, replace the filters with
#'   #       annee >= Tmin and annee <= Tmax.
#' @param domaine sf POLYGON/MULTIPOLYGON defining the study area in km-based LCC.
#' @param file character path to an .rda file that contains objects named
#'   `roadkill_sp` and `datVIV` with coordinates as columns X,Y (collisions)
#'   and x,y (live), respectively, in meter-based LCC (PROJ_m).
#'
#' @details
#' The function assumes that two pixel grids already exist in the calling
#' environment: `pxlVIV$pixels_sf` (for live observations) and
#' `pxlMOR$pixels_sf` (for collisions/mortality). These grids are typically
#' produced from meshes (see main script). For manuscript reproducibility,
#' consider passing these grids explicitly as function arguments instead of
#' relying on globals.
#'
#' @return A named list with elements:
#'   - sfviv   : sf POINTS of live observations (filtered & clipped)
#'   - sfcoll  : sf POINTS of collisions (filtered & clipped)
#'   - NobsVIV : integer vector, counts per pixel in `pxlVIV$pixels_sf`
#'   - NobsMOR : integer vector, counts per pixel in `pxlMOR$pixels_sf`
#'
#' @seealso inla.mesh.2d, pixels(), st_intersects()
#'
chargeDON = function(espece, Tmin, Tmax, domaine, file) {
  # -------------------------------
  # Coordinate reference systems
  # -------------------------------
  PROJ_m  = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  PROJ_km = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"

  # -------------------------------
  # Input validation (lightweight)
  # -------------------------------
  stopifnot(is.character(espece), length(espece) >= 1)
  stopifnot(inherits(domaine, "sf"))
  if (!all(sf::st_is(domaine, c("POLYGON", "MULTIPOLYGON")))) {
    warning("`domaine` is not a polygonal sf object; attempting st_cast().")
    domaine <- sf::st_cast(domaine, "MULTIPOLYGON")
  }
  if (!file.exists(file)) stop("`file` does not exist: ", file)

  # -------------------------------
  # Load data objects (expects roadkill_sp, datVIV)
  # -------------------------------
  loaded <- load(file)
  if (!all(c("roadkill_sp", "datVIV") %in% loaded | c("roadkill_sp", "datVIV") %in% ls())) {
    stop("`file` must provide objects named `roadkill_sp` and `datVIV`.")
  }

  # -------------------------------
  # Filter collisions by years and species; build sf layer
  # -------------------------------
  collision <- roadkill_sp
  collision <- collision[(collision$annee > Tmin) & (collision$annee < Tmax), ]  # exclusive bounds
  collision <- collision[collision$spNEW %in% espece, ]
  sfcoll <- sf::st_as_sf(collision, coords = c("X", "Y"))
  sf::st_crs(sfcoll) <- PROJ_m
  sfcoll <- sf::st_transform(sfcoll, PROJ_km)
  # Clip to study domain; st_make_valid can help if geometries are invalid
  sfcoll <- sf::st_intersection(sfcoll, domaine)

  # -------------------------------
  # Filter live observations by years and species; build sf layer
  # -------------------------------
  vivant <- datVIV
  vivant <- vivant[(vivant$annee > Tmin) & (vivant$annee < Tmax), ]  # exclusive bounds
  vivant <- vivant[vivant$nom_vern %in% espece, ]
  sfviv <- sf::st_as_sf(vivant, coords = c("x", "y"))
  sf::st_crs(sfviv) <- PROJ_m
  sfviv <- sf::st_transform(sfviv, PROJ_km)
  sfviv <- sf::st_intersection(sfviv, domaine)

  # -------------------------------
  # Per-pixel counts on predefined grids
  # -------------------------------
  # Expect these to exist in the parent environment (see main script section 6).
  if (!exists("pxlVIV", inherits = TRUE) || !exists("pxlMOR", inherits = TRUE)) {
    stop("`pxlVIV` and `pxlMOR` must exist (with $pixels_sf) in the calling environment.")
  }

  # Live observations per pixel
  rSF <- pxlVIV$pixels_sf
  test <- sf::st_intersects(rSF, sfviv, sparse = FALSE)
  NobsVIV <- rowSums(test)

  # Collisions per pixel
  rSF <- pxlMOR$pixels_sf
  test <- sf::st_intersects(rSF, sfcoll, sparse = FALSE)
  NobsMOR <- rowSums(test)

  # Named return for clarity
  return(list(
    sfviv   = sfviv,
    sfcoll  = sfcoll,
    NobsVIV = NobsVIV,
    NobsMOR = NobsMOR
  ))
}
