# ----------------------------------------------------------------------------
# Title: Spatial covariates, PCA construction, and mesh generation for INLA
# Author: <Your name>
# Date: <YYYY-MM-DD>
# Description:
#   This script prepares landscape covariates, computes a PCA to derive
#   composite axes, and builds spatial meshes (global and road-buffer meshes)
#   for point process models (e.g., INLA/inlabru). The workflow targets a
#   regional case study in western France and produces rasterized PCA axes and
#   auxiliary covariates (e.g., road traffic, hydrology). All operations are
#   performed in a Lambert Conformal Conic system equivalent to Lambert-93,
#   with coordinates expressed in kilometers for convenience.
#
# Reproducibility tips for publication:
#   - Record package versions with sessionInfo() at the end of the analysis.
#   - Set a deterministic seed if any stochastic step is added.
#   - Provide data provenance (sources, resolution, year, processing).
#   - Consider archiving input rasters/shapefiles or providing a data DOI.
# ----------------------------------------------------------------------------

# -------------------------------
# Required packages (versions)
# -------------------------------
# sf (>= 1.0), raster (>= 3.6), terra (>= 1.7), ade4, INLA, inlabru/fmesher
# NOTE: The code mixes 'raster' and 'terra'. This is fine if you are careful
#       with conversions (raster<->terra) and object classes.
# library(sf)
# library(raster)
# library(terra)
# library(ade4)
# library(INLA)      # or R-INLA
# library(inlabru)   # for fm_mesh_2d, fm_contains, pixels, projectors, etc.

# -------------------------------
# Coordinate reference systems
# -------------------------------
# PROJ_m : Lambert Conformal Conic (Lambert-93 equivalent) with metric units
# PROJ_km: Same parameters but units in kilometers (convenience for window sizes)
PROJ_m  = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
PROJ_km = "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"

# ----------------------------------------------------------------------------
# 1) Study domain (Brittany departments)
# ----------------------------------------------------------------------------
# Read French departments, transform to km-based CRS, then subset to
# Finistère (29), Côtes-d'Armor (22), Morbihan (56), Ille-et-Vilaine (35).
# NOTE: If possible, reference an authoritative source and version date for
#       the shapefile. Consider using st_read options(stringsAsFactors=FALSE).

domaine = st_read("data/vector/Dep_france_wgs84.shp")

# Work in kilometers to make window sizes (e.g., 11x11, 21x21) interpretable.
domaine = st_transform(domaine, PROJ_km)

# Keep only Brittany departments; non-numeric codes become 99 and are dropped.
numDep = as.numeric(domaine$CODE_DEPT)
numDep[is.na(numDep)] = 99

# Subset to the four departments and dissolve to one geometry
# (union + polygon cast). The [5,] selection assumes a particular order
# after st_cast; this can be brittle—consider selecting by area or attribute.
domaine = domaine[numDep %in% c(29, 22, 56, 35), ]
domaine = st_union(domaine)
domaine = st_cast(domaine, "POLYGON")
domaine = domaine[5, ]   # TODO: ensure this selects the intended polygon

# Return to sf with CRS attached
domaine = st_as_sf(domaine)

# ----------------------------------------------------------------------------
# 2) Covariates for PCA (land cover & linear features at 1 km)
# ----------------------------------------------------------------------------
# Each raster is:
#   - read from disk
#   - extent is rescaled by 1/1000 to match km units (assumes original was m)
#   - projection text is updated from units=m to units=km
#   - standardized (scale) and masked to the study domain
# CAUTION: Directly rewriting the extent and CRS string assumes the only
#          difference is the unit. When available, projecting with projectRaster
#          (raster) or project() (terra) is safer/clearer.

foretR  = raster("data/raster/foret_landcover_1km.tif")
extent(foretR)      = extent(c(xmin(foretR), xmax(foretR), ymin(foretR), ymax(foretR)) / 1000)
projection(foretR)  = gsub("units=m", "units=km", projection(foretR))
foretR = mask(scale(foretR), domaine)

cultureR = raster("data/raster/culture_landcover_1km.tif")
extent(cultureR)     = extent(c(xmin(cultureR), xmax(cultureR), ymin(cultureR), ymax(cultureR)) / 1000)
projection(cultureR) = gsub("units=m", "units=km", projection(cultureR))
cultureR = mask(scale(cultureR), domaine)

haieR = raster("data/raster/densite_haie_1km.tif")
extent(haieR)     = extent(c(xmin(haieR), xmax(haieR), ymin(haieR), ymax(haieR)) / 1000)
projection(haieR) = gsub("units=m", "units=km", projection(haieR))
haieR = mask(scale(haieR), domaine)

hydroR = raster("data/raster/densite_hydro_1km.tif")
extent(hydroR)     = extent(c(xmin(hydroR), xmax(hydroR), ymin(hydroR), ymax(hydroR)) / 1000)
projection(hydroR) = gsub("units=m", "units=km", projection(hydroR))
hydroR = mask(scale(hydroR), domaine)

urbainR = raster("data/raster/urbain_landcover_1km.tif")
extent(urbainR)     = extent(c(xmin(urbainR), xmax(urbainR), ymin(urbainR), ymax(urbainR)) / 1000)
projection(urbainR) = gsub("units=m", "units=km", projection(urbainR))
urbainR = mask(scale(urbainR), domaine)

# Assemble PCA input table from standardized rasters; remove rows with NA
# across all covariates to avoid dropping cells inconsistently.
datACP = data.frame(
  foret  = values(foretR),
  culture= values(cultureR),
  hydro  = values(hydroR),
  haie   = values(haieR),
  urbain = values(urbainR)
)
rownames(datACP) = paste("L", 1:nrow(datACP))
datACP = datACP[!is.na(apply(datACP, 1, sum)), ]

# Principal Component Analysis (ade4)
# Data are already standardized, so center=FALSE, scale=FALSE is appropriate.
pc = dudi.pca(datACP, scannf = FALSE, nf = 4, center = FALSE, scale = FALSE)

# Optionally visualize variable loadings on axes 3 & 4 with s.arrow (commented)
# s.arrow(pc$c1, xax = 3, yax = 4, clabel = 1.2,
#         sub = "Représentation des covariables", csub = 1.5)

# Attach site scores (pc$li) to the corresponding non-NA rows, then expand back
# to the full raster grid so we can map each PCA axis as a raster layer.
datACP = cbind(datACP, pc$li)

datRASTER = data.frame(
  foret  = values(foretR),
  culture= values(cultureR),
  hydro  = values(hydroR),
  haie   = values(haieR),
  urbain = values(urbainR)
)
rownames(datRASTER) = paste("L", 1:nrow(datRASTER))

datNEW = matrix(NA, ncol = ncol(datACP), nrow = nrow(datRASTER))
linonNA = (1:nrow(datRASTER))[!is.na(apply(datRASTER, 1, sum))]
datNEW[linonNA, ] = as.matrix(datACP)
colnames(datNEW) = colnames(datACP)

# Rasterize PCA axes by copying the foretR template and assigning cell values
Axis1R = foretR; values(Axis1R) = datNEW[, "Axis1"]
Axis2R = foretR; values(Axis2R) = datNEW[, "Axis2"]
Axis3R = foretR; values(Axis3R) = datNEW[, "Axis3"]
Axis4R = foretR; values(Axis4R) = datNEW[, "Axis4"]

# ----------------------------------------------------------------------------
# 3) Additional covariates (sampling, speed/traffic, distances)
# ----------------------------------------------------------------------------
# 'sampling' provides point-level PCA scores (Axis1/Axis2) to be spatially
# averaged on a focal window, producing smoothed fields.

sampling = st_read("data/sampling/new_sampling_raw.shp")
sampling = st_transform(sampling, PROJ_km)

# Template raster for sampling-based layers (Axis1/Axis2). As above, extent &
# projection are rewritten to km units; see caution note in Section 2.
samplingR = raster("data/raster/sampling_axis1_1km.tif")
extent(samplingR)    = extent(c(xmin(samplingR), xmax(samplingR), ymin(samplingR), ymax(samplingR)) / 1000)
projection(samplingR) = gsub("units=m", "units=km", projection(samplingR))
samplingR = mask(samplingR, domaine)

# Switch to terra for rasterize operations; convert sf->SpatVector, Raster->SpatRaster
samplingR      = rast(samplingR)
points_vect    = vect(sampling)
points_rasterA1 = rasterize(points_vect, samplingR, field = "Axis1")
points_rasterA2 = rasterize(points_vect, samplingR, field = "Axis2")

# Combine layers, then extract the two that hold Axis1/Axis2 fields
raster_with_new_layer = c(samplingR, points_rasterA1)
raster_with_new_layer = c(raster_with_new_layer, points_rasterA2)
Axis1 = raster_with_new_layer[[2]]
Axis2 = raster_with_new_layer[[3]]

# Back to 'raster' package for focal smoothing on a square window.
# Window sizes (21x21) are in pixels; given 1-km grid and km-based CRS,
# this corresponds to ~21 km moving average.
raster_layerA1 = raster::raster(Axis1)
raster_layerA2 = raster::raster(Axis2)
samplingRA1 = focal(raster_layerA1, w = matrix(1, 21, 21), fun = mean, pad = TRUE, na.rm = TRUE)
samplingRA2 = focal(raster_layerA2, w = matrix(1, 21, 21), fun = mean, pad = TRUE, na.rm = TRUE)

# Road speed (vitesse) and traffic covariates; zeros set to NA, then smoothed
# with an 11x11 (~11 km) mean filter, standardized, and masked.
vitesseR = raster("data/raster/vitesse_1km.tif")
extent(vitesseR)     = extent(c(xmin(vitesseR), xmax(vitesseR), ymin(vitesseR), ymax(vitesseR)) / 1000)
projection(vitesseR) = gsub("units=m", "units=km", projection(vitesseR))
vitesseR[vitesseR == 0] = NA
vitesseR = focal(vitesseR, w = matrix(1, 11, 11), fun = mean, pad = TRUE, na.rm = TRUE)
vitesseR = mask(scale(vitesseR), domaine)

traficR = raster("data/raster/trafic_1km.tif")
extent(traficR)     = extent(c(xmin(traficR), xmax(traficR), ymin(traficR), ymax(traficR)) / 1000)
projection(traficR) = gsub("units=m", "units=km", projection(traficR))
traficR[traficR == 0] = NA
traficR = focal(traficR, w = matrix(1, 11, 11), fun = mean, pad = TRUE, na.rm = TRUE)
traficR = mask(scale(traficR), domaine)

# Distance rasters (hydrography, hedgerows). Apply log(x+1) to dampen
# long tails, then standardize and mask.
hydrodistR = raster("data/raster/raster_hydro_dist_50m.tif")
extent(hydrodistR)     = extent(c(xmin(hydrodistR), xmax(hydrodistR), ymin(hydrodistR), ymax(hydrodistR)) / 1000)
projection(hydrodistR) = gsub("units=m", "units=km", projection(hydrodistR))
hydrodistR = mask(scale(log(hydrodistR + 1)), domaine)

haiedistR = raster("data/raster/raster_veg_dist_1km.tif")
extent(haiedistR)     = extent(c(xmin(haiedistR), xmax(haiedistR), ymin(haiedistR), ymax(haiedistR)) / 1000)
projection(haiedistR) = gsub("units=m", "units=km", projection(haiedistR))
haiedistR = mask(scale(log(haiedistR + 1)), domaine)

# ----------------------------------------------------------------------------
# 4) Mesh construction (INLA / inlabru)
# ----------------------------------------------------------------------------
# Load point data for collision (roadkill) and live observations; convert to sf,
# bring to km-based CRS, and clip to the domain.
load("data/donBON.rda")
collision = roadkill_sp
sfcoll = sf::st_as_sf(collision, coords = c("X", "Y"))
st_crs(sfcoll) = PROJ_m
sfcoll = st_transform(sfcoll, PROJ_km)
sfcoll = st_intersection(sfcoll, domaine)

vivant = datVIV
sfviv = sf::st_as_sf(vivant, coords = c("x", "y"))
st_crs(sfviv) = PROJ_m
sfviv = st_transform(sfviv, PROJ_km)
sfviv = st_intersection(sfviv, domaine)

# Build inner/outer non-convex hulls around domain + observation coordinates
# to control mesh boundaries; negative convex values yield slight erosions.
coordMESH = rbind(as.matrix(st_coordinates(domaine)[, 1:2]),
                  as.matrix(st_coordinates(sfviv)[, 1:2]))
bndint = inla.nonconvex.hull(coordMESH, convex = -0.03)
bndext = inla.nonconvex.hull(coordMESH, convex = -0.2)

# Main mesh (coarser far from data, finer near data); cutoff avoids near-duplicate
# points; max.edge controls triangle sizes (first for inner, second for outer).
mesh = inla.mesh.2d(
  loc      = rbind(as.matrix(st_coordinates(sfviv)[, 1:2]),
                   as.matrix(st_coordinates(sfcoll)[, 1:2])),
  boundary = list(bndint, bndext),
  max.edge = c(0.5, 15),
  cutoff   = c(1.5),
  crs      = PROJ_km
)

# Voronoi (dual) mesh cells for integration weights, converted to sf
dmesh = book.mesh.dual(mesh, Ncores = 20)
dmesh = st_as_sf(dmesh)
st_crs(dmesh) = PROJ_km

# ----------------------------------------------------------------------------
# 5) Road-buffer mesh (for mortality around roads)
# ----------------------------------------------------------------------------
# Read road network, transform, and restrict to domain. Build a 2-km outer
# buffer, then shrink by 1 km to create an annulus (routeBUF) around roads.
# routeBAR = domain minus that annulus (where we will place mesh triangles).
routes = st_read("data/routes_diro_fill.shp")
routes = st_transform(routes, PROJ_km)
routes = st_intersection(routes, domaine)
routeBUFtmp = st_buffer(routes, dist = 2)
routeBUFtmp = st_union(routeBUFtmp)
routeBUF    = st_buffer(routeBUFtmp, dist = -1)
routeBUF    = st_sf(geometry = routeBUF)
routeBAR    = st_difference(domaine, routeBUF)

# Build a finer mesh over the full domain; then identify triangles whose
# centroids lie in routeBAR. Edge/offset are fractions of the domain range.
r = mean(apply(as.matrix(st_coordinates(domaine)[, 1:2]), 2, max) -
         apply(as.matrix(st_coordinates(domaine)[, 1:2]), 2, min))
meshROUTE = fm_mesh_2d(
  loc.domain = as.matrix(st_coordinates(domaine)[, 1:2]),
  max.edge   = c(0.008, 0.1) * r,
  offset     = c(0.01, 0.3) * r,
  cutoff     = 0.005 * r,
  crs        = PROJ_km
)

route.tri = unlist(fm_contains(routeBAR, y = meshROUTE, type = "centroid"))

# Voronoi dual for the road-buffer mesh
dmeshROUTE = book.mesh.dual(meshROUTE, Ncores = 20)
dmeshROUTE = st_as_sf(dmeshROUTE)
st_crs(dmeshROUTE) = PROJ_km

# ----------------------------------------------------------------------------
# 6) Prediction grids and projectors
# ----------------------------------------------------------------------------
# Create pixel grids for (i) live observations model (domain) and (ii) road
# mortality model (route buffer). Compute centroids for projector construction.
# The 'resopred' parameter controls the raster dimensions (~96 x aspect).

toto   = extent(domaine)
xrange = toto[2] - toto[1]
yrange = toto[4] - toto[3]

resopred = 96
pxlVIV = pixels(mesh, mask = domaine, nx = resopred, ny = round(resopred * yrange / xrange))
PROJgridVIV = inla.mesh.projector(mesh, st_coordinates(st_centroid(pxlVIV$pixels_sf)))

resopred = 96
pxlMOR = pixels(meshROUTE, mask = routeBUF, nx = resopred, ny = round(resopred * yrange / xrange))
# Area of intersection with the route buffer (integration weights)
pxlMOR$area_inter = as.vector(st_area(st_intersection(pxlMOR$pixels_sf, routeBUF)))
PROJgridMOR = inla.mesh.projector(meshROUTE, st_coordinates(st_centroid(pxlMOR$pixels_sf)))

# ----------------------------------------------------------------------------
# 7) Species grouping helper
# ----------------------------------------------------------------------------
# Build a mapping between species names present in roadkill and live datasets,
# then collapse selected species into a 'mustelidae' group for analyses.

listeSP = intersect(unique(roadkill_sp$spNEW), unique(datVIV$nom_vern))
grpSP   = listeSP
grpSP[c(3, 7, 8, 9, 11)] = "mustelidae"   # document the indices/species in MS
listeSP = cbind(grpSP, listeSP)

# ----------------------------------------------------------------------------
# 8) Cleanup
# ----------------------------------------------------------------------------
# Free memory by removing large intermediate objects and triggering GC.
rm("pc", "datRASTER", "datACP", "datNEW", "datVIV",
   "roadkill_sp", "vivant", "routes", "sfviv", "sfcoll")
gc()

# ----------------------------------------------------------------------------
# Optional: capture session info for reproducibility in the manuscript
# ----------------------------------------------------------------------------
# sessionInfo()
