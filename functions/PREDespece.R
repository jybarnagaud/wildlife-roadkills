# prediction_commented.R
# -----------------------------------------------------------------------------
# Title: Prediction workflow (held-out buffer validation & barrier SPDE model)
# Author: <Your name>
# Date: <YYYY-MM-DD>
#
# Description
# -----------
# For each species group, this script:
#   1) Fits a live-fauna spatial point-process model (LGCP style) on a domain
#      mesh to estimate the latent intensity surface (with PC priors on the
#      Matérn field hyperparameters). It refits with fixed hyperparameters and
#      keeps a single posterior sample of the latent field for prediction.
#   2) Repeats Nrep times: draws a random circular buffer within the domain,
#      removes collision points within the buffer, (re)computes integration
#      weights for the route-buffer mesh excluding the buffer, fits a barrier
#      SPDE model for road mortality on the remaining data, and predicts the
#      mortality intensity inside the held-out buffer. Results are saved per
#      repetition to assess out-of-sample predictive capacity.
#
# Notes
# -----
# - Coordinates and meshes are prepared upstream (covariables_mesh.R).
# - `chargeDON()` loads/filters observations and returns counts per prediction
#   pixels; it relies on global pixel grids (`pxlVIV`, `pxlMOR`).
# - The live-fauna stage uses one posterior sample (`inla.posterior.sample(1)`)
#   to build a covariate raster for the road stage. Adjust as needed if you
#   prefer to propagate uncertainty by taking >1 samples (as in main fitting).
# - Random buffer placement uses `st_sample(domaine, size = 1)`. For strict
#   reproducibility, set a seed before the foreach loop.
# -----------------------------------------------------------------------------

# -------------------------------
# 0) Setup
# -------------------------------
setwd("/mnt/c/Users/jpapaix/Documents/modelisation/mortaliteroutiere/codes_propres")

library(raster)
library(ade4)
library(sf)
library(INLA)
library(terra)
library(stars)
library(doParallel)
library(foreach)
library(qs)
library(dplyr)
library(ggplot2)
inla.setOption(num.threads = 1)       # Global default; each inla() call may override

# Helpers (utilities & grid generator)
source("inlabookfunctions.R")
source("pixels.R")

# Spatial modelling extras
library(INLAspacetime)
library(inlabru)

# -------------------------------
# 1) Covariates, meshes, observations helpers
# -------------------------------
source("covariables_mesh.R")  # builds meshes, projectors, covariate rasters
source("observations.R")      # provides chargeDON()

# -------------------------------
# 2) Loop over species groups
# -------------------------------
for (gsp in unique(listeSP[, "grpSP"])) {  # e.g., "Chevreuil europeen"
  # Species names mapped to the current group
  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]

  # Load observations (2015–2019 inclusive/exclusive per chargeDON's filter)
  DON     <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  datviv  <- DON[[1]]  # sf POINTS: live fauna
  datcoll <- DON[[2]]  # sf POINTS: collisions
  NobsVIV <- DON[[3]]  # counts per pxlVIV pixel (for diagnostics)
  NobsMOR <- DON[[4]]  # counts per pxlMOR pixel (for diagnostics)

  # -----------------------------
  # 3) Live fauna (LGCP-style) model
  # -----------------------------
  # Extract covariates at mesh node locations and presence points by stacking
  coords <- rbind(mesh$loc[, 1:2], st_coordinates(datviv))

  CAxis1      <- extract(Axis1R,     coords)
  CAxis2      <- extract(Axis2R,     coords)
  CAxis3      <- extract(Axis3R,     coords)
  CAxis4      <- extract(Axis4R,     coords)
  CsamplingA1 <- extract(samplingRA1, coords)
  CsamplingA2 <- extract(samplingRA2, coords)  # available if needed later

  # --- 3.1 SPDE prior (PC priors) and precision ---
  matern <- inla.spde2.pcmatern(
    mesh,
    alpha = 2,                       # SPDE order (nu = 1)
    prior.sigma = c(1.6, 0.1),       # Example: P(sigma > 1.6) = 0.1
    prior.range = c(55, 0.9)         # Example: P(range < 55) = 0.9
  )

  # Integration weights (areas of dual cells intersecting the domain)
  domaine_outside_buffer <- domaine  # placeholder for possible future exclusion
  intersect_idx <- st_intersects(dmesh, domaine_outside_buffer, sparse = FALSE)
  intersections <- st_intersection(dmesh[intersect_idx, ], domaine_outside_buffer)
  w <- numeric(nrow(dmesh))
  w[intersect_idx] <- as.numeric(st_area(intersections))

  # Build point-process stack (mesh integration points + events)
  nv  <- mesh$n
  n   <- nrow(datviv)
  y.pp <- rep(0:1, c(nv, n))
  e.pp <- c(w, rep(0, n))
  imat <- Diagonal(nv, rep(1, nv))
  lmat <- inla.spde.make.A(mesh, st_coordinates(datviv))
  A.pp <- rbind(imat, lmat)

  stk.pp <- inla.stack(
    data = list(y = y.pp, e = e.pp),
    A    = list(1, A.pp),
    effects = list(
      list(
        b0 = 1,
        Axis1      = CAxis1,
        Axis2      = CAxis2,
        Axis3      = CAxis3,
        Axis4      = CAxis4,
        samplingA1 = CsamplingA1
      ),
      list(i = 1:nv)
    ),
    tag = 'pp'
  )

  # --- 3.2 Initial fit to learn hyperparameters ---
  ppVIV <- inla(
    y ~ 0 + b0 + f(i, model = matern),
    family = 'poisson', data = inla.stack.data(stk.pp),
    control.predictor = list(A = inla.stack.A(stk.pp)),
    E = inla.stack.data(stk.pp)$e,
    control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
    control.inla     = list(strategy = 'simplified.laplace', int.strategy = 'eb'),
    num.threads = 20
  )

  sumsum   <- summary(ppVIV)
  rangeFIX <- sumsum$hyperpar$mean[1]
  sigmaFIX <- sumsum$hyperpar$mean[2]

  # --- 3.3 Refit with fixed range/sigma and covariates ---
  matern <- inla.spde2.pcmatern(
    mesh,
    alpha = 2,
    prior.sigma = c(sigmaFIX, NA),
    prior.range = c(rangeFIX, NA)
  )

  ppVIV <- inla(
    y ~ 0 + b0 + Axis1 + Axis2 + Axis3 + Axis4 + samplingA1 + f(i, model = matern),
    family = 'poisson', data = inla.stack.data(stk.pp),
    control.predictor = list(A = inla.stack.A(stk.pp)),
    E = inla.stack.data(stk.pp)$e,
    control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
    control.inla     = list(strategy = 'simplified.laplace', int.strategy = 'eb'),
    num.threads = 20
  )

  # --- 3.4 Build a fauna-intensity covariate raster from one posterior sample ---
  Nrep <- 100                # repetitions for held-out buffer predictions
  num_cores <- 15            # parallel workers for road stage

  pr.int.tot <- inla.posterior.sample(1, ppVIV)  # single sample of latent field
  predTMP <- inla.mesh.project(PROJgridVIV, pr.int.tot[[1]]$latent[1:mesh$n])

  # Rasterise projected intensity to extract values at mesh nodes/points later
  predTMP_sf <- st_sf(pxlVIV$pixels_sf)
  predTMP_sf$predTMP <- predTMP
  predTMP_R <- st_rasterize(predTMP_sf)
  predTMP_R <- as(predTMP_R, "Raster")
  predTMP_R <- scale(predTMP_R)  # standardise

  # -----------------------------
  # 4) Held-out buffer predictions for road mortality (parallel)
  # -----------------------------
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  results <- foreach(i = 1:Nrep, .packages = c("sf", "INLA", "raster", "stars", "INLAspacetime", "inlabru")) %dopar% {
    # ---- 4.1 Random buffer and data partition ----
    # One random point inside the domain; draw a 30-km (?) radius buffer
    random_point  <- st_sample(domaine, size = 1)
    buffer_radius <- 30
    buffer_area   <- st_buffer(random_point, dist = buffer_radius)

    # Remove collision points inside buffer; keep the rest for fitting
    datcoll_inside_buffer <- st_intersects(datcoll, buffer_area, sparse = FALSE)
    datcoll_remove    <- datcoll[datcoll_inside_buffer, ]     # not used directly here (but could be for inspection)
    datcoll_remaining <- datcoll[!datcoll_inside_buffer, ]

    # Exclude the buffer from the route-buffer domain when computing weights
    routeBUF_outside_buffer <- st_difference(routeBUF, buffer_area)
    intersect_idx <- st_intersects(dmeshROUTE, routeBUF_outside_buffer, sparse = FALSE)
    intersections <- st_intersection(dmeshROUTE[intersect_idx, ], routeBUF_outside_buffer)
    wROUTE <- numeric(nrow(dmeshROUTE))
    wROUTE[intersect_idx] <- as.numeric(st_area(intersections))

    # Coordinates for covariate extraction (mesh nodes + retained events)
    coordsMOR <- rbind(meshROUTE$loc[, 1:2], st_coordinates(datcoll_remaining))

    vitesse   <- extract(vitesseR,   coordsMOR)
    trafic    <- extract(traficR,    coordsMOR)
    disthydro <- extract(hydrodistR, coordsMOR)
    distveg   <- extract(haiedistR,  coordsMOR)
    Nbebete   <- extract(predTMP_R,  coordsMOR)

    # ---- 4.2 Barrier SPDE definition (road effects) ----
    maternROUTE <- barrierModel.define(
      mesh = meshROUTE,
      barrier.triangles = route.tri,
      prior.range = c(2, 0.01),   # Example PC prior: P(range < 2) = 0.01
      prior.sigma = c(2, 0.01),   # Example PC prior: P(sigma > 2) = 0.01
      range.fraction = 0.1
    )

    # Point-process stack on road mesh
    nv <- meshROUTE$n
    n  <- nrow(datcoll_remaining)
    y.pp <- rep(0:1, c(nv, n))
    e.pp <- c(wROUTE, rep(0, n))
    imat <- Diagonal(nv, rep(1, nv))
    lmat <- inla.spde.make.A(meshROUTE, st_coordinates(datcoll_remaining))
    A.pp <- rbind(imat, lmat)

    stk.pp <- inla.stack(
      data = list(y = y.pp, e = e.pp),
      A    = list(1, A.pp),
      effects = list(
        list(b0 = 1, Nbebete = Nbebete, vitesse = vitesse,
             trafic = trafic, disthydro = disthydro, distveg = distveg),
        list(i = 1:nv)
      ),
      tag = 'pp'
    )

    # ---- 4.3 Fit model & predict inside buffer ----
    ppMOR <- inla(
      y ~ 0 + b0 + Nbebete + vitesse + trafic + disthydro + distveg + f(i, model = maternROUTE),
      family = 'poisson', data = inla.stack.data(stk.pp),
      control.predictor = list(A = inla.stack.A(stk.pp)),
      E = inla.stack.data(stk.pp)$e,
      control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
      control.inla     = list(strategy = 'simplified.laplace', int.strategy = 'eb'),
      num.threads = 1
    )

    fitval <- inla.mesh.project(PROJgridMOR, ppMOR$summary.fitted.values$`0.5quant`[1:meshROUTE$n])
    F.effects <- ppMOR$summary.fixed
    rm("ppMOR"); gc()

    # Convert intensities to expected counts/AUC for each pixel
    predMORT_val <- fitval * pxlMOR$area_inter
    predAUC_val  <- 1 - exp(-as.numeric(fitval * pxlMOR$area_inter))

    # Package predictions limited to the held-out buffer for downstream assessment
    pred_remove <- st_sf(as.data.frame(cbind(predMORT_val = predMORT_val, predAUC_val = predAUC_val)),
                         geometry = st_geometry(pxlMOR$pixels_sf))
    pred_remove <- st_intersection(pred_remove, buffer_area)

    list(pred_remove = pred_remove, buffer_area = buffer_area, F.effects = F.effects)
  }

  # Stop parallel backend
  stopCluster(cl)

  # Persist results for this species group
  qs::qsave(results, file = paste0("PRED/PRED_", gsp, ".qs"), preset = "high")

  rm("results", "ppVIV", "pr.int.tot"); gc()
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
