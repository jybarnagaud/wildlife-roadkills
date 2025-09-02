# main_estimation_commented.R
# -----------------------------------------------------------------------------
# Title: Two-stage spatial point-process estimation (fauna presence & road mortality)
# Author: <Your name>
# Date: <YYYY-MM-DD>
#
# Description
# -----------
# Main workflow for estimating (i) the spatial intensity of live fauna presence
# over a study domain and (ii) road mortality intensity using a barrier SPDE
# model around roads. The analysis follows an LGCP-like Poisson point-process
# formulation using R-INLA with SPDE meshes, with a first stage to estimate
# spatial range/sigma (and covariate effects) for live fauna, then a second
# stage that uses those predictions (with the sampling effect removed) as a
# covariate in the road-mortality model.
#
# Provenance
# ----------
# This script adapts code patterns from the SPDE GitBook by V. Gómez-Rubio:
#   https://becarioprecario.bitbucket.io/spde-gitbook/index.html
# Please cite that resource and foundational references (Lindgren, Rue &
# Lindström, 2011) in your manuscript.
#
# Reproducibility & Computing
# ---------------------------
# - Record package versions via sessionInfo() in the manuscript supplement.
# - Paths are relative to the working directory set below.
# - Parallelisation: barrier and road-stage fits are parallelised using
#   doParallel/foreach. Adjust cores according to your machine.
# - INLA thread options are set explicitly in each call.
# -----------------------------------------------------------------------------

# -------------------------------
# 0) Setup
# -------------------------------
setwd("/mnt/c/Users/jpapaix/Documents/modelisation/mortaliteroutiere/codes_propres")

# Core packages
library(raster)       # legacy raster operations (used alongside terra)
library(ade4)         # PCA (if needed upstream)
library(sf)           # vector GIS
library(INLA)         # R-INLA core
library(terra)        # modern raster package
library(stars)        # rasterisation of sf
library(doParallel)   # parallel backend
library(foreach)      # foreach loops
library(qs)           # fast serialization of results
library(dplyr)        # data manipulation
library(ggplot2)      # plotting
inla.setOption(num.threads = 1)  # global default; overridden in inla() calls

# Helper functions (from your local utility files)
source("inlabookfunctions.R")  # helpers from SPDE book (e.g., barrier tools)
source("pixels.R")             # pixel grid helper (sf polygon grid)

# Additional spatial modelling packages
library(INLAspacetime)  # barrierModel.define, etc.
library(inlabru)        # mesh helpers & projectors

# -------------------------------
# 1) Covariates, meshes, and observation helpers
# -------------------------------
# These scripts are expected to:
# - build the study domain (sf), route buffer polygons, and meshes (mesh, meshROUTE)
# - construct dual meshes (dmesh, dmeshROUTE) and projector grids (pxlVIV, pxlMOR)
# - load or compute environmental rasters (Axis1R..Axis4R, vitesseR, traficR,
#   hydrodistR, haiedistR) and smoothed sampling rasters (samplingRA1, etc.)
# - define `listeSP` (mapping group -> species)
source("covariables_mesh.R")

# The observations helper is expected to provide `chargeDON()` which returns
# sf layers of live observations and collisions, and per-pixel counts
# (NobsVIV, NobsMOR) on the two projector grids.
source("observations.R")

# -------------------------------
# 2) Species-group loop
# -------------------------------
for (gsp in unique(listeSP[, "grpSP"])) {  # e.g., "Renard roux"
  # Species list for the current group
  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]

  # Load & filter observations for 2015–2019 (exclusive bounds in chargeDON)
  DON     <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  datviv  <- DON[[1]]  # sf POINTS: live observations
  datcoll <- DON[[2]]  # sf POINTS: collisions (roadkill)
  NobsVIV <- DON[[3]]  # counts per pxlVIV pixel
  NobsMOR <- DON[[4]]  # counts per pxlMOR pixel

  # Quick diagnostic map of data and meshes
  png(file = paste0("RES/don_", gsp, ".png"), width = 300, height = 200, res = 300, units = "mm")
  par(mar = rep(1, 4), mfrow = c(1, 1))
  plot(mesh, main = "", asp = 1)
  plot(routeBUF, add = TRUE, col = 2, border = 2, lwd = 2)
  plot(st_geometry(domaine), add = TRUE, border = 1, lwd = 2)
  plot(st_geometry(datviv),  add = TRUE, col = 3, pch = 16, cex = 0.8)
  plot(st_geometry(datcoll), add = TRUE, col = 4, pch = 16, cex = 0.8)
  dev.off()

  # -----------------------------
  # 3) Live fauna (LGCP-style) model
  # -----------------------------
  # Extract covariates at both mesh nodes and presence points. We do this by
  # stacking mesh node locations with point coordinates, then splitting in A-matrix.
  coords <- rbind(mesh$loc[, 1:2], st_coordinates(datviv))

  CAxis1      <- extract(Axis1R,      coords)
  CAxis2      <- extract(Axis2R,      coords)
  CAxis3      <- extract(Axis3R,      coords)
  CAxis4      <- extract(Axis4R,      coords)
  CsamplingA1 <- extract(samplingRA1,  coords)
  CsamplingA2 <- extract(samplingRA2,  coords)  # available if needed

  # --- 3.1 Prior and precision for SPDE field (PC priors) ---
  matern <- inla.spde2.pcmatern(
    mesh,
    alpha = 2,                       # SPDE order (nu = 1)
    prior.sigma = c(1.6, 0.1),       # P(sigma > 1.6) = 0.1 (interpret per your choice)
    prior.range = c(55, 0.9)         # P(range < 55) = 0.9
  )

  # Integration weights over dual mesh cells (intersection with domain)
  intersect_idx <- st_intersects(dmesh, domaine, sparse = FALSE)
  intersections <- st_intersection(dmesh[intersect_idx, ], domaine)
  w <- numeric(nrow(dmesh))
  w[intersect_idx] <- as.numeric(st_area(intersections))

  # Point-process stack: Append mesh integration points (offset e = w) with events
  nv  <- mesh$n
  n   <- nrow(datviv)
  y.pp <- rep(0:1, c(nv, n))
  e.pp <- c(w, rep(0, n))
  imat <- Diagonal(nv, rep(1, nv))
  lmat <- inla.spde.make.A(mesh, st_coordinates(datviv))
  A.pp <- rbind(imat, lmat)

  # Split covariate vectors for mesh nodes vs. points (first nv entries)
  stk.pp <- inla.stack(
    data = list(y = y.pp, e = e.pp),
    A    = list(1, A.pp),
    effects = list(
      list(
        b0 = 1,
        Axis1      = CAxis1,      # length = nv + n
        Axis2      = CAxis2,
        Axis3      = CAxis3,
        Axis4      = CAxis4,
        samplingA1 = CsamplingA1  # include sampling bias proxy in stage-1
        # samplingA2 = CsamplingA2 # available if used
      ),
      list(i = 1:nv)
    ),
    tag = 'pp'
  )

  # --- 3.2 Initial fit to learn range/sigma ---
  ppVIV <- inla(
    y ~ 0 + b0 + f(i, model = matern),
    family = 'poisson',
    data = inla.stack.data(stk.pp),
    control.predictor = list(A = inla.stack.A(stk.pp)),
    E = inla.stack.data(stk.pp)$e,
    control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
    control.inla     = list(strategy = 'simplified.laplace', int.strategy = 'eb'),
    num.threads = 20
  )

  # Extract posterior means for hyperparameters (range, sigma)
  sumsum   <- summary(ppVIV)
  rangeFIX <- sumsum$hyperpar$mean[1]
  sigmaFIX <- sumsum$hyperpar$mean[2]

  # --- 3.3 Refit with fixed hyperparameters and covariates ---
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

  # Project posterior fitted intensity (mesh nodes only) to prediction grid
  fitvalmin <- inla.mesh.project(PROJgridVIV, ppVIV$summary.fitted.values$`0.025quant`[1:mesh$n])
  fitvalmed <- inla.mesh.project(PROJgridVIV, ppVIV$summary.fitted.values$`0.5quant`  [1:mesh$n])
  fitvalmax <- inla.mesh.project(PROJgridVIV, ppVIV$summary.fitted.values$`0.975quant`[1:mesh$n])
  predVIV   <- rbind(fitvalmin, fitvalmed, fitvalmax) * pxlVIV$area_pix
  predVIVAUC <- 1 - exp(-as.numeric(fitvalmed * pxlVIV$area_pix))

  save(predVIVAUC, ppVIV, predVIV, NobsVIV, NobsMOR, file = paste0("RES/fitVIV_obs_", gsp, ".rda"))

  # -----------------------------
  # 4) Road mortality model (barrier SPDE)
  # -----------------------------
  # Compute integration weights for the road-buffer dual mesh
  intersect_idx <- st_intersects(dmeshROUTE, routeBUF, sparse = FALSE)
  intersections <- st_intersection(dmeshROUTE[intersect_idx, ], routeBUF)
  wROUTE <- numeric(nrow(dmeshROUTE))
  wROUTE[intersect_idx] <- as.numeric(st_area(intersections))

  # Number of posterior samples carried from stage-1 into stage-2
  Nrep <- 100
  num_cores <- 15

  # Draw posterior samples of the latent field from the live-fauna model
  pr.int.tot <- inla.posterior.sample(Nrep, ppVIV)

  # Combine mesh nodes and collision point coordinates for covariate extraction
  coordsMOR <- rbind(meshROUTE$loc[, 1:2], st_coordinates(datcoll))

  vitesse   <- extract(vitesseR,   coordsMOR)
  trafic    <- extract(traficR,    coordsMOR)
  disthydro <- extract(hydrodistR, coordsMOR)
  distveg   <- extract(haiedistR,  coordsMOR)

  # Parallel backend for repeated barrier fits
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  results <- foreach(i = 1:Nrep, .packages = c("sf", "INLA", "raster", "stars", "INLAspacetime", "inlabru")) %dopar% {
    # ---- 4.1 Build fauna-intensity covariate without sampling bias ----
    eta_with  <- pr.int.tot[[i]]$latent[1:mesh$n]
    beta_samplingA1 <- pr.int.tot[[i]]$latent[grep("^samplingA1", rownames(pr.int.tot[[i]]$latent))]
    samplingA1_covariate_values <- CsamplingA1[1:mesh$n]
    # Remove the contribution of samplingA1 from the latent field
    eta_without <- eta_with - beta_samplingA1 * samplingA1_covariate_values

    predTMP <- inla.mesh.project(PROJgridVIV, eta_without)

    # Rasterise projected intensity to later extract as covariate at coordsMOR
    predTMP_sf <- st_sf(pxlVIV$pixels_sf)
    predTMP_sf$predTMP <- predTMP
    predTMP_R <- st_rasterize(predTMP_sf)
    predTMP_R <- as(predTMP_R, "Raster")
    predTMP_R <- scale(predTMP_R)  # standardise for numerical stability
    Nbebete   <- extract(predTMP_R, coordsMOR)

    # ---- 4.2 Define barrier SPDE model over road-buffer mesh ----
    maternROUTE <- barrierModel.define(
      mesh = meshROUTE,
      barrier.triangles = route.tri,
      prior.range = c(2, 0.01),  # Example PC prior: P(range < 2) = 0.01
      prior.sigma = c(2, 0.01),  # Example PC prior: P(sigma > 2) = 0.01
      range.fraction = 0.1
    )

    # Poisson point-process stack on the road mesh
    nv <- meshROUTE$n
    n  <- nrow(datcoll)
    y.pp <- rep(0:1, c(nv, n))
    e.pp <- c(wROUTE, rep(0, n))
    imat <- Diagonal(nv, rep(1, nv))
    lmat <- inla.spde.make.A(meshROUTE, st_coordinates(datcoll))
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

    # ---- 4.3 Fit INLA barrier model ----
    ppMOR <- inla(
      y ~ 0 + b0 + Nbebete + vitesse + trafic + disthydro + distveg + f(i, model = maternROUTE),
      family = 'poisson', data = inla.stack.data(stk.pp),
      control.predictor = list(A = inla.stack.A(stk.pp)),
      E = inla.stack.data(stk.pp)$e,
      control.compute = list(dic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
      control.inla     = list(strategy = 'simplified.laplace', int.strategy = 'eb'),
      num.threads = 1
    )

    # Median fitted intensity projected to road prediction grid
    fitval <- inla.mesh.project(PROJgridMOR, ppMOR$summary.fitted.values$`0.5quant`[1:meshROUTE$n])
    F.effects <- ppMOR$summary.fixed

    # Posterior sign probabilities for fixed effects (based on 500 samples)
    post <- inla.posterior.sample(500, ppMOR)
    postsamp <- lapply(post, function(x) x$latent[c("Nbebete:1", "vitesse:1", "trafic:1", "disthydro:1", "distveg:1"), ])
    postsamp <- matrix(unlist(postsamp), ncol = 5, byrow = TRUE)
    proba <- apply(postsamp > 0, 2, sum) / 500

    rm("ppMOR", "post"); gc()

    # Convert intensities to expected counts / AUC per pixel
    predMORT_val <- fitval * pxlMOR$area_inter
    predAUC_val  <- 1 - exp(-as.numeric(fitval * pxlMOR$area_inter))

    list(predMORT = predMORT_val, predAUC = predAUC_val, F.effects = F.effects, probpost = proba)
  }

  stopCluster(cl)

  # Persist outputs for the current species group
  qs::qsave(results, file = paste0("RES/RES_withoutsampling_", gsp, ".qs"), preset = "high")

  rm("results", "ppVIV", "pr.int.tot"); gc()
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
