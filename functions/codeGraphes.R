# =============================================================================
# Title:     Spatio-temporal wildlife–road mortality analysis and visualization
# Purpose:   Prepare inputs, compute model diagnostics/metrics, and produce
#            publication-quality figures (mesh+points, fit & CV histograms,
#            observed vs. predicted plots, and forest/violin-style summaries).
# Authors:   <Add full names + ORCIDs>
# Affiliation: <Add lab/institute>
# Contact:   <corresponding.author@institution>
# License:   <MIT/CC-BY/… per journal policy>
# Reproducibility:
#   - INLA threads fixed to 1 (deterministic Laplace approximations).
#   - Set a random seed if your sourced scripts simulate (e.g., set.seed(123)).
#   - Record session info at the end (sessionInfo()).
# Data:
#   - Assumes project tree with `data/`, `RES/`, `PRED/`, and helper scripts.
#   - Observations are loaded through chargeDON() (from observations.R).
# Figures:
#   - PNG outputs written to the working directory with 300 dpi (journal-grade).
# Notes:
#   - This file is a script (not a package). For a code supplement, consider
#     wrapping utilities as functions and documenting them with roxygen2.
# =============================================================================


# -----------------------------------------------------------------------------
# 0) Working directory
#    In production, prefer project-rooted paths via here::here(). Kept as-is
#    to match your current environment. Consider making this configurable.
# -----------------------------------------------------------------------------
setwd("/mnt/c/Users/jpapaix/Documents/modelisation/mortaliteroutiere/codes_propres")
# Alternative (recommended):
# if (requireNamespace("here", quietly = TRUE)) setwd(here::here())


# -----------------------------------------------------------------------------
# 1) Libraries
#    Spatial: sf, terra, stars, raster (legacy).
#    Modeling: INLA, inlabru, INLAspacetime.
#    Plotting: ggplot2, cowplot, patchwork, gridExtra/grid.
#    Utilities: dplyr, tidyr, purrr, tibble, forcats, qs.
#    Diagnostics: pROC.
# -----------------------------------------------------------------------------
library(raster)        # legacy raster support; be mindful of classes vs. terra
library(ade4)          # (kept—used by sourced scripts)
library(sf)            # simple features (vector)
library(INLA)          # Bayesian inference via INLA
library(terra)         # modern raster/vector (successor to raster)
library(stars)         # spatiotemporal arrays
library(doParallel)    # foreach backend (kept if used downstream)
library(foreach)       # foreach (kept if used in sourced code)
library(qs)            # fast serialization (.qs)
library(pROC)          # ROC/AUC
library(ggplot2)       # plotting
library(patchwork)     # compose ggplots
library(dplyr)         # data manipulation
library(tidyr)         # tidy reshaping
library(INLAspacetime) # helpers for INLA spatiotemporal models
library(inlabru)       # spatial components for INLA
library(gridExtra)     # arrange grobs
library(grid)          # low-level grid graphics
library(forcats)       # factor helpers
library(tibble)        # tidy tibbles
library(ggpubr)        # publication helpers for ggplot
library(cowplot)       # themes and plot grids
library(purrr)         # functionals
library(colorspace)    # color manipulation


# -----------------------------------------------------------------------------
# 2) INLA computational settings
#    Deterministic behavior: 1 thread. Consider also setting set.seed().
# -----------------------------------------------------------------------------
inla.setOption(num.threads = 1)
# set.seed(123)  # uncomment if any stochasticity appears in sourced code


# -----------------------------------------------------------------------------
# 3) Source project-specific helper scripts
#    - inlabookfunctions.R : custom INLA/inlabru helpers
#    - pixels.R            : grid/pixel utilities for rasterization/meshes
# -----------------------------------------------------------------------------
source("inlabookfunctions.R")
source("pixels.R")


# -----------------------------------------------------------------------------
# 4) Load covariates and mesh definition
#    Expected side effects (objects in env):
#      - mesh      : SPDE mesh
#      - domaine   : sf polygon (study domain)
#      - routeBUF  : roads or buffered roads for context plotting
#      - covariates as required elsewhere (e.g., pxlMOR in later sections)
# -----------------------------------------------------------------------------
source("covariables_mesh.R")


# -----------------------------------------------------------------------------
# 5) Load or construct observation data (species-list & loader)
#    Must define:
#      - listeSP : data.frame with columns incl. "grpSP" and species name
#      - chargeDON(espece, annee_min, annee_max, domaine, file=...):
#          returns list:
#             [[1]] live (sf POINT), [[2]] dead (sf POINT),
#             [[3]] NobsVIV (int),   [[4]] NobsMOR (int)
# -----------------------------------------------------------------------------
source("observations.R")


# -----------------------------------------------------------------------------
# 6) Count observations by species group (years 2015–2020)
#    Rationale: transparent summary of live vs mortality sample sizes used in
#    modeling/diagnostics, by species group.
# -----------------------------------------------------------------------------
Nmort   <- NULL  # mortality counts by group
Nvivant <- NULL  # live-sighting counts by group

for (gsp in unique(listeSP[, "grpSP"])) {
  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]

  # Returns live, dead, counts; limited to [2015, 2020] and to 'domaine'
  DON <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")

  # Keep raw row counts (rather than reported integers) for full transparency
  Nvivant <- c(Nvivant, nrow(DON[[1]]))
  Nmort   <- c(Nmort,   nrow(DON[[2]]))
}

Ndon <- data.frame(Nvivant = Nvivant, Nmort = Nmort)
rownames(Ndon) <- unique(listeSP[, "grpSP"])
# Optional QA:
# print(Ndon); colSums(Ndon)


# -----------------------------------------------------------------------------
# 7) Figure: mesh + raw points (one PNG per species group)
#    Plots:
#      - mesh edges
#      - road buffer/network context
#      - study domain boundary
#      - live (purple) vs dead (orange) points
#    Colors: Okabe–Ito (colorblind-safe).
#    Output: don_<grpSP>.png (300mm × 200mm @ 300 dpi).
# -----------------------------------------------------------------------------

# Okabe–Ito palette (colorblind-friendly)
cols <- c(
  black      = "#000000",
  orange     = "#E69F00",
  skyblue    = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  purple     = "#CC79A7"
)

for (gsp in listeSP[, "grpSP"]) {

  espece  <- listeSP[listeSP[, "grpSP"] == gsp, 2]
  DON     <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  datviv  <- DON[[1]]  # sf points (live)
  datcoll <- DON[[2]]  # sf points (dead)
  NobsVIV <- DON[[3]]
  NobsMOR <- DON[[4]]

  # Device: journal-friendly size & resolution
  png(file = paste0("don_", gsp, ".png"),
      width = 300, height = 200, res = 300, units = "mm")

  par(mar = rep(1, 4), mfrow = c(1, 1))

  # Mesh (ensure CRS matches data/routeBUF; transform if needed)
  plot(mesh, main = "", asp = 1, edge.color = grey(0.8))

  # Roads/buffers layer
  plot(routeBUF, add = TRUE, col = cols["skyblue"], border = cols["skyblue"], lwd = 2)

  # Study domain
  plot(st_geometry(domaine), add = TRUE, border = cols["black"], lwd = 2)

  # Points: live vs dead
  plot(st_geometry(datviv),  add = TRUE, col = cols["purple"], pch = 16, cex = 0.7)
  plot(st_geometry(datcoll), add = TRUE, col = cols["orange"], pch = 16, cex = 0.7)

  # Optional legend (commented for cleaner figure)
  # legend("topleft",
  #        legend = c(paste0("Live (n=", NobsVIV, ")"),
  #                   paste0("Dead (n=", NobsMOR, ")")),
  #        pch = 16, col = c(cols["purple"], cols["orange"]),
  #        pt.cex = 0.8, bty = "n")

  dev.off()
}


# =============================================================================
# MODEL METRICS & EVALUATION PLOTS — DEAD ANIMALS
# =============================================================================

# -----------------------------------------------------------------------------
# 8) Model fit metrics on dead animals (AUC, TSS, R², COR; also store preds)
#    Input: RES/RES_withoutsampling_<grpSP>.qs (list of length Nrep=100)
#           Each element contains:
#             - $predAUC: numeric vector of probabilities for ROC/AUC
#             - $predMORT: numeric vector of predicted counts (fit)
#    Output: metricsFIT.rda with AUC, TSS, PRED, R2, COR, NobsMOR_res0
# -----------------------------------------------------------------------------
AUCtot <- list(); TSStot <- list(); PREDtot <- list()
R2tot  <- list(); CORtot <- list(); NobsMOR_res0 <- list()

i <- 0
for (gsp in unique(listeSP[, "grpSP"])) {
  i <- i + 1
  results <- qread(paste0("RES/RES_withoutsampling_", gsp, ".qs"))

  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]
  DON    <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  NobsMOR <- DON[[4]]
  NobsMOR_bis <- as.numeric(NobsMOR > 0)  # presence/absence for AUC/TSS

  Nrep <- 100
  AUC <- TSS <- R2 <- COR <- NULL
  PRED <- NULL

  for (j in 1:Nrep) {
    predAUC <- results[[j]]$predAUC
    # ROC; handle cases where predAUC may be constant by pROC’s internal checks
    testroc <- roc(NobsMOR_bis, predAUC)
    AUC <- c(AUC, testroc$auc)
    TSS <- c(TSS, max(testroc$sensitivities + testroc$specificities - 1))

    # Predicted vs observed (fit on training data)
    PRED <- rbind(PRED, results[[j]]$predMORT)
    mod  <- lm(results[[j]]$predMORT ~ -1 + NobsMOR)  # slope-only fit
    R2   <- c(R2, summary(mod)$adj.r.squared)
    COR  <- c(COR, cor(results[[j]]$predMORT, NobsMOR))
  }

  AUCtot[[i]] <- AUC
  TSStot[[i]] <- TSS
  PREDtot[[i]] <- PRED
  R2tot[[i]] <- R2
  CORtot[[i]] <- COR
  NobsMOR_res0[[i]] <- NobsMOR
}

names(AUCtot) <- names(TSStot) <- names(PREDtot) <- names(R2tot) <-
  names(CORtot) <- names(NobsMOR_res0) <- unique(listeSP[, "grpSP"])

save(NobsMOR_res0, AUCtot, TSStot, PREDtot, R2tot, CORtot,
     file = "metricsFIT.rda")


# -----------------------------------------------------------------------------
# 9) Cross-validated prediction metrics on dead animals
#    Logic: remove-buffer evaluation, compute RMSE (relative), AUC, TSS, COR.
#    Inputs: PRED/PRED_<grpSP>.qs (list length 100) with elements:
#            - $buffer_area: sf polygon of removal buffer
#            - $pred_remove$predMORT_val: numeric predictions on held-out area
#            - $pred_remove$predAUC_val: probabilities for ROC
#    Output: metricsPRED.rda with AUC, TSS, totopred, err_pred, COR, totoobs
# -----------------------------------------------------------------------------
AUCtot <- list(); TSStot <- list(); totopredtot <- list()
err_predtot <- list(); CORtot <- list(); totoobs_tot <- list()

i <- 0
for (gsp in unique(listeSP[, "grpSP"])) {
  i <- i + 1

  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]
  DON    <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  Nobs_MORT <- DON[[4]]

  # Tie observations to a pixel grid (assumes pxlMOR$pixels_sf exists)
  obs_MORT <- st_sf(as.data.frame(cbind(Nobs_MORT = Nobs_MORT)),
                    geometry = st_geometry(pxlMOR$pixels_sf))

  results <- qread(paste0("PRED/PRED_", gsp, ".qs"))

  err_pred <- NULL; totoobs <- NULL; totopred <- NULL
  AUC <- NULL; TSS <- NULL; COR <- NULL

  for (j in 1:100) {
    # Intersect held-out area and compute errors/metrics
    obs_remove <- st_intersection(obs_MORT, results[[j]]$buffer_area)

    # Relative RMSE: sqrt(MSE)/mean(obs)
    predTMP <- sqrt(mean((obs_remove$Nobs_MORT -
                          results[[j]]$pred_remove$predMORT_val)^2)) /
               mean(obs_remove$Nobs_MORT)

    err_pred <- c(err_pred, predTMP)
    totoobs  <- rbind(totoobs, cbind(as.numeric(rownames(obs_remove)),
                                     obs_remove$Nobs_MORT))
    totopred <- c(totopred, results[[j]]$pred_remove$predMORT_val)

    predAUC      <- results[[j]]$pred_remove$predAUC_val
    NobsMOR_bis  <- as.numeric(obs_remove$Nobs_MORT > 0)

    if (length(unique(NobsMOR_bis)) > 1) {
      testroc <- roc(NobsMOR_bis, predAUC)
      AUC <- c(AUC, testroc$auc)
      TSS <- c(TSS, max(testroc$sensitivities + testroc$specificities - 1))
    } else {
      AUC <- c(AUC, NA); TSS <- c(TSS, NA)
    }

    COR <- c(COR, cor(results[[j]]$pred_remove$predMORT_val, obs_remove$Nobs_MORT))
  }

  AUCtot[[i]]       <- AUC
  TSStot[[i]]       <- TSS
  totopredtot[[i]]  <- totopred
  err_predtot[[i]]  <- err_pred
  CORtot[[i]]       <- COR
  totoobs_tot[[i]]  <- totoobs
}

names(AUCtot) <- names(TSStot) <- names(totopredtot) <-
  names(err_predtot) <- names(CORtot) <- names(totoobs_tot) <-
  unique(listeSP[, "grpSP"])

save(totoobs_tot, AUCtot, TSStot, totopredtot, err_predtot, CORtot,
     file = "metricsPRED.rda")


# -----------------------------------------------------------------------------
# 10) FIT & CV plots (DEAD animals): histograms (AUC/TSS) + Obs vs Pred
#     Output: eval_<grpSP>.png (3×3 grid with titles and facets)
# -----------------------------------------------------------------------------
for (gsp in unique(listeSP[, "grpSP"])) {

  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]
  DON    <- chargeDON(espece, 2015, 2020, domaine, file = "data/donBON.rda")
  NobsMOR <- DON[[4]]

  # --- Fit metrics & predicted summaries (within-sample) ----------------------
  load("metricsFIT.rda")
  AUC <- AUCtot[[gsp]]
  TSS <- TSStot[[gsp]]
  R2  <- R2tot[[gsp]]
  COR <- CORtot[[gsp]]
  PRED <- PREDtot[[gsp]]  # matrix [Nrep × Ncells]

  PRED.stat <- apply(PRED, 2, quantile, probs = c(0.025, 0.5, 0.975))

  df_pred <- data.frame(
    Observed    = NobsMOR,
    Pred_median = PRED.stat[2, ],
    Pred_lower  = PRED.stat[1, ],
    Pred_upper  = PRED.stat[3, ]
  )

  fit <- lm(Pred_median ~ Observed, data = df_pred)

  p1FIT <- ggplot(data.frame(AUC = AUC), aes(x = AUC)) +
    geom_histogram(binwidth = 0.0005, fill = "gray80", color = "black") +
    theme_classic() + ggtitle("(a)") +
    xlab("AUC") + ylab("Counts")

  p2FIT <- ggplot(data.frame(TSS = TSS), aes(x = TSS)) +
    geom_histogram(binwidth = 0.004, fill = "gray80", color = "black") +
    theme_classic() + ggtitle("(b)") +
    xlab("TSS") + ylab("Counts")

  p3FIT <- ggplot(df_pred, aes(x = Observed, y = Pred_median)) +
    geom_point(size = 2, color = "black") +
    geom_errorbar(aes(ymin = Pred_lower, ymax = Pred_upper), width = 0, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "blue", linewidth = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
                color = "goldenrod1", linewidth = 1) +
    theme_classic() + ggtitle("(c)") +
    xlab("Observed") + ylab("Predicted")

  # --- Cross-validated metrics -----------------------------------------------
  load("metricsPRED.rda")
  AUC     <- AUCtot[[gsp]]
  TSS     <- TSStot[[gsp]]
  totoobs <- totoobs_tot[[gsp]]
  totopred <- totopredtot[[gsp]]

  p1PRED <- ggplot(data.frame(AUC = AUC), aes(x = AUC)) +
    geom_histogram(binwidth = 0.0218, fill = "gray80", color = "black") +
    theme_classic() + ggtitle("(d)") +
    xlab("AUC") + ylab("Counts")

  p2PRED <- ggplot(data.frame(TSS = TSS), aes(x = TSS)) +
    geom_histogram(binwidth = 0.062, fill = "gray80", color = "black") +
    theme_classic() + ggtitle("(e)") +
    xlab("TSS") + ylab("Counts")

  PRED <- cbind(totoobs, totopred) |> as.data.frame()
  colnames(PRED) <- c("ID", "OBS", "PRED")

  summary_pred <- PRED |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      OBS       = unique(OBS),
      PRED_2.5  = quantile(PRED, probs = 0.025),
      PRED_50   = quantile(PRED, probs = 0.5),
      PRED_97.5 = quantile(PRED, probs = 0.975)
    ) |>
    as.data.frame()

  p3PRED <- ggplot(summary_pred, aes(x = OBS, y = PRED_50)) +
    geom_point(size = 2, color = "black") +
    geom_errorbar(aes(ymin = PRED_2.5, ymax = PRED_97.5), width = 0, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "blue", linewidth = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
                color = "goldenrod1", linewidth = 1) +
    theme_classic() + ggtitle("(f)") +
    xlab("Observed") + ylab("Predicted")

  # --- Arrange as 3×3 grid with row/column labels ----------------------------
  col1_title <- textGrob("AUC", gp = gpar(fontface = "bold", fontsize = 12))
  col2_title <- textGrob("TSS", gp = gpar(fontface = "bold", fontsize = 12))
  col3_title <- textGrob("Observed vs. Predicted", gp = gpar(fontface = "bold", fontsize = 12))

  row1_title <- textGrob("Model fit", rot = 90, gp = gpar(fontface = "bold", fontsize = 12))
  row2_title <- textGrob("Cross-validation", rot = 90, gp = gpar(fontface = "bold", fontsize = 12))

  png(file = paste0("eval_", gsp, ".png"), width = 300, height = 200, res = 300, units = "mm")
  grid.arrange(
    nullGrob(), col1_title, col2_title, col3_title,
    row1_title, p1FIT, p2FIT, p3FIT,
    row2_title, p1PRED, p2PRED, p3PRED,
    ncol = 4, nrow = 3,
    widths = c(0.5, 4, 4, 4),
    heights = c(0.5, 4, 4)
  )
  dev.off()
}


# =============================================================================
# LIVE ANIMALS — FIT DIAGNOSTICS & FOREST PLOT
# =============================================================================

# -----------------------------------------------------------------------------
# 11) Live animals: AUC/TSS and Observed vs Predicted (per group)
#     Input: RES/fitVIV_obs_<grpSP>.rda with objects:
#            - NobsVIV (numeric), predVIV (3×N matrix, [2,] = median),
#              predVIVAUC (probabilities for ROC), ppVIV (INLA result)
#     Output: evalVIV_<grpSP>.png
# -----------------------------------------------------------------------------
for (gsp in unique(listeSP[, "grpSP"])) {
  espece <- listeSP[listeSP[, "grpSP"] == gsp, 2]
  load(paste0("RES/fitVIV_obs_", gsp, ".rda"))

  NobsVIV01 <- as.numeric(NobsVIV > 0)
  testroc   <- roc(NobsVIV01, predVIVAUC)
  AUCviv    <- testroc$auc
  TSSviv    <- max(testroc$sensitivities + testroc$specificities - 1)

  df <- data.frame(
    Nobs  = NobsVIV,
    pred  = as.numeric(predVIV[2, ]),
    lower = as.numeric(predVIV[1, ]),
    upper = as.numeric(predVIV[3, ])
  )

  mod <- lm(pred ~ Nobs, data = df)  # regression (golden line in plot)

  cols <- c(
    black      = "#000000", orange     = "#E69F00", skyblue    = "#56B4E9",
    green      = "#009E73", yellow     = "#F0E442", blue       = "#0072B2",
    vermillion = "#D55E00", purple     = "#CC79A7", lightgrey  = "#B0B0B0"
  )

  pp <- ggplot(df, aes(x = Nobs, y = pred)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0,
                  linewidth = 0.5, color = cols["lightgrey"]) +
    geom_point(size = 2.2, alpha = 0.95, color = cols["green"]) +
    geom_abline(slope = 1, intercept = 0, linetype = "longdash",
                linewidth = 0.8, color = cols["blue"]) +
    geom_abline(slope = coef(mod)[2], intercept = coef(mod)[1],
                linewidth = 1.1, color = cols["vermillion"]) +
    labs(
      x = "Observed", y = "Predicted",
      title = paste0("AUC = ", round(AUCviv, 3),
                     "    TSS = ", round(TSSviv, 3)),
      subtitle = sprintf("Regression: y = %.3f + %.3f × x",
                         unname(coef(mod)[1]), unname(coef(mod)[2]))
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  png(file = paste0("evalVIV_", gsp, ".png"), width = 150, height = 200, res = 300, units = "mm")
  print(pp)
  dev.off()
}


# -----------------------------------------------------------------------------
# 12) Live animals: forest plot (posterior sign probabilities & 95% CIs)
#     Method:
#       - Draw 500 posterior samples (inla.posterior.sample) from ppVIV.
#       - Extract latent effects (b0, Axis1–Axis4, samplingA1).
#       - Flip signs for Axis3/4/samplingA1 to align interpretation (as in orig).
#       - Compute posterior P(effect > 0) and 2.5/50/97.5% summaries.
#     Output: effectVIV.png
# -----------------------------------------------------------------------------
resVIV <- list(); i <- 0
for (gsp in unique(listeSP[, "grpSP"])) {
  i <- i + 1
  load(paste0("RES/fitVIV_obs_", gsp, ".rda"))

  post <- inla.posterior.sample(500, ppVIV)
  postsamp <- lapply(post, function(x) {
    x$latent[c("b0:1", "Axis1:1", "Axis2:1", "Axis3:1", "Axis4:1", "samplingA1:1"), ]
  })
  postsamp <- matrix(unlist(postsamp), ncol = 6, byrow = TRUE)

  # Sign flips (document rationale in manuscript)
  postsamp[, c(3, 4, 6)] <- -postsamp[, c(3, 4, 6)]

  proba <- apply(postsamp > 0, 2, sum) / 500
  toto  <- summary(ppVIV)$fixed
  toto[c(3, 4, 6), c(1, 3, 4, 5, 6)] <- -toto[c(3, 4, 6), c(1, 5, 4, 3, 6)]
  toto <- cbind(toto, proba)

  rownames(toto) <- c("Intercept","Agricultural landscapes vs. forests",
                      "Rural vs. urban","Density of watercourses",
                      "Bocages vs. open fields","Sampling effort")
  resVIV[[i]] <- toto
}

names(resVIV) <- c("Wild boar","Red fox","Other mustelids","European badger",
                   "Roe deer","European hare","Red deer")

# Tidy into a single df with probability bins for coloring
df_all <- bind_rows(lapply(names(resVIV), function(species) {
  df <- as.data.frame(resVIV[[species]])
  df <- tibble::rownames_to_column(df, "Parameter") |>
    dplyr::filter(Parameter %in% c("Agricultural landscapes vs. forests",
                                   "Rural vs. urban",
                                   "Density of watercourses",
                                   "Bocages vs. open fields",
                                   "Sampling effort")) |>
    dplyr::mutate(Species = species)
}))

df_all$Species <- factor(
  df_all$Species,
  levels = c("Wild boar","Red fox","Other mustelids","European hare","Roe deer","Red deer","European badger")
)

df_all$Parameter <- factor(
  df_all$Parameter,
  levels = c("Agricultural landscapes vs. forests","Rural vs. urban","Density of watercourses","Bocages vs. open fields","Sampling effort"),
  labels = c("(a) Agricultural landscapes vs. forests","(b) Rural vs. urban","(c) Density of watercourses","(d) Bocages vs. open fields","(e) Sampling effort")
)

breaks <- c(0, 0.05, 0.1, 0.9, 0.95, 1)
labels <- c("[0–0.05)", "[0.05–0.1)", "[0.1–0.9)", "[0.9–0.95)", "[0.95–1]")
df_all$pp_bin <- cut(df_all$proba, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)

bin_colors <- c(
  "[0–0.05)"   = "#E41A1C",
  "[0.05–0.1)" = darken("firebrick", amount = 0.3),
  "[0.1–0.9)"  = "black",
  "[0.9–0.95)" = darken("steelblue", amount = 0.3),
  "[0.95–1]"   = "#1E90FF"
)

pp <- ggplot(df_all, aes(x = mean, y = Species, color = pp_bin)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`), height = 0) +
  geom_point(size = 2) +
  scale_color_manual(values = bin_colors,
                     name = "Posterior Probability\n(effect > 0)") +
  facet_wrap(~ Parameter, scales = "free_x", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10))) +
  labs(x = "Estimate — 95% CI", y = "Species")

legend <- get_legend(pp + theme(legend.position = "right"))
pp_nolegend <- pp + theme(legend.position = "none")
final_plot <- pp_nolegend + inset_element(legend, left = 0.7, bottom = 0.05, right = 0.8, top = 0.2)

png(file = "effectVIV.png", width = 300, height = 200, units = "mm", res = 200)
final_plot
dev.off()


# =============================================================================
# DEAD ANIMALS — FOREST/“BAND” PLOTS OF FIXED EFFECTS
# =============================================================================

# -----------------------------------------------------------------------------
# 13) Aggregate fixed effects across CV replicates for each species & parameter
#     Inputs: RES/RES_withoutsampling_<grpSP>.qs (Nrep = 100)
#             F.effects columns include 0.025quant, 0.5quant, 0.975quant
#             probpost = posterior sign probability per replicate
#     Output: effectMORT.png + estim_effects.rda for reuse
# -----------------------------------------------------------------------------
Efaune <- list(); Evitesse <- list(); Etrafic <- list()
Edisthydro <- list(); Edistveg <- list()

il <- 0; Nrep <- 100
for (gsp in unique(listeSP[, "grpSP"])) {
  il <- il + 1
  results <- qread(paste0("RES/RES_withoutsampling_", gsp, ".qs"))
  espece  <- listeSP[listeSP[, "grpSP"] == gsp, 2]

  Efaunetmp <- Evitessetmp <- Etrafictmp <- Edisthydrotmp <- Edistvegtmp <- NULL

  for (j in 1:Nrep) {
    if (nrow(results[[j]]$F.effects) > 0) {
      Efaunetmp     <- rbind(Efaunetmp,     c(unlist(results[[j]]$F.effects["Nbebete",   3:5]), results[[j]]$probpost[1]))
      Evitessetmp   <- rbind(Evitessetmp,   c(unlist(results[[j]]$F.effects["vitesse",   3:5]), results[[j]]$probpost[2]))
      Etrafictmp    <- rbind(Etrafictmp,    c(unlist(results[[j]]$F.effects["trafic",    3:5]), results[[j]]$probpost[3]))
      Edisthydrotmp <- rbind(Edisthydrotmp, c(unlist(results[[j]]$F.effects["disthydro", 3:5]), results[[j]]$probpost[4]))
      Edistvegtmp   <- rbind(Edistvegtmp,   c(unlist(results[[j]]$F.effects["distveg",   3:5]), results[[j]]$probpost[5]))
    } else {
      Efaunetmp     <- rbind(Efaunetmp,     rep(NA, 4))
      Evitessetmp   <- rbind(Evitessetmp,   rep(NA, 4))
      Etrafictmp    <- rbind(Etrafictmp,    rep(NA, 4))
      Edisthydrotmp <- rbind(Edisthydrotmp, rep(NA, 4))
      Edistvegtmp   <- rbind(Edistvegtmp,   rep(NA, 4))
    }
  }
  rm(results)

  colnames(Efaunetmp)[4]     <- "probpost"
  colnames(Evitessetmp)[4]   <- "probpost"
  colnames(Etrafictmp)[4]    <- "probpost"
  colnames(Edisthydrotmp)[4] <- "probpost"
  colnames(Edistvegtmp)[4]   <- "probpost"

  Efaune[[il]]     <- as.data.frame(Efaunetmp)
  Evitesse[[il]]   <- as.data.frame(Evitessetmp)
  Etrafic[[il]]    <- as.data.frame(Etrafictmp)
  Edisthydro[[il]] <- as.data.frame(Edisthydrotmp)
  Edistveg[[il]]   <- as.data.frame(Edistvegtmp)
}

names(Efaune) <- names(Evitesse) <- names(Etrafic) <-
  names(Edisthydro) <- names(Edistveg) <- unique(listeSP[, "grpSP"])

save(Efaune, Evitesse, Etrafic, Edisthydro, Edistveg, file = "estim_effects.rda")
load("estim_effects.rda")

# Helper to create filled band polygons from quantiles for each species
make_ci_polygon_df <- function(d) {
  stopifnot(all(c("0.025quant","0.5quant","0.975quant") %in% names(d)))
  n <- nrow(d)
  data.frame(
    x = c(d$`0.975quant`, rev(d$`0.025quant`)),
    y = c(seq_len(n),     rev(seq_len(n))) / 100  # scaled to (0,1)
  )
}

# Build polygon coordinates & effect summaries for all parameters
all_sp_polys <- list(
  lapply(Efaune,     make_ci_polygon_df),
  lapply(Evitesse,   make_ci_polygon_df),
  lapply(Etrafic,    make_ci_polygon_df),
  lapply(Edisthydro, make_ci_polygon_df),
  lapply(Edistveg,   make_ci_polygon_df)
)

all_sp <- list(
  lapply(Efaune, function(x) {
    qEff <- quantile(x$`0.5quant`, probs = c(0.025, 0.5, 0.975))
    proba <- mean(x$probpost, na.rm = TRUE); c(qEff, proba)
  }),
  lapply(Evitesse, function(x) {
    qEff <- quantile(x$`0.5quant`, probs = c(0.025, 0.5, 0.975))
    proba <- mean(x$probpost, na.rm = TRUE); c(qEff, proba)
  }),
  lapply(Etrafic, function(x) {
    qEff <- quantile(x$`0.5quant`, probs = c(0.025, 0.5, 0.975))
    proba <- mean(x$probpost, na.rm = TRUE); c(qEff, proba)
  }),
  lapply(Edisthydro, function(x) {
    qEff <- quantile(x$`0.5quant`, probs = c(0.025, 0.5, 0.975))
    proba <- mean(x$probpost, na.rm = TRUE); c(qEff, proba)
  }),
  lapply(Edistveg, function(x) {
    qEff <- quantile(x$`0.5quant`, probs = c(0.025, 0.5, 0.975))
    proba <- mean(x$probpost, na.rm = TRUE); c(qEff, proba)
  })
)

common_names <- names(all_sp[[1]])
EffectsMORT <- lapply(common_names, function(name) {
  do.call(rbind, lapply(all_sp, function(lst) lst[[name]]))
})

names(EffectsMORT) <- c("Wild boar","Red fox","Other mustelids",
                        "European badger","Roe deer","European hare","Red deer")

EffectsMORT <- lapply(EffectsMORT, function(x) {
  rownames(x) <- c("Exposure","Speed","Traffic","Distance to water","Distance to vegetation")
  colnames(x) <- c("0.025quant","median","0.975quant","Posterior probability")
  x
})

# Tidy polygons for plotting
if (is.null(names(all_sp_polys)) || any(names(all_sp_polys) == "")) {
  names(all_sp_polys) <- paste0("Param ", seq_along(all_sp_polys))
}

poly_df <- purrr::imap_dfr(all_sp_polys, function(sp_list, param_name) {
  purrr::imap_dfr(sp_list, function(df_xy, species_name) {
    df_xy %>% mutate(Species = species_name, Parameter = param_name)
  })
})

poly_df <- poly_df %>%
  mutate(Species = forcats::fct_rev(factor(Species)),
         Parameter = factor(Parameter, levels = names(all_sp_polys)))

param_labels   <- c("Exposure","Speed","Traffic","Distance to water","Distance to vegetation")
species_labels <- c("Wild boar","Red fox","Other mustelids","European hare","Roe deer","Red deer","European badger")

poly_df <- poly_df %>%
  mutate(Parameter = factor(Parameter, labels = param_labels),
         Species   = factor(Species,   labels = species_labels))

# Vertical offsets so each species sits on a categorical row
band <- 0.8
sp_levels <- levels(poly_df$Species)
sp_index  <- setNames(seq_along(sp_levels), sp_levels)

poly_df <- poly_df %>%
  group_by(Species, Parameter) %>%
  mutate(y_std = (y - min(y, na.rm = TRUE)) /
                 (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) - 0.5) %>%
  ungroup() %>%
  mutate(sp_row = sp_index[as.character(Species)],
         y_plot = sp_row + band * y_std)

# Posterior-probability bins for fill/color
pp_df <- do.call(rbind, lapply(names(EffectsMORT), function(sp) {
  df <- EffectsMORT[[sp]]
  data.frame(
    Species = sp,
    Parameter = if (!is.null(rownames(df))) rownames(df) else rownames(df <- transform(df, Parameter = rownames(df))),
    PosteriorProbability = df[, "Posterior probability", drop = TRUE],
    row.names = NULL
  )
}))

bin_colors <- c(
  "[0–0.05)"   = "#E41A1C",
  "[0.05–0.1)" = darken("firebrick", amount = 0.3),
  "[0.1–0.9)"  = "black",
  "[0.9–0.95)" = darken("steelblue", amount = 0.3),
  "[0.95–1]"   = "#1E90FF"
)

pp_df <- pp_df %>%
  mutate(pp_bin = cut(PosteriorProbability,
                      breaks = c(0, 0.05, 0.1, 0.9, 0.95, 1),
                      include.lowest = TRUE, right = TRUE,
                      labels = names(bin_colors)))

poly_df <- poly_df %>% left_join(pp_df, by = c("Species", "Parameter"))

# Order and facet labels
poly_df$Species <- factor(poly_df$Species, levels = rev(unique(poly_df$Species)))
poly_df$Parameter <- factor(poly_df$Parameter,
                            levels = c("Exposure","Speed","Traffic","Distance to water","Distance to vegetation"),
                            labels = c("(a) Exposure","(b) Speed","(c) Traffic","(d) Distance to water","(e) Distance to vegetation"))

# Plot
p <- ggplot(poly_df, aes(x = x, y = y_plot, group = interaction(Species, Parameter))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_polygon(aes(fill = pp_bin, color = pp_bin), alpha = 0.35, linewidth = 0.2) +
  facet_wrap(~ Parameter, scales = "free_x", ncol = 3) +
  scale_y_continuous(breaks = sp_index, labels = names(sp_index),
                     expand = expansion(mult = c(0.05, 0.08))) +
  labs(x = "Estimate (x)", y = "Species",
       fill = "Posterior\nProbability", color = "Posterior Probability\n(effect > 0)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.title.x       = element_text(margin = margin(t = 10)))

if ("pp_bin" %in% names(poly_df)) {
  p <- p + scale_fill_manual(values = bin_colors) +
           scale_color_manual(values = bin_colors)
}

legend <- get_legend(p + theme(legend.position = "right"))
p_nolegend <- p + theme(legend.position = "none")
final_plot <- p_nolegend + inset_element(legend, left = 0.7, bottom = 0.05, right = 0.8, top = 0.2)

png(file = "effectMORT.png", width = 300, height = 200, units = "mm", res = 200)
final_plot
dev.off()


# -----------------------------------------------------------------------------
# 14) Session info for reproducibility (add to supplement)
# -----------------------------------------------------------------------------
# capture.output(sessionInfo(), file = "SESSION_INFO.txt")
# writeLines(capture.output(sessionInfo()), "SESSION_INFO.txt")
