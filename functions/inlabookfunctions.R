# spde_utils_commented.R
# -----------------------------------------------------------------------------
# Title: Utilities for simulating and visualising SPDE-based Gaussian fields
# Author: <Your name>
# Date: <YYYY-MM-DD>
#
# Description
# -----------
# This file collects helper functions for simulating Matérn Gaussian fields,
# sampling SPDE fields on a mesh, constructing dual (Voronoi) cells from an
# INLA mesh, visualising sparse matrices and spatial fields, simple colour
# utilities, and constructing barrier polygons in parallel.
#
# Provenance & Acknowledgements
# -----------------------------
# The original versions of several functions are adapted from the SPDE book by
# Virgilio Gómez-Rubio ("INLA and SPDEs" GitBook):
#   https://becarioprecario.bitbucket.io/spde-gitbook/index.html
# Please cite that resource (and the Lindgren, Rue & Lindström 2011 paper) in
# your manuscript when using this code. Modifications here include extended
# documentation, minor input checks, comments, and parallel variants.
#
# Dependencies
# ------------
# - INLA / R-INLA (inla.mesh.*, inla.spde2.*, inla.qsample)
# - inlabru/fmesher (optional, for some mesh tools)
# - Matrix (for dgTMatrix class used by plot.dgTMatrix)
# - fields (image.plot), viridisLite, RColorBrewer
# - sp and sf (geometries), parallel (clusters)
#
# Licensing
# ---------
# If you intend to redistribute, ensure that this file and your manuscript
# comply with licenses of upstream sources (SPDE GitBook and package licenses).
# -----------------------------------------------------------------------------

# Packages used directly below
library(fields)       # image.plot
library(viridisLite)  # viridis(), magma()
library(RColorBrewer) # brewer.pal()
# Note: Functions below also use INLA, sp, sf, Matrix, parallel; they are loaded
# lazily within functions to avoid hard dependencies at import time.

# -----------------------------------------------------------------------------
# Simulate Matérn fields using the covariance function
# -----------------------------------------------------------------------------
#' Simulate a Matérn Gaussian random field at given coordinates
#'
#' @param n Integer. Number of field realisations to generate.
#' @param coords Numeric matrix (n_locations x 2) of coordinates.
#' @param sigma Numeric. Standard deviation parameter (\u03C3); see also `variance`.
#' @param range Numeric. Practical correlation range parameter.
#' @param kappa Numeric. Equivalent SPDE scale parameter; default linked to
#'   `range` and smoothness `nu` via kappa = sqrt(8*nu)/range.
#' @param variance Numeric. Field variance (\u03C3^2). Default `sigma^2`.
#' @param nu Numeric. Matérn smoothness parameter.
#'
#' @return A numeric matrix with `nrow(coords)` rows and `n` columns containing
#'   realisations of the field.
#' @details Uses the closed-form Matérn covariance (with Bessel K) to assemble
#'   the covariance matrix, then samples via the Cholesky factor.
book.rMatern <- function(n, coords, sigma = 1, range,
                         kappa = sqrt(8*nu)/range,
                         variance = sigma^2, nu = 1) {
  # Pairwise Euclidean distances
  m <- as.matrix(dist(coords))
  # Matérn covariance kernel (no nugget). For zero distance, we set diag=1 below.
  m <- exp((1 - nu) * log(2) + nu * log(kappa * m) - lgamma(nu)) * besselK(m * kappa, nu)
  diag(m) <- 1
  # Draw n realisations using Cholesky factorisation
  return(drop(crossprod(chol(variance * m), matrix(rnorm(nrow(coords) * n), ncol = n))))
}

# -----------------------------------------------------------------------------
# Simulate SPDE fields on an INLA mesh and project to arbitrary coordinates
# -----------------------------------------------------------------------------
#' Sample a Gaussian field from an SPDE (Matérn) model and project to `coords`
#'
#' @param coords Numeric matrix (n_locations x 2) of coordinates for projection.
#' @param sigma,range,variance As in `book.rMatern`.
#' @param alpha Integer >= 1. SPDE order; alpha=2 corresponds to Matérn with nu=1.
#' @param kappa Numeric. SPDE scale; default linked to `alpha` and `range`.
#' @param n Integer. Number of field realisations to draw.
#' @param mesh Optional `inla.mesh` object. If missing, a mesh is created
#'   heuristically around the convex hull of `coords`.
#' @param verbose Logical. Print diagnostic information.
#' @param seed Optional integer seed passed to INLA's sampler.
#' @param return.attributes Logical. If TRUE, attach mesh, SPDE, Q, A, timings.
#'
#' @return A numeric vector/matrix of projected field values at `coords`.
#' @note Requires R-INLA. The internal mesh creation line mirrors the original
#'   book code; check syntax and adjust to your R-INLA version if needed.
book.rspde <- function(coords, sigma = 1, range,
                       variance = sigma^2, alpha = 2,
                       kappa = sqrt(8*(alpha - 1))/range,
                       n = 1, mesh,
                       verbose = FALSE, seed, return.attributes = FALSE) {
  # Timing start
  t0 <- Sys.time()

  # Parameterisation used by inla.spde2.precision()
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (verbose) cat('theta =', theta, '\n')

  # Mesh handling: use provided mesh or build a quick mesh around coords
  if (missing(mesh)) {
    # Heuristic mesh edge/cutoff parameters akin to the book defaults
    mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords) / 2) / kappa
    if (verbose) cat('mesh.pars =', mesh.pars, '\n')
    # NOTE: The original line in the book may read differently depending on INLA
    # versions. The trailing comma in `inla.mesh.2d(, coords[...])` in some
    # snippets is a typo. Replace with a proper call for your setup, e.g.:
    #   inla.mesh.2d(loc = coords[chull(coords), ], max.edge = mesh.pars[1:2],
    #                cutoff = mesh.pars[3], offset = mesh.pars[4:5])
    attributes <- list(
      mesh = inla.mesh.2d(,
                          coords[chull(coords), ],
                          max.edge = mesh.pars[1:2],
                          cutoff   = mesh.pars[3],
                          offset   = mesh.pars[4:5])
    )
    if (verbose) cat('n.mesh =', attributes$mesh$n, '\n')
  } else {
    attributes <- list(mesh = mesh)
  }

  # SPDE model and precision matrix
  attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q    <- inla.spde2.precision(attributes$spde, theta = theta)
  # Projector from mesh to requested coordinates
  attributes$A    <- inla.mesh.project(mesh = attributes$mesh, loc = coords)$A

  # If n==1, we could sample once and project; kept for backward-compatibility.
  if (n == 1)
    result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                 constr = attributes$spde$f$extraconstr))
  t1 <- Sys.time()

  # General sampling path (overwrites the above if executed)
  result <- inla.qsample(n, attributes$Q,
                         seed   = ifelse(missing(seed), 0, seed),
                         constr = attributes$spde$f$extraconstr)

  # Ensure the returned object has one row per coordinate after projection
  if (nrow(result) < nrow(attributes$A)) {
    # Pad and project each draw
    result <- rbind(result, matrix(NA, nrow(attributes$A) - nrow(result), ncol(result)))
    dimnames(result)[[1]] <- paste('x', 1:nrow(result), sep = '')
    for (j in 1:ncol(result))
      result[, j] <- drop(attributes$A %*% result[1:ncol(attributes$A), j])
  } else {
    for (j in 1:ncol(result))
      result[1:nrow(attributes$A), j] <- drop(attributes$A %*% result[, j])
    result <- result[1:nrow(attributes$A), ]
  }
  t2 <- Sys.time()
  attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 - t0)

  if (return.attributes)
    attributes(result) <- c(attributes(result), attributes)

  return(drop(result))
}

# -----------------------------------------------------------------------------
# Dual mesh (Voronoi cells) construction with parallelisation
# -----------------------------------------------------------------------------
#' Construct dual (Voronoi-like) polygons for an INLA 2D mesh
#'
#' @param mesh An `inla.mesh` object on manifold 'R2'.
#' @param Ncores Integer. Number of workers for parallel polygon construction.
#'
#' @return A `sp::SpatialPolygons` object with one polygon per mesh node.
#' @details Parallel variant of the original helper; uses `parallel::parLapply`.
book.mesh.dual <- function(mesh, Ncores = 1) {
  if (mesh$manifold != 'R2') stop("It only works for R2!")

  # Triangle centroids for later polygon assembly
  ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
    colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))

  # Create a PSOCK cluster
  parallel::cluster <- NULL
  cl <- parallel::makeCluster(Ncores)

  # Export necessary data and attach packages on workers
  parallel::clusterExport(cl, varlist = c("mesh", "ce"), envir = environment())
  parallel::clusterEvalQ(cl, { library(sp) })

  # Build one polygon per mesh node in parallel
  pls <- parallel::parLapply(cl, 1:mesh$n, function(i) {
    p <- unique(Reduce('rbind', lapply(1:3, function(k) {
      j <- which(mesh$graph$tv[, k] == i)
      if (length(j) > 0) {
        return(rbind(
          ce[j, , drop = FALSE],
          cbind(
            (mesh$loc[mesh$graph$tv[j, k], 1] + mesh$loc[mesh$graph$tv[j, c(2, 3, 1)[k]], 1]) / 2,
            (mesh$loc[mesh$graph$tv[j, k], 2] + mesh$loc[mesh$graph$tv[j, c(2, 3, 1)[k]], 2]) / 2
          )
        ))
      } else {
        return(ce[j, , drop = FALSE])
      }
    })))

    # Boundary handling
    j1 <- which(mesh$segm$bnd$idx[, 1] == i)
    j2 <- which(mesh$segm$bnd$idx[, 2] == i)

    if ((length(j1) > 0) | (length(j2) > 0)) {
      p <- unique(rbind(
        mesh$loc[i, 1:2], p,
        mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2,
        mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2
      ))
      yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
      xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
    } else {
      yy <- p[, 2] - mesh$loc[i, 2]
      xx <- p[, 1] - mesh$loc[i, 1]
    }

    sp::Polygon(p[order(atan2(yy, xx)), ])
  })

  # Clean up
  parallel::stopCluster(cl)

  # Wrap in SpatialPolygons (IDs as character)
  return(sp::SpatialPolygons(lapply(1:mesh$n, function(i)
    sp::Polygons(list(pls[[i]]), ID = as.character(i)))))
}

# -----------------------------------------------------------------------------
# Colour utilities
# -----------------------------------------------------------------------------
#' Generate a custom diverging colour ramp (red/green/blue families)
#'
#' @param n Integer (> 2). Number of colours.
#' @param type One of 'red', 'green', 'blue'.
#' @param u Optional numeric in [0,1] at which to evaluate the ramp.
#' @return A vector of RGB colours (character hex strings).
genColor <- function(n, type = c('red', 'green', 'blue'), u = NULL) {
  cbp <- list(
    red = list(c(255, 254, 252, 252, 251, 239, 203, 165, 103),
               c(245, 224, 187, 146, 106, 59, 24, 15, 0),
               c(240, 210, 161, 114, 74, 44, 29, 21, 13)),
    green = list(c(247, 229, 199, 161, 116, 65, 35, 0, 0),
                 c(252, 245, 233, 217, 196, 171, 139, 109, 68),
                 c(245, 224, 192, 155, 118, 93, 69, 44, 27)),
    blue = list(c(247, 222, 198, 158, 107, 66, 33, 8, 8),
                c(251, 235, 219, 202, 174, 146, 113, 81, 48),
                c(255, 247, 239, 225, 214, 198, 181, 156, 107)))
  if (n < 2) stop("Works for 'n>2'!")
  if (is.null(u)) u <- 0:(n - 1) / (n - 1)
  u0 <- 0:8 / 8
  i <- findInterval(u, u0, TRUE)
  k <- pmatch(match.arg(type), c('red', 'green', 'blue'))
  w1 <- 8 * (u0[i + 1] - u) / 255; w2 <- 8 * (u - u0[i]) / 255
  grDevices::rgb(cbp[[k]][[1]][i] * w1 + cbp[[k]][[1]][i + 1] * w2,
                 cbp[[k]][[2]][i] * w1 + cbp[[k]][[2]][i + 1] * w2,
                 cbp[[k]][[3]][i] * w1 + cbp[[k]][[3]][i + 1] * w2)
}

# -----------------------------------------------------------------------------
# Visualisation helpers
# -----------------------------------------------------------------------------
#' Plot method for sparse matrices (dgTMatrix) with a colour scale
#'
#' @param x A Matrix::dgTMatrix object.
#' @param y Ignored.
#' @param ... Additional arguments passed to `image`.
#' @details Rounds values to `digits` unique levels and builds a diverging scale.
plot.dgTMatrix <- function(x, y, ...) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Matrix package required for dgTMatrix.")
  cl <- match.call()
  digits <- if (is.null(cl$digits)) 2 else eval(cl$digits)
  z <- sort(unique(round(x@x, digits)))
  nz <- length(z)
  n1 <- sum(z < 0)
  n2 <- sum(z > 0)
  if (is.null(cl$colors)) {
    if (any(c(n1, n2) == 0)) {
      colors <- gray(0.9 * (1 - (z - min(z)) / diff(range(z))))
    } else {
      colors <- c(genColor(n1, 'red',  z[z < 0] / min(z)),
                  rep('white', nz - n1 - n2),
                  genColor(n2, 'blue', z[z > 0] / max(z)))
    }
  } else {
    colors <- eval(cl$colors)
  }
  z.breaks <- c(z[1] - diff(z[1:2]) / 2,
                z[-nz] / 2 + z[-1] / 2,
                z[nz] + diff(z[nz - 1:0]) / 2)
  x@x <- round(x@x, digits)
  image(x, at = z.breaks, col.regions = colors, ...)
}

#' Project and plot a spatial field using an INLA projector
#'
#' @param field Vector of values at mesh nodes, or a list with x,y,z as for `image`.
#' @param mesh Optional INLA mesh; if provided and `projector` missing, one is built.
#' @param projector Optional projector list from `inla.mesh.projector`.
#' @param xlim,ylim Plot extents.
#' @param dims Dimensions passed to projector creation when `mesh` is given.
#' @param poly Optional polygon(s) to overlay.
#' @param asp,axes,xlab,ylab,col Graphical parameters.
book.plot.field <- function(field, mesh, projector, xlim, ylim,
                            dims = c(300, 300), poly, asp = 1,
                            axes = FALSE, xlab = '', ylab = '',
                            col = book.color.c(), ...) {
  # If no mesh, either plot provided matrix/list directly or project via projector
  if (missing(mesh)) {
    if (missing(projector)) {
      if (missing(xlim) | missing(ylim)) {
        fields::image.plot(field, asp = asp, axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
      } else {
        fields::image.plot(field, xlim = xlim, ylim = ylim, asp = asp, axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
      }
    } else {
      if (missing(xlim)) xlim <- range(projector$x)
      if (missing(ylim)) ylim <- range(projector$y)
      field.proj <- inla.mesh.project(projector, field)
      fields::image.plot(x = projector$x, y = projector$y, z = field.proj,
                         asp = asp, axes = axes, xlab = xlab, ylab = ylab,
                         col = col, xlim = xlim, ylim = ylim, ...)
    }
  } else {
    # Build projector from mesh and project the field
    if (missing(xlim)) xlim <- range(mesh$loc[, 1])
    if (missing(ylim)) ylim <- range(mesh$loc[, 2])
    projector <- inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, dims = dims)
    field.proj <- inla.mesh.project(projector, field)
    fields::image.plot(x = projector$x, y = projector$y, z = field.proj,
                       asp = asp, axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
  }
  if (!missing(poly)) plot(poly, add = TRUE, col = 'grey')
}

# -----------------------------------------------------------------------------
# Barrier models: correlation and polygon construction
# -----------------------------------------------------------------------------
#' Compute spatial correlation from a precision matrix at a given location
#'
#' @param Q Sparse precision matrix (e.g., from `inla.spde2.precision`).
#' @param location Numeric length-2 vector (x,y).
#' @param mesh INLA mesh used to define `Q`.
#' @return Numeric vector of correlations at mesh nodes with the closest node to
#'   `location`.
book.spatial.correlation <- function(Q, location, mesh) {
  # Marginal SDs from the (approximate) covariance
  sd <- sqrt(diag(inla.qinv(Q)))

  # Identify nearest mesh node via a temporary A-matrix
  A.tmp <- inla.spde.make.A(mesh = mesh, loc = matrix(c(location[1], location[2]), 1, 2))
  id.node <- which.max(A.tmp[1, ])

  # Solve for the covariance column associated with id.node
  Inode <- rep(0, dim(Q)[1]); Inode[id.node] <- 1
  covar.column <- solve(Q, Inode)
  corr <- drop(matrix(covar.column)) / (sd * sd[id.node])
  return(corr)
}

# Colour palettes for continuous/discrete scales --------------------------------
book.color.c  <- function(n = 201) viridis(n)
book.color.c2 <- function(n = 201) magma(n)
book.color.dc <- function(n = 11)  viridis(n)
book.color.d  <- function(n = 4)   brewer.pal(n = n, name = "Paired")

# -----------------------------------------------------------------------------
# Barrier polygon construction in parallel
# -----------------------------------------------------------------------------
#' Build barrier polygons (Omega, Gamma) for inla.barrier models in parallel
#'
#' @param mesh INLA mesh.
#' @param barrier.triangles Optional integer indices of triangles forming the
#'   barrier set; if provided, overrides `Omega`.
#' @param Omega Optional list of length 2: indices for non-barrier and barrier
#'   triangle sets (Omega and Gamma).
#' @param Ncores Integer. Number of workers.
#' @return If `barrier.triangles` missing: list of two `sp::SpatialPolygons` for
#'   Omega and Gamma. Otherwise, the barrier polygon (`Gamma`) only.
barrier_polygon_parallel <- function(mesh, barrier.triangles, Omega = NULL, Ncores = 1) {
  stopifnot(inherits(mesh, "inla.mesh"))

  if (!requireNamespace("sf", quietly = TRUE))
    stop("The 'sf' package is required but not installed.")

  # If barrier.triangles provided, compute Omega from it (overwrite user Omega)
  if (!missing(barrier.triangles)) {
    barrier.triangles <- unique(barrier.triangles)
    t <- nrow(mesh$graph$tv)
    remaining <- setdiff(1:t, barrier.triangles)
    if (!is.null(Omega)) warning("Omega is replaced by barrier.triangles")
    Omega <- list(remaining, barrier.triangles)
  }

  # Start cluster
  cl <- parallel::makeCluster(Ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Export and attach packages
  parallel::clusterExport(cl, varlist = c("mesh", "Omega"), envir = environment())
  parallel::clusterEvalQ(cl, { library(sp); library(sf) })

  # Build polygons per subset in parallel
  Omega.SP.list <- parallel::parLapply(cl, 1:length(Omega), function(j) {
    tri_ids <- Omega[[j]]

    # Build triangle polygons
    poly.list <- lapply(tri_ids, function(tri) {
      px <- mesh$graph$tv[tri, ]
      temp <- mesh$loc[px, ]
      sp::Polygon(rbind(temp[c(3, 2, 1), 1:2], temp[3, 1:2]), hole = FALSE)
    })
    mesh.polys <- sp::SpatialPolygons(list(sp::Polygons(poly.list, ID = "noid")))

    # Convert to sf and union to a single polygon set
    mesh.polys <- sf::st_as_sfc(mesh.polys)
    sf::as_Spatial(sf::st_union(mesh.polys))
  })

  if (missing(barrier.triangles)) {
    return(Omega.SP.list)
  } else {
    return(Omega.SP.list[[2]])
  }
}

# -----------------------------------------------------------------------------
# End of file
# -----------------------------------------------------------------------------
