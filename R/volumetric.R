# X must be centetered and rotated
.radius.cylinder <- function(X, side = c("inner", "outer")) {
  side <- match.arg(side)
  x <- X[, 1]
  R <- sqrt(X[, 2]^2 + X[, 3]^2)
  phi <- atan2(X[, 2], X[, 3])

  Xi <- switch(side,
               outer = X,
               inner = cbind(X[, 1], X[, 2] / R^2, X[, 3] / R^2))
  Xi.with.borders <- rbind(Xi, c(min(X[, 1]), 0, 0), c(max(X[, 1]), 0, 0))

  hull <- convhulln(Xi.with.borders, options = "Pp")

  idx.with.borders <- unique(as.vector(hull))
  idx <- setdiff(idx.with.borders,
                 c(nrow(Xi.with.borders), nrow(Xi.with.borders) - 1))

  # TODO use proper interpolation using hull-induced tesselation instead of naive one
  x.base <- c(min(x), max(x), x[idx])
  phi.base <- c(0, 0, phi[idx])
  R.base <- c(0, 0, R[idx])

  xphi.base <- cbind(x.base, c(phi.base, phi.base + 2*pi, phi.base - 2*pi))
  R.base <- rep(R.base, 3)

  R <- linear.interpolate(cbind(x, phi),
                          xphi.base,
                          R.base)

  R
}

# X must be centetered and rotated
.radius.sphere <- function(X, side = c("inner", "outer")) {
  side <- match.arg(side)
  R <- sqrt(X[, 1]^2 + X[, 2]^2 + X[, 3]^2)

  Rt <- 1 / X[, 3] # MB use atan2(R, X[, 3]) / sqrt(X[, 1]^2 + sqrt(X[, 2]^2))
  x <- X[, 1] * Rt
  y <- X[, 2] * Rt

  Xi <- switch(side,
               outer = X,
               inner = X / R^2)
  Xi.with.borders <- rbind(Xi, c(0, 0, 0))

  hull <- convhulln(Xi.with.borders, options = "Pp")

  idx.with.borders <- unique(as.vector(hull))
  idx <- setdiff(idx.with.borders, nrow(Xi.with.borders))

  # TODO use proper interpolation using hull-induced tesselation instead of naive one
  x.base <- x[idx]
  y.base <- y[idx]
  R.base <- R[idx]

  xy.base <- cbind(x.base, y.base)

  R <- linear.interpolate(cbind(x, y),
                          xy.base,
                          R.base)

  R
}

.unfold3d.sphere <- function(X) {
  # X <- rotate.sphere(X) must be done

  R <- sqrt(X[, 1]^2 + X[, 2]^2 + X[, 3]^2)
  Rt <- 1 / X[, 3] # MB use atan2(R, X[, 3]) / sqrt(X[, 1]^2 + sqrt(X[, 2]^2))
  x <- X[, 1] * Rt
  y <- X[, 2] * Rt

  R.inner <- .radius.sphere(X, side = "inner")
  R.outer <- .radius.sphere(X, side = "outer")

  depth <- (R - R.inner) / (R.outer - R.inner)
  depth[depth > 1] <- 1
  depth[depth < 0] <- 0

  cbind(x = x, y = y, depth = depth) # All of them are normalized
}

.unfold3d.cylinder <- function(X) {
  # X <- rotate.cylinder(X) must be done

  R <- sqrt(X[, 2]^2 + X[, 3]^2)
  x <- X[, 1]
  phi <- atan2(X[, 2], X[, 3])

  R.inner <- .radius.cylinder(X, side = "inner")
  R.outer <- .radius.cylinder(X, side = "outer")

  depth <- (R - R.inner) / (R.outer - R.inner)
  depth[depth > 1] <- 1
  depth[depth < 0] <- 0

  cbind(x = x, depth = depth, phi = phi) # phi is 2pi-periodic, x is NOT normalized
}

interpolate2grid.embryo3d.cylinder <- function(x, ...,
                                               cuts = c(x = 200, depth = 10, phi = 200)) {
  stopifnot(.is.unfolded(x))

  uX <- cbind(x = x$x, depth = x$depth, phi = x$phi)
  v <- x$values

  # Omit NAs FIXME
  mask <- !is.na(rowSums(uX))

  uX <- uX[mask,, drop = FALSE]
  v <- v[mask]

  eps <- 1e-5
  ox <- seq(min(uX[, "x"]) + eps, max(uX[, "x"]) - eps, length.out = cuts["x"])
  oy <- seq(min(uX[, "depth"]) + eps, max(uX[, "depth"]) - eps, length.out = cuts["depth"])
  oz <- seq(0, 2*pi, length.out = cuts["phi"] + 1); oz <- oz[-length(oz)]

  grid <- as.matrix(expand.grid(x = ox, y = oy, z = oz))

  uX.up <- uX.down <- uX
  uX.up[, "phi"] <- uX[, "phi"] + 2*pi
  uX.down[, "phi"] <- uX[, "phi"] - 2*pi

  f <- linear.interpolate(grid, rbind(uX, uX.down, uX.up), rep(v, 3))

  dim(f) <- sapply(list(ox, oy, oz), length)
  field <- list(x = ox, y = oy, z = oz, f = f)

  x$field <- field

  x
}

interpolate2grid.embryo3d.sphere <- function(x, ...,
                                             cuts = c(x = 200, y = 200, depth = 10)) {
  stopifnot(.is.unfolded(x))

  uX <- cbind(x = x$x, y = x$y, depth = x$depth)
  v <- x$values

  # Omit NAs
  mask <- !is.na(rowSums(uX))

  uX <- uX[mask,, drop = FALSE]
  v <- v[mask]

  eps <- 1e-5
  ox <- seq(min(uX[, 1]) + eps, max(uX[, 1]) - eps, length.out = cuts["x"])
  oy <- seq(min(uX[, 2]) + eps, max(uX[, 2]) - eps, length.out = cuts["y"])
  oz <- seq(min(uX[, 3]) + eps, max(uX[, 3]) - eps, length.out = cuts["depth"])

  grid <- as.matrix(expand.grid(x = ox, y = oy, z = oz))

  f <- linear.interpolate(grid, uX, v)

  dim(f) <- sapply(list(ox, oy, oz), length)
  field <- list(x = ox, y = oy, z = oz, f = f)

  x$field <- field

  x
}

.is.interpolated <- function(x, ...) {
  !is.null(x["field"])
}

.is.unfolded <- function(x, ...) {
  !is.null(x["x"]) || !is.null(x["x2d"])
}

update.field <- function(x, ...)
  UseMethod("update.field")

update.field.embryo3d <- function(x, ...)
  stop("Abstract method")


shrink <- function(what, to) {
  to <- range(to)
  eps <- 1e-5
  what[what <= to[1]] <- to[1] + eps
  what[what >= to[2]] <- to[2] - eps

  what
}

update.field.embryo3d.sphere <- function(x, newvalues = x$field$f, ...) {
  stopifnot(.is.interpolated(x))

  x$field$f[] <- as.numeric(newvalues)
  ox <- x$x; ox <- shrink(ox, x$field$x)
  oy <- x$y; oy <- shrink(oy, x$field$y)
  oz <- x$depth; oz <- shrink(oz, x$field$z)

  x$values <- approx3d(x$field$x, x$field$y, x$field$z, x$field$f,
                       ox, oy, oz)

  x
}

approx3d.cycled <- function(x, y, z, f, xout, yout, zout) {
  z <- c(z - 2*pi, z, z + 2*pi)
  f <- rep(f, 3); dim(f) <- sapply(list(x, y, z), length)

  zout <- zout - 2*pi * floor(zout / (2*pi))

  approx3d(x, y, z, f, xout, yout, zout)
}

update.field.embryo3d.cylinder <- function(x, newvalues = x$field$f, ...) {
  stopifnot(.is.interpolated(x))

  x$field$f[] <- as.numeric(newvalues)

  ox <- x$x; ox <- shrink(ox, x$field$x)
  oy <- x$depth; oy <- shrink(oy, x$field$y)
  oz <- x$phi; oz <- shrink(oz, x$field$z)

  x$values <- approx3d.cycled(x$field$x, x$field$y, x$field$z, x$field$f,
                              ox, oy, oz)

  x
}


find.center.sphere <- function(X, initial = rep(0, 3), ..., degree = 1) {
  d2 <- function(center)
    (X[, 1] - center[1])^2 + (X[, 2] - center[2])^2 + (X[, 3] - center[3])^2

  optim(initial, function(center) sd(d2(center)^degree), ..., method = "BFGS")$par
}

rotate.sphere <- function(X, center = find.center.sphere(X)) {
  X <- scale(X, center = center, scale = FALSE)
  cX <- colMeans(X)
  cX <- cX / sqrt(sum(cX^2))

  M <- cbind(cX, c(1, 0, 0), c(0, 1, 0))
  U <- qr.Q(qr(M))
  U <- U[, c(2, 3, 1)]

  if (U[2, 2] < 0)
    U[, 2] <- -U[, 2]

  if (det(U) < 0)
    U[, 2] <- -U[, 2]

  X %*% U
}

unfold.embryo3d.sphere <- function(x, ...) {
  X <- cbind(x$x3d, x$y3d, x$z3d)

  X <- rotate.sphere(X)

  uX <- .unfold3d.sphere(X)
  x$x <- uX[, "x"]
  x$y <- uX[, "y"]
  x$depth <- uX[, "depth"]

  class(x) <- c("embryo3d.sphere", "embryo3d")
  x
}

unfold.embryo3d.cylynder <- function(x, ...) {
  X <- cbind(x$x3d, x$y3d, x$z3d)

  X <- rotate.cylinder(X)

  uX <- .unfold3d.cylinder(X)
  x$x <- uX[, "x"]
  x$phi <- uX[, "phi"]
  x$depth <- uX[, "depth"]

  class(x) <- c("embryo3d.cylinder", "embryo3d")
  x
}

BioSSAv <- function(x, ...) {
  UseMethod("BioSSAv")
}

interpolate2grid <- function(x, ...)
  UseMethod("interpolate2grid")

BioSSAv.formula <- function(x, data = NULL, ...,
                            cuts = c(x = 100, y = 100, phi = 100, depth = 10),
                            kind = c("sphere", "cylinder")) {
  kind <- match.arg(kind)

  emb3 <- embryo3d(x, data = data)
  emb3 <- switch(kind,
                 sphere = unfold.embryo3d.sphere(emb3),
                 cylinder = unfold.embryo3d.cylynder(emb3))

  emb3 <- interpolate2grid(emb3, cuts = cuts)

  BioSSAv(emb3, ...)
}

BioSSAv.embryo3d <- function(x, L, ...) {
  stopifnot(.is.interpolated(x))
  f <- x$field$f

  if (inherits(x, "embryo3d.sphere")) {
    circular <- c(FALSE, FALSE, FALSE)
  } else if (inherits(x, "embryo3d.cylinder")) {
    circular <- c(FALSE, FALSE, TRUE)
  }

  dec <- ssa(f, L = L, ..., kind = "nd-ssa", circular = circular)

  res <- list(emb3 = x,
              ssa = dec)
  class(res) <- "BioSSAv" # "v" for volumetric =)

  res
}

decompose.BioSSAv <- function(x, ...) {
  x$ssa <- decompose(x$ssa, ...)

  x
}

reconstruct.BioSSAv <- function(x, groups, ...) {
  rec <- reconstruct(x$ssa, groups = groups, ...)
  res <- lapply(rec,
                function(component) {
                  update.field(x$emb3, newvalues = component)
                })
  attr(res, "series") <- x$emb3

  residuals <- update.field(x$emb3, attr(rec, "residuals"))
  residuals$values <- x$emb3$values
  for (j in seq_along(res)) {
    residuals$values <- residuals$values - res[[j]]$values
  }
  attr(res, "residuals") <- residuals

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSAv.reconstruction"
  invisible(res)
}


residuals.BioSSAv <- residuals.BioSSA2d
residuals.BioSSAv.reconstruction <- residuals.BioSSA2d.reconstruction
plot.BioSSAv <- plot.BioSSA2d





field.section.embryo3d <- function(emb3, slice = list(), units = c("percent", "original")) {
  stopifnot(.is.interpolated(emb3))

  units <- match.arg(units)

  field <- emb3$field

  # if (identical(units, "percent")) {
  #   field$x <- (field$x - emb2$x0perc) / emb2$x1perc
  #   field$y <- (field$y - emb2$y0perc) / emb2$y1perc
  # }

  # Determine free coord
  free.coord <- setdiff(names(field), c("f", names(slice)))

  slice.idx <- lapply(names(slice),
                      function(name) {
                        s <- slice[[name]]
                        coord <- field[[name]]
                        name <- names(s)
                        which.min(abs(coord - s))
                      })
  names(slice.idx) <- names(slice) # FIXME MB this is not needed
  stopifnot(length(slice.idx) == 2)
  tmp <- list(TRUE)
  names(tmp) <- free.coord
  slice.idx <- c(slice.idx, tmp)
  stopifnot(length(slice.idx) == 3)

  slice.idx <- slice.idx[order(names(slice.idx))]
  names(slice.idx) <- NULL

  values <- do.call("[", c(list(field$f, drop = TRUE), slice.idx))
  list(x = field[[free.coord]], y = values, free.coord = free.coord)
}

subset.embryo3d <- function(x, subset = list(), tolerance, ...) {
  x$field <- NULL

  mask <- rep(TRUE, length(x$x3d)) # TODO make proper nrow method
  for (name in names(subset)) {
    interval <- subset[[name]]
    coord <- x[[name]]

    if (length(interval) == 2) {
      interval <- sort(interval)
      cmask <- (coord >= interval[1]) & (coord <= interval[2])
      # TODO add support of circular topology
      # TODO add tolerance support
    } else {
      cmask <- coord %in% interval
    }

    mask <- mask & cmask
  }

  for (name in names(x)) {
    x[[name]] <- x[[name]][mask]
  }

  x
}
#
# nuclei.stripe.embryo3d <- function(emb3, slice = list(x = c(), y = c()),
#                                    units = c("percent", "original"), tolerance = 0.1, coord = "y") {
#   units <- match.arg(units)
#
#   if (inherits(emb2, "embryo3d")) emb2 <- embryo2d(emb2)
#
#   if (identical(units, "percent")) {
#     emb2$x2d <- (emb2$x2d - emb2$x0perc) / emb2$x1perc
#     emb2$topology[1] <- (emb2$topology[1] - emb2$x0perc) / emb2$x1perc
#     emb2$y2d <- (emb2$y2d - emb2$y0perc) / emb2$y1perc
#     emb2$topology[2] <- (emb2$topology[2] - emb2$y0perc) / emb2$y1perc
#   }
#
#   axe <- emb2[[paste0(coord, "2d")]]
#
#   along.coord.index <- if (coord == "x") 1 else 2
#   topology <- emb2$topology[along.coord.index]
#   idx <- (abs(axe - at) %% min(topology, 10e10)) < tolerance # TODO Get rid off this
#
#   out <- list(x2d = emb2$x2d[idx],
#               y2d = emb2$y2d[idx],
#               values = emb2$values[idx])
#
#   attr(out, "topology") <- attr(emb2, "topology")
#   attr(out, "topology")[-along.coord.index] <- FALSE
#
#   class(out) <- c("embryo2d", class(out))
#   out
# }
#
.plot1d.embryo3d.section.field <- function(x,
                                           slice,
                                           units = c("percent", "original"),
                                           ...,
                                           ref = FALSE) {
  dots <- list(...)
  units <- match.arg(units)

  section <- field.section.embryo3d(x, slice, units = units)
  free.coord <- section$free.coord

  dots <- .defaults(dots,
                    type = "l",
                    col = "black",
                    ylab = "Expression",
                    xlab = sprintf("Spatial coordinate: %s%s",
                                   free.coord,
                                   if (identical(units, "percent")) ", %" else ""))

  res <- do.call("xyplot", c(list(y ~ x, data = section), dots))

  if (ref) {
    res <- res + layer(panel.abline(h = 0, reference = TRUE))
  }

  res
}

.plot1d.embryo3d.section.nuclei <- function(x,
                                            slice,
                                            tolerance = 0.1,
                                            units = c("percent", "original"),
                                            ...,
                                            ref = FALSE) {
  dots <- list(...)
  units <- match.arg(units)

  stripe <- subset(x, subset = slice, tolerance = tolerance)
  free.coord <- setdiff(names(x), c(names(slice), "field", "values", "x3d", "y3d", "z3d"))

  dots <- .defaults(dots,
                    type = "p",
                    col = "green",
                    pch = 18,
                    ylab = "Expression",
                    xlab = sprintf("Spatial coordinate: %s%s",
                                   free.coord,
                                   if (identical(units, "percent")) ", %" else ""))

  stripe$xvalues <- stripe[[free.coord]]
  res <- do.call("xyplot", c(list(values ~ xvalues, data = stripe), dots))

  if (ref) {
    res <- res + layer(panel.abline(h = 0, reference = TRUE))
  }

  res
}


plot.embryo2d <- function(x, type = c("nuclei-2d", "field-2d",
                                      "field-section", "nuclei-section",
                                      "field-3d", "nuclei-3d"),
                          ...) {
  type <- match.arg(type)

  if (identical(type, "field-section")) {
    .plot1d.embryo2d.section.field(x, ...)
  } else if (identical(type, "nuclei-section")) {
    .plot1d.embryo2d.section.nuclei(x, ...)
  } else if (identical(type, "field-3d")) {
    .plot3d.embryo2d.field(x, ...)
  } else if (identical(type, "nuclei-3d")) {
    .plot3d.embryo2d.nuclei(x, ...)
  } else if (identical(type, "field-2d")) {
    .plot2d.embryo2d.field(x, ...)
  } else if (identical(type, "nuclei-2d")) {
    .plot2d.embryo2d.nuclei(x, ...)
  } else {
    stop("Unknown `type'")
  }
}
