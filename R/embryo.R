embryo2d <- function(x, ...)
  UseMethod("embryo2d")

embryo2d.formula <- function(x,
                             data = NULL,
                             topology = c(Inf, Inf),
                             ...) {
  mf <- model.frame(x, data = data, ..., na.action = NULL)

  emb2 <- list(x2d = mf[, 2],
               y2d = mf[, 3],
               values = mf[, 1])
  attr(emb2, "topology") <- rep(topology, 2)[1:2]

  class(emb2) <- "embryo2d"
  emb2
}

embryo3d <- function(x, ...)
  UseMethod("embryo3d")

embryo3d.formula <- function(x,
                             data = NULL,
                             ...) {
  mf <- model.frame(x, data = data, ..., na.action = NULL)

  emb3 <- list(x3d = mf[, 2],
               y3d = mf[, 3],
               z3d = mf[, 4],
               values = mf[, 1])

  class(emb3) <- "embryo3d"
  emb3
}

interpolate <- function(x, step) {
  # data == data.frame(x, y) + topology attr

  # TODO Add `units`

  if (length(step) != 2)
    step <- rep(step, 2)[1:2]

  topology <- attr(x, "topology")

  X <- as.matrix(cbind(x$x2d, x$y2d, x$values))
  meshes <- list()
  for (i in 1:2) {
    stopifnot(diff(range(X[, i])) <= topology[i])

    if (is.finite(topology[i])) {
      count <- round(topology[i] / step[i])
      # Fixup step
      step[i] <- topology[i] / count
      meshes[[i]] <- seq(from = min(X[, i]),
                         by = step[i],
                         length.out = count)

      # Recycle
      # TODO Move it to interpolation procedure
      X.up <- X.down <- X
      X.up[, i] <- X.up[, i] + topology[i]
      X.down[, i] <- X.down[, i] - topology[i]
      X <- rbind(X, X.up, X.down)
    } else {
      meshes[[i]] <- seq(from = min(X[, i]) + step[i] / 2,
                         to = max(X[, i]) - step[i] / 2,
                         by = step[i])
    }
  }

  grid <- interp.with.NAs(X[, 2], X[, 1], X[, 3],
                          xo = rev(meshes[[2]]),
                          yo = meshes[[1]])
  grid[c("x", "y")] <- grid[c("y", "x")]

  attr(grid, "circular") <- rev(is.finite(topology))
  attr(grid, "step") <- step

  x$field <- grid
  invisible(x)
}

deinterpolate <-  function(x, newdata) {
  x$field$z <- newdata
  x$values <- interp2(x$field$x, rev(x$field$y), x$field$z[rev(seq_len(nrow(x$field$z))), ], x$x2d, x$y2d)

  x
}

find.center.cylinder <- function(X, ...) {
  colMeans(X)
}

rotate.cylinder <- function(X, center = find.center.cylinder(X)) {
  X <- scale(X, center = center, scale = FALSE)

  U <- svd(X, nv = 3)$v

  U[, 2:3] <- cbind(c(0, 1, 0), c(0, 0, 1))
  U <- qr.Q(qr(U))

  if (U[1, 1] < 0) U[, 1] <- -U[, 1]
  if (U[2, 2] > 0) U[, 2] <- -U[, 2]
  if (det(U) < 0) U[, 3] <- -U[, 3]

  X %*% U
}

sweep.cylinder <- function(X) {
  X <- rotate.cylinder(X)
  r <- sqrt(max(X[, 2]^2 + X[, 3]^2))

  out <- cbind(X[, 1], r * atan2(X[, 2], X[, 3]))
  attr(out, "r") <- r
  attr(out, "xlim") <- range(X[, 1])

  invisible(out)
}

sweep.function <- function(data, ...) {
  X <- cbind(data$x3d, data$y3d, data$z3d)

  cyl <- sweep.cylinder(X)
  dy <- attr(cyl, "r") * 2 * pi

  res <- list(x2d = cyl[, 1],
              y2d = cyl[, 2],
              values = data$values)

  attr(res, "topology") <- c(Inf, dy)

  class(res) <- c("embryo2d", class(res))
  res
}

desweep <- function(emb3, emb2) {
  # TODO Implement this for multiple covers
  #  bad <- !is.finite(values) | !is.finite(weights) | (weights == 0)

  # res <- values
  # res[bad] <- NA

  emb3$values <- emb2$values

  emb3
}

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

  Rt <- 1 / X[, 3]
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

  R.inner <- .radius.sphere(X, side = "inner")
  R.outer <- .radius.sphere(X, side = "outer")

  depth <- (R - R.inner) / (R.outer - R.inner)
  depth[depth > 1] <- 1
  depth[depth < 0] <- 0

  x <- X[, 1] / R
  y <- X[, 2] / R
  z <- X[, 3] / R

  R.median <- median(R)
  dR.median <- median(R.outer - R.inner, na.rm = TRUE)
  x.size <- diff(asin(range(x))) * R.median
  y.size <- diff(asin(range(y))) * R.median

  N <- nrow(X)
  dx <- x.size; dy <- y.size; dd <- dR.median
  vol <- volume_ashape3d(ashape3d(X, alpha = 250)) # TODO Use proper volume estimation
  rvol <- dx * dy * dd
  dense <- N / vol
  xN <- dense^(1/3) * dx
  yN <- dense^(1/3) * dy
  dN <- dense^(1/3) * dd
  print(sprintf("dx = %f, dy = %f, dd = %f, N = %d, dense = %f, dN = %f, vol = %f, rvol = %f, R = %f",
                dx, dy, dd, N, dense, dN, vol, rvol, R.median))

  x <- x / (z + 1)
  y <- y / (z + 1)

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

.unfold3d.sphere.cylinder <- function(X) {
  # X <- rotate.sphere(X) must be done

  R <- sqrt(X[, 1]^2 + X[, 2]^2 + X[, 3]^2)

  R.inner <- .radius.sphere(X, side = "inner")
  R.outer <- .radius.sphere(X, side = "outer")

  depth <- (R - R.inner) / (R.outer - R.inner)
  depth[depth > 1] <- 1
  depth[depth < 0] <- 0

  x <- X[, 1] / R
  y <- X[, 2] / R
  z <- X[, 3] / R

  R.median <- median(R)
  dR.median <- median(R.outer - R.inner, na.rm = TRUE)
  x.size <- diff(asin(range(x))) * R.median
  y.size <- diff(asin(range(y))) * R.median

  N <- nrow(X)
  dx <- x.size; dy <- y.size; dd <- dR.median
  vol <- volume_ashape3d(ashape3d(X, alpha = 75)) # TODO Use proper volume estimation
  rvol <- dx * dy * dd
  dense <- N / vol
  xN <- dense^(1/3) * dx
  yN <- dense^(1/3) * dy
  dN <- dense^(1/3) * dd
  print(sprintf("dx = %f, dy = %f, dd = %f, N = %d, dense = %f, dN = %f, vol = %f, rvol = %f, R = %f",
                dx, dy, dd, N, dense, dN, vol, rvol, R.median))


  phi <- atan2(y, x)
  r <- sqrt(x^2 + y^2)
  x <- atan2(z, r) # Equidistant projection
  # May be use proper conformal mercator projection
  # see: http://mathworld.wolfram.com/MercatorProjection.html

  cbind(x = x, depth = depth, phi = phi) # phi is 2pi-periodic, x is NOT normalized
}

interpolate2grid.embryo3d.cylinder <- function(x, ...,
                                               cuts = c(x = 200, depth = 10, phi = 200),
                                               na.impute = TRUE,
                                               circular = FALSE) {
  stopifnot(.is.unfolded(x))

  uX <- cbind(x = x$x, depth = x$depth, phi = x$phi)
  v <- x$values

  # Omit NAs FIXME
  mask <- !is.na(rowSums(uX))

  if (na.impute) {
    # Omit NAs in values
    mask <- mask & !is.na(v)
  }

  uX <- uX[mask,, drop = FALSE]
  v <- v[mask]

  eps <- 1e-5
  ox <- seq(min(uX[, "x"]) + eps, max(uX[, "x"]) - eps, length.out = cuts["x"])
  odepth <- seq(min(uX[, "depth"]) + eps, max(uX[, "depth"]) - eps, length.out = cuts["depth"])

  if (circular) {
    ophi <- seq(0, 2*pi, length.out = cuts["phi"] + 1); ophi <- ophi[-length(ophi)]

    uX.up <- uX.down <- uX
    uX.up[, "phi"] <- uX[, "phi"] + 2*pi
    uX.down[, "phi"] <- uX[, "phi"] - 2*pi
    uX <- rbind(uX.down, uX, uX.up)
    v <- rep(v, 3)
  } else {
    phi <- uX[, "phi"]
    # FIXME FIXME FIXME fix backward (trilinear) interpolation
    zphi <- atan2(mean(sin(phi)), mean(cos(phi)))
    uX[, "phi"] <- (phi - zphi + 3*pi) %% (2*pi) - pi + zphi

    ophi <- seq(min(uX[, "phi"]) + eps, max(uX[, "phi"]) - eps, length.out = cuts["phi"])
  }

  grid <- as.matrix(expand.grid(x = ox, depth = odepth, phi = ophi))
  f <- linear.interpolate(grid, uX, v, scale = c(1, 200, 1))

  dim(f) <- sapply(list(ox, odepth, ophi), length)
  field <- list(x = ox, depth = odepth, phi = ophi, f = f)

  attr(field, "circular") <- circular

  x$field <- field

  x
}

interpolate2grid.embryo3d.sphere <- function(x, ...,
                                             cuts = c(x = 200, y = 200, depth = 10),
                                             na.impute = TRUE) {
  stopifnot(.is.unfolded(x))

  uX <- cbind(x = x$x, y = x$y, depth = x$depth)
  v <- x$values

  # Omit NAs
  mask <- !is.na(rowSums(uX))

  if (na.impute) {
    # Omit NAs in values
    mask <- mask & !is.na(v)
  }

  uX <- uX[mask,, drop = FALSE]
  v <- v[mask]

  eps <- 1e-5
  ox <- seq(min(uX[, 1]) + eps, max(uX[, 1]) - eps, length.out = cuts["x"])
  oy <- seq(min(uX[, 2]) + eps, max(uX[, 2]) - eps, length.out = cuts["y"])
  odepth <- seq(min(uX[, 3]) + eps, max(uX[, 3]) - eps, length.out = cuts["depth"])

  grid <- as.matrix(expand.grid(x = ox, y = oy, depth = odepth))

  f <- linear.interpolate(grid, uX, v)

  dim(f) <- sapply(list(ox, oy, odepth), length)
  field <- list(x = ox, y = oy, depth = odepth, f = f)

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


shrink <- function(what, to, eps = 1e-5) {
  to <- range(to)
  what[what <= to[1]] <- to[1] + eps
  what[what >= to[2]] <- to[2] - eps

  what
}

update.field.embryo3d.sphere <- function(x, newvalues = x$field$f, ...,
                                         restore.na = TRUE) {
  stopifnot(.is.interpolated(x))

  x$field$f[] <- as.numeric(newvalues)
  ox <- x$x; ox <- shrink(ox, x$field$x)
  oy <- x$y; oy <- shrink(oy, x$field$y)
  odepth <- x$depth; odepth <- shrink(odepth, x$field$depth)

  na.mask <- is.na(x$values)
  x$values <- approx3d(x$field$x, x$field$y, x$field$depth, x$field$f,
                       ox, oy, odepth)

  if (!restore.na) {
    x$values[na.mask] <- NA
  }

  x
}

approx3d.cycled <- function(x, y, z, f, xout, yout, zout) {
  z <- c(z - 2*pi, z, z + 2*pi)
  f <- rep(f, 3); dim(f) <- sapply(list(x, y, z), length)

  zout <- zout - 2*pi * floor(zout / (2*pi))

  approx3d(x, y, z, f, xout, yout, zout)
}

update.field.embryo3d.cylinder <- function(x, newvalues = x$field$f, ...,
                                           restore.na = TRUE) {
  stopifnot(.is.interpolated(x))

  x$field$f[] <- as.numeric(newvalues)

  ox <- x$x; ox <- shrink(ox, x$field$x)
  odepth <- x$depth; odepth <- shrink(odepth, x$field$depth)
  ophi <- x$phi; ophi <- shrink(ophi, x$field$phi)

  na.mask <- is.na(x$values)
  if (attr(x$field, "circular")) {
    x$values <- approx3d.cycled(x$field$x, x$field$depth, x$field$phi, x$field$f,
                                ox, odepth, ophi)
  } else {
    x$values <- approx3d(x$field$x, x$field$depth, x$field$phi, x$field$f,
                         ox, odepth, ophi)
  }

  if (!restore.na) {
    x$values[na.mask] <- NA
  }

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

  if (U[3, 3] < 0)
    U[, 3] <- -U[, 3]

  if (U[1, 1] < 0)
    U[, 1] <- -U[, 1]

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

unfold.embryo3d.sphere.cylynder <- function(x, ...) {
  X <- cbind(x$x3d, x$y3d, x$z3d)

  X <- rotate.sphere(X)

  uX <- .unfold3d.sphere.cylinder(X)
  x$x <- uX[, "x"]
  x$phi <- uX[, "phi"]
  x$depth <- uX[, "depth"]

  class(x) <- c("embryo3d.sphere.cylinder", "embryo3d.cylinder", "embryo3d")
  x
}
