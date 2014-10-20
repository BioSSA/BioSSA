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

  Rt <- 1 / sqrt(R^2 - X[, 1]^2 - X[, 2]^2)
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
  Rt <- 1 / sqrt(R^2 - X[, 1]^2 - X[, 2]^2)
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
  field <- list(x = ox, y = oy, oz = oz, f = f)

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
  ox <- seq(min(uX[, 1]) + eps, max(uX[, 1]) - eps, length.out = cuts[1])
  oy <- seq(min(uX[, 2]) + eps, max(uX[, 2]) - eps, length.out = cuts[2])
  oz <- seq(min(uX[, 3]) + eps, max(uX[, 3]) - eps, length.out = cuts[3])

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
  z <- c(z - 2*pi, z, 2*pi)
  f <- rep(f, 3); dim(f) <- sapply(list(x, y, z), length)

  zout <- zout - 2*pi * floor(zout / (2*pi))

  approx3d(x, y, z, f, xout, yout, zout)
}

update.field.embryo3d.cylynder <- function(x, newvalues = x$field$f, ...) {
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

BioSSAv.formula <- function(x, data = NULL, ..., kind = c("sphere", "cylinder")) {
  kind <- match.arg(kind)

  emb3 <- embryo3d(x, data = data)
  emb3 <- switch(kind,
                 sphere = unfold.embryo3d.sphere(emb3),
                 cylinder = unfold.embryo3d.cylynder(emb3))

  emb3 <- interpolate2grid(emb3)

  BioSSAv(emb3, ...)
}

BioSSAv.embryo3d <- function(x, ...) {
  stopifnot(.is.interpolated(x))
  f <- x$field$f

  if (inherits(x, "embryo3d.sphere")) {
    circular <- c(FALSE, FALSE, FALSE)
  } else if (inherits(x, "embryo3d.cylinder")) {
    circular <- c(FALSE, FALSE, TRUE)
  }

  dec <- ssa(f, ..., kind = "nd-ssa", circular = circular)

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
  attr(res, "residuals") <- update.field(x$emb3, attr(rec, "residuals"))

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSAv.reconstruction"
  invisible(res)
}


residuals.BioSSAv <- residuals.BioSSA2d
residuals.BioSSAv.reconstruction <- residuals.BioSSA2d.reconstruction
plot.BioSSAv <- plot.BioSSA2d
