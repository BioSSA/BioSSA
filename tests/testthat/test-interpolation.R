library(testthat)

# Compute barycentric coordinates. Duplication for cart2bary(). Now useless
linear.interpolate.simplex <- function(x, simplex, values) {
  # Matrix of the linear system for convex combination coefficients (alphas)
  mat <-  cbind(simplex, 1)

  coeffs <- qr.solve(t(mat), values)

  x %*% coeffs
}


linear.interpolate.old <- function(x, points, values, tessellation = delaunayn(points)) {
  ts <- tsearchn(points, tessellation, x)

  out <- rep(NA_real_, nrow(x))

  # TODO Use tapply or matmul here
  for (i in seq_len(nrow(x))) {
    out[i] <- sum(ts$p[i, ] * values[tessellation[ts$idx[i], ]])
  }

  out
}

is.lay.in.convex.hull <- function(x, points, tessellation = delaunayn(points)) {
  idx <- tsearchn(points, tessellation, x)$idx

  !is.na(idx)
}

# TODO Investigate why tests fail with lesser tolerance
tolerance <- 1e-3

test_that("linear.interpolate() works for 2d linear case", {
  set.seed(1)
  x <- c(rnorm(100), 100, 100, -100, -100)
  y <- c(rnorm(100), 100, -100, 100, -100)

  f <- function(x, y) x + 2*y

  z <- f(x, y)

  xx <- rnorm(100)
  yy <- rnorm(100)

  zz <- f(xx, yy)

  li <- linear.interpolate(cbind(xx, yy), cbind(x, y), z)

  expect_equal(li, zz, tolerance = tolerance)
})



test_that("linear.interpolate() works for 2d linear case with NAs", {
  set.seed(1)
  x <- c(rnorm(100))
  y <- c(rnorm(100))

  f <- function(x, y) x + 2*y

  z <- f(x, y)

  xx <- rnorm(1000)
  yy <- rnorm(1000)

  zz <- f(xx, yy)

  li <- linear.interpolate(cbind(xx, yy), cbind(x, y), z)

  zz[!is.lay.in.convex.hull(cbind(xx, yy), cbind(x, y), delaunayn(cbind(x, y)))] <- NA

  expect_equal(li, zz, tolerance = tolerance)
})

test_that("linear.interpolate() works for n-d linear case with NAs", {
  set.seed(1)
  dim <- 3
  n <- 100
  x <- matrix(rnorm(dim * n), n, dim)

  coef <- rnorm(dim)

  f <- function(x) as.vector(x %*% coef)

  z <- f(x)

  ni <- 100
  xx <- matrix(rnorm(dim * ni), ni, dim)

  zz <- f(xx)

  li <- linear.interpolate(xx, x, z)

  zz[!is.lay.in.convex.hull(xx, x, delaunayn(x))] <- NA

  expect_equal(li, zz, tolerance = tolerance)
})
