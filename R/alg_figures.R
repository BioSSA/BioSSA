rgl.kill.text <- function() {
  ids <- rgl.ids()

  rgl.pop(id = ids$id[ids$type == "text"])
}


save.rgl.plot <- function(filename, make.png = TRUE, make.vec = TRUE, fmts = c("eps"), do.pause = TRUE) {
  figs <- readLines("fig_list.txt")
  if (!filename %in% figs && !grepl("bbox", filename)) return(0)

  print(sprintf("3d Plotting %s .......", filename))

  # complicated <- filename %in% readLines("complicated.txt")
  complicated <- FALSE


  if (grepl("profile", filename) && !("svg" %in% fmts)) {
    fmts <- c(fmts, "svg")
  }

  if (complicated) make.vec <- FALSE

  if (do.pause) Sys.sleep(1.5)

  if (make.png) {
    rgl.snapshot(filename)
    system(sprintf("convert %s -transparent white %s", filename, filename))
  }

  if (make.vec) {
    for (fmt in fmts) {
      # rgl.postscript(sub("\\.png", sprintf(".%s", fmt), filename), fmt = fmt)
    }
  }

  print(sprintf("3d Plotted %s", filename))
}


plot.em <- function(x, lit = TRUE, type = c("hull", "nuclei"),
                    alpha = 50,
                    add.points = 25000, cos.trashold = 0, ..., plot.empty.nuclei = FALSE,
                    plot.axes = TRUE, plot.plot = TRUE) {
  em <- x
  type <- match.arg(type)

  data <- unclass(em)[c("x3d", "y3d", "z3d", "values")]

  if (add.points > 0) {
    rot <- attr(em, "rotated")
    unf <- attr(em, "unfolded")

    R <- attr(unf, "R")
    d <- attr(unf, "d")
    center <- attr(rot, "center")
    U <- attr(rot, "U")

    set.seed(1)
    X <- matrix(rnorm(3 * add.points), ncol = 3)
    X <- X / sqrt(rowSums(X^2))
    coss <- X %*% U[, 3]
    X <- X[coss > cos.trashold,, drop = FALSE]

    Rs <- runif(nrow(X), min = R, max = R + d)
    X <- X * Rs
    X <- X + rep(center, each = nrow(X))

    data.add <- data.frame(x3d = X[, 1], y3d = X[, 2], z3d = X[, 3], values = NA)
    data.add <- data.add[names(data)]
    data <- if (plot.empty.nuclei) rbind(data, data.add) else data.add
  }

  if (plot.axes) {
    with(data, plot3d(x3d, y3d, z3d, type = "n", col = "gray", aspect = "iso",
                      xlab = "x", ylab = "y", zlab = "z"))
  } else {
    open3d()
    with(data, plot3d(x3d, y3d, z3d, type = "n", col = "gray", aspect = "iso",
                      add = TRUE))
  }

  if (!plot.plot) return(0)


  material3d(alpha = 0.25, col = "grey", lit = lit)
  with(data, points3d(x3d, y3d, z3d))

  material3d(alpha = 0.99, col = "grey", size = if (identical(type, "hull")) 0 else 1, lit = lit)
  plot(em, type = type, alpha = alpha, add = TRUE, aspect = "iso", ...)
}


.pl <- function(X, hull, values = 0, ..., add = TRUE) {
  mesh <- tmesh3d(t(X), t(hull), homogeneous = FALSE)

  data <- data.frame(x3d = X[, 1], y3d = X[, 2], z3d = X[, 3], values = values)


  dots <- list(...)
  dots <- .defaults(dots,
                    aspect = "iso",
                    col = "red",
                    xlab = "x",
                    ylab = "y",
                    override = FALSE,
                    zlab = "z")

  if (!add) {
    plot3d(data$x3d, data$y3d, data$z3d, type = "n",
           aspect = dots$aspect,
           xlab = dots$xlab,
           ylab = dots$ylab,
           zlab = dots$zlab)
  }
  do.call("shade3d", c(list(mesh), dots))
}


plot.bb <- function(x, what = c("original", "rotated", "equalized"), ...) {
  em <- x
  what <- match.arg(what)

  rot <- attr(em, "rotated")
  unf <- attr(em, "unfolded")
  eqv <- attr(unf, "X.eq")

  h.outer <- attr(unf, "outer.surface")
  h.inner <- attr(unf, "inner.surface")

  X <- switch(what,
              original = cbind(em$x3d, em$y3d, em$z3d),
              rotated = rot,
              equalized = eqv)
  plot3d(X[, 1], X[, 2], X[, 3], aspect = "iso")

  .pl(X, h.inner, values = em$values)
  .pl(X, h.outer, values = em$values)
}


.plot.embryo3d.nuclei.flattened <- function(data, ..., plot.type = "s",
                                            scale.range = range(data$values, na.rm = TRUE),
                                            plot.triangilation = FALSE,
                                            original.units = TRUE) {
  un <- attr(data, "units")

  data <- as.data.frame(data[c("psi", "phi", "depth", "values")])

  if (original.units) {
    for (unn in names(un)) {
      data[[unn]] <- data[[unn]] * un[unn]
    }
  }

  mask <- !is.na(data$values)
  data <- data[mask, ]

  dots <- list(...)
  dots <- .defaults(dots,
                    aspect = "iso",
                    type = plot.type,
                    col = c("blue", "red"),
                    xlab = "psi",
                    ylab = "phi",
                    zlab = "depth")
  dots$col <- .get.colors(data$values, col = dots$col, scale.range = scale.range)

  do.call("plot3d", c(list(data$psi, data$phi, data$depth), dots))

  if (plot.triangilation) {
    X <- as.matrix(data[, 1:3])
    ash <- ashape3d(X, alpha = Inf, pert = TRUE)
    tetra <- ash$tetra
    triang <- ash$triang
    indices <- tetra[, 1:4, drop = FALSE]
    indices <- triang[, 1:3, drop = FALSE]

    mesh <- tmesh3d(t(X), t(indices), homogeneous = FALSE)

    wire3d(mesh, col = "grey")
  }
}


field2emb <- function(em) {
  field <- em$field
  data <- expand.grid(psi = field$psi, depth = field$depth, phi = field$phi)
  data$values <- as.vector(field$f)

  attr(data, "units") <- attr(em, "units")

  data
}


make.gl <- function(filename, expr, write.webgl = FALSE, lit = TRUE, zoom = 0.725, FOV = 0, cex = 1.5) {
  windowRect <- c(500, 0, 500 + 1000, 700)

  um <- structure(c(-0.577869534492493, -0.595892071723938, 0.557655572891235,
                    0, 0.815322816371918, -0.451870083808899, 0.362023651599884,
                    0, 0.0362609140574932, 0.663872301578522, 0.74696671962738, 0,
                    0, 0, 0, 1), .Dim = c(4L, 4L))

  try(rgl.close())
  par3d(userMatrix = um, zoom = zoom, FOV = FOV, windowRect = windowRect, cex = cex)
  material3d(alpha = 0.99, col = "grey", size = 1, lit = lit)
  eval(expr, envir = parent.frame())
  par3d(userMatrix = um, zoom = zoom, FOV = FOV, windowRect = windowRect, cex = cex)

  save.rgl.plot(filename)

  if (write.webgl) {
    writeWebGL(sub("\\.png", "-webGL", filename))
  }
}
