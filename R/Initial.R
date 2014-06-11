interpolate <- function(x, step) {
  if (any(is.na(x$values))) {
    mask <- !is.na(x$values)
    xup <- xdown <- x
    M <- 10 * sum(abs(x$values[mask]))
    xup$values[!mask] <- M
    xdown$values[!mask] <- -M
    iup <- interpolate(xup, step)
    idown <- interpolate(xdown, step)

    rmask <- iup$field$z == idown$field$z
    x$field <- iup$field
    x$field$z[!rmask] <- NA

    return(invisible(x))
  }

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

  grid <- interp(X[, 2], X[, 1], X[, 3],
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

rotate.ellipsoid <- function(X) {
  X <- scale(X, center = TRUE, scale = FALSE)
  # TODO We should use svd() here
  U <- eigen(crossprod(X), symmetric = TRUE)$vectors

  # TODO Identify head and tail, up and down (m.b. copypaste?)
  if (det(U) < 0) U[, 2] <- -U[, 2]

  invisible(X %*% U)
}

sweep.cylinder <- function(X) {
  X <- rotate.ellipsoid(X)
  r <- sqrt(max(X[, 2]^2 + X[, 3]^2))

  out <- cbind(X[, 1], r * atan2(X[, 3], X[, 2]))
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

BioSSA <- function(x, ...)
  UseMethod("BioSSA")

BioSSA.embryo2d <- function(x,
                            ...,
                            step = 0.5,
                            L = c(50, 50),
                            xlim = c(-Inf, Inf),
                            ylim = c(-Inf, Inf),
                            xperc = c(0, 100),
                            yperc = c(0, 100),
                            units = c("percent", "original")) {
  emb2 <- x
  units <- match.arg(units)


  xlims.real <- range(emb2$x2d)
  ylims.real <- range(emb2$y2d)
  x1perc <- diff(xlims.real) / diff(xperc)
  x0perc <- xlims.real[1] - x1perc * xperc[1]
  y1perc <- diff(ylims.real) / diff(yperc)
  y0perc <- ylims.real[1] - y1perc * yperc[1]


  emb2[c("x0perc", "x1perc", "y0perc", "y1perc")] <- list(x0perc, x1perc, y0perc, y1perc)

  if (identical(units, "percent")) {
    if (xlim[1] < 0) xlim[1] = 0
    if (xlim[2] > 100) xlim[2] = 100
    if (ylim[1] < 0) ylim[1] = 0
    if (ylim[2] > 100) ylim[2] = 100

    stopifnot(all(xlim <= max(xperc)),
              all(xlim >= min(xperc)),
              all(ylim <= max(yperc)),
              all(ylim >= min(yperc)))

    xlim <- xlim * x1perc + x0perc
    ylim <- ylim * y1perc + y0perc
    step <- rep(step, 2)[1:2]
    step.x <- step[1] * x1perc
    step.y <- step[2] * y1perc
    step <- c(step.x, step.y)
  }

  emb2 <- interpolate(emb2,
                      step = step)

  ### Cutoff
  xmask <- (emb2$field$x >= xlim[1]) & (emb2$field$x <= xlim[2])
  ymask <- (emb2$field$y >= ylim[1]) & (emb2$field$y <= ylim[2])
  mask <- outer(ymask, xmask, "&") #TODO MB, speedup?

  if (identical(units, "percent")) {
    L <- rep(L, 2)[1:2]
    L.x <- round(L[1] / 100 * length(emb2$field$x))
    L.y <- round(L[2] / 100 * length(emb2$field$y))
    L <- c(L.y, L.x)
  }

  dec <- ssa(emb2$field$z,
             L = L,
             ...,
             kind = "2d-ssa",
             mask = mask,
             circular = attr(emb2$field, "circular"))

  res <- list(emb2 = emb2,
              ssa = dec)
  class(res) <- "BioSSA2d"

  res
}

BioSSA.embryo3d <- function(x, ..., sweep = sweep.function) {
  emb3 <- x
  emb2 <- sweep(emb3)

  bssa2d <- BioSSA(emb2, ...)

  res <- list(emb3 = emb3,
              bssa2d = bssa2d)

  class(res) <- "BioSSA3d"
  res
}

decompose.BioSSA2d <- function(x, ...) {
  decompose(x$ssa, ...)

  x
}

decompose.BioSSA3d <- function(x, ...) {
  decompose(x$bssa2d, ...)

  x
}

BioSSA.formula <- function(x, data = NULL, ...) {
  dim <- length(all.vars(x, unique = FALSE)) - 1

  if (dim == 2) {
    BioSSA(embryo2d(x, data = data), ...)
  } else if (dim == 3) {
    BioSSA(embryo3d(x, data = data), ...)
  } else {
    stop("Incorrect `formula' argument dimension")
  }
}

reconstruct.BioSSA2d <- function(x, groups, ...) {
  rec <- reconstruct(x$ssa, groups = groups, ...)
  res <- lapply(rec,
                function(component) {
                  deinterpolate(x$emb2, component)
                })
  attr(res, "series") <- x$emb2
  attr(res, "residuals") <- deinterpolate(x$emb2, attr(rec, "residuals"))

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSA2d.reconstruction"
  invisible(res)
}

reconstruct.BioSSA3d <- function(x, groups, ...) {
  rec <- reconstruct(x$bssa2d, groups = groups, ...)
  res <- lapply(rec,
                function(component) {
                  desweep(x$emb3, component)
                })
  attr(res, "series") <- x$emb3
  attr(res, "residuals") <- desweep(x$emb3, attr(rec, "residuals"))

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSA3d.reconstruction"
  invisible(res)
}

embryo2d <- function(x, ...)
  UseMethod("embryo2d")

embryo2d.formula <- function(x,
                             data = NULL,
                             topology = c(Inf, Inf),
                             ...) {
  mf <- model.frame(x, data = data, ...)

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
  mf <- model.frame(x, data = data, ...)

  emb3 <- list(x3d = mf[, 2],
               y3d = mf[, 3],
               z3d = mf[, 4],
               values = mf[, 1])

  class(emb3) <- "embryo3d"
  emb3
}

residuals.BioSSA2d.reconstruction <- function(object,
                                              model = "additive", offset = 0,
                                              ...) {
  rec <- object # Fuck the system!!!

  res <- attr(object, "residuals")
  series <- attr(rec, "series")

  if (is.character(model)) {
    model <- match.arg(model, c("additive", "multiplicative", "poisson"))

    alpha <- switch(model,
                    additive = 0,
                    multiplicative = 1,
                    poisson = 0.5)
  } else {
    alpha <- model
  }

  res$values <- res$values / (series$values - res$values + offset) ^ alpha

  res
}

residuals.BioSSA3d.reconstruction <- residuals.BioSSA2d.reconstruction

residuals.BioSSA2d <- function(object, groups,
                               model = "additive", offset = 0,
                               ...) {

  rec <- reconstruct(object, groups = groups, ...)
  residuals(rec, model = model, offset = offset)
}

residuals.BioSSA3d <- residuals.BioSSA2d

wcor.BioSSA2d <- function(x, ...) {
  wcor(x$ssa, ...)
}

read.emb.data <- function(file) {
  lines <- readLines(file)

  pattern <- "channel \\d - ([a-z]*)"
  re <- gregexpr(pattern, lines[5])
  matches <- regmatches(lines[5], re)[[1]]
  extract <- function(m) {
    re <- regexec(pattern, m)
    regmatches(m, re)[[1]][2]
  }

  gene.names <- sapply(matches, extract)

  df <- read.table(file, skip = 5, header = FALSE)
  names(df) <- c("N", "AP", "DV", gene.names)

  df
}
