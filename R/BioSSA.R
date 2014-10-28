BioSSA <- function(x, ...)
  UseMethod("BioSSA")

BioSSA2d3d <- function(x, ...)
  UseMethod("BioSSA2d3d")

BioSSA2.5d <- BioSSA2d3d

BioSSA2d <- function(x, ...)
  UseMethod("BioSSA2d")

BioSSA3d <- function(x, ...)
  UseMethod("BioSSA3d")

BioSSA2d.embryo2d <- function(x,
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

BioSSA2d3d.embryo3d <- function(x, ..., sweep = sweep.function) {
  emb3 <- x
  emb2 <- sweep(emb3)

  bssa2d <- BioSSA(emb2, ...)

  res <- list(emb3 = emb3,
              bssa2d = bssa2d)

  class(res) <- "BioSSA2d3d"
  res
}

BioSSA.embryo2d <- BioSSA2d.embryo2d
BioSSA.embryo3d <- BioSSA2d3d.embryo3d

decompose.BioSSA2d <- function(x, ...) {
  decompose(x$ssa, ...)

  x
}

decompose.BioSSA2d3d <- function(x, ...) {
  decompose(x$bssa2d, ...)

  x
}

BioSSA.formula <- function(x, data = NULL, ...) {
  dim <- length(all.vars(x, unique = FALSE)) - 1

  if (dim == 2) {
    BioSSA2d(embryo2d(x, data = data), ...)
  } else if (dim == 3) {
    BioSSA2d3d(embryo3d(x, data = data), ...)
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

  residuals <- deinterpolate(x$emb2, attr(rec, "residuals"))
  residuals$values <- x$emb2$values
  for (j in seq_along(res)) {
    residuals$values <- residuals$values - res[[j]]$values
  }
  attr(res, "residuals") <- residuals

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSA2d.reconstruction"
  invisible(res)
}

reconstruct.BioSSA2d3d <- function(x, groups, ...) {
  rec <- reconstruct(x$bssa2d, groups = groups, ...)
  res <- lapply(rec,
                function(component) {
                  desweep(x$emb3, component)
                })
  attr(res, "series") <- x$emb3
  attr(res, "residuals") <- desweep(x$emb3, attr(rec, "residuals"))

  names(res) <- names(rec)
  attr(res, "rec") <- rec

  class(res) <- "BioSSA2d3d.reconstruction"
  invisible(res)
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

residuals.BioSSA2d3d.reconstruction <- residuals.BioSSA2d.reconstruction

residuals.BioSSA2d <- function(object, groups,
                               model = "additive", offset = 0,
                               ...) {

  rec <- reconstruct(object, groups = groups, ...)
  residuals(rec, model = model, offset = offset)
}

residuals.BioSSA2d3d <- residuals.BioSSA2d

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
