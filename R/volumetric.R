BioSSA3d <- function(x, ...) {
  UseMethod("BioSSA3d")
}

interpolate2grid <- function(x, ...)
  UseMethod("interpolate2grid")

BioSSA3d.formula <- function(x, data = NULL, ...,
                            cuts = c(x = 100, y = 100, phi = 100, depth = 10),
                            kind = c("sphere", "cylinder", "sphere.cylinder")) {
  kind <- match.arg(kind)

  emb3 <- embryo3d(x, data = data)
  emb3 <- switch(kind,
                 sphere = unfold.embryo3d.sphere(emb3),
                 cylinder = unfold.embryo3d.cylynder(emb3),
                 sphere.cylinder = unfold.embryo3d.sphere.cylynder(emb3))

  emb3 <- interpolate2grid(emb3, cuts = cuts)

  BioSSA3d(emb3, ...)
}

BioSSA3d.embryo3d <- function(x, L, ...) {
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
  class(res) <- "BioSSA3d"

  res
}

decompose.BioSSA3d <- function(x, ...) {
  x$ssa <- decompose(x$ssa, ...)

  x
}

reconstruct.BioSSA3d <- function(x, groups, ...) {
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

  class(res) <- "BioSSA3d.reconstruction"
  invisible(res)
}


plot.BioSSA3d <- plot.BioSSA2d





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

  slice.idx <- slice.idx[order(match(names(slice.idx), names(field)))] # KILLMEPLS
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
