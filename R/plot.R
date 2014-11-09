.defaults <- function(x, ...) {
  dots <- list(...)
  modifyList(dots, x)
}

.get.colors <- function(values, indices = seq_along(values), col = c("blue", "red")) {
  values[!is.finite(values)] <- mean(values[is.finite(values)])
  grade <- (values - min(values)) / diff(range(values))
  cl <- colorRamp(col)(grade[indices]) / 256
  rgb(cl[, 1], cl[, 2], cl[, 3])
}

.plot.embryo3d.nuclei <- function(data, ...) {
  data <- as.data.frame(data[c("x3d", "y3d", "z3d", "values")])
  mask <- !is.na(data$values)
  data <- data[mask, ]

  dots <- list(...)
  dots <- .defaults(dots,
                    aspect = FALSE,
                    type = "s",
                    col = c("blue", "red"),
                    xlab = "x",
                    ylab = "y",
                    zlab = "z")
  dots$col <- .get.colors(data$values, col = dots$col)

  do.call("plot3d", c(list(data$x3d, data$y3d, data$z3d), dots))
}

.plot.embryo3d.hull <- function(data, ..., alpha = Inf) {
  data <- as.data.frame(data[c("x3d", "y3d", "z3d", "values")])
  mask <- !is.na(data$values)
  data <- data[mask, ]

  X <- cbind(data$x3d, data$y3d, data$z3d)
  ash <- ashape3d(X, alpha = alpha, pert = TRUE)
  triang <- ash$triang
  ach <- triang[triang[, 9] == 2, 1:3] # Extract regular triangles (0 --- not in, 1 --- interior, 2 --- regular, 3 --- singular)

  mesh <- tmesh3d(t(cbind(X, 0.5)), t(ach))

  dots <- list(...)
  dots <- .defaults(dots,
                    aspect = FALSE,
                    col = c("blue", "red"),
                    xlab = "x",
                    ylab = "y",
                    zlab = "z")
  dots$col <- .get.colors(data$values, t(ach), col = dots$col)

  plot3d(data$x3d, data$y3d, data$z3d, type = "n",
         aspect = dots$aspect,
         xlab = dots$xlab,
         ylab = dots$ylab,
         zlab = dots$zlab)
  do.call("shade3d", c(list(mesh), dots))
}

# TODO Merge with field.section.embryo3d
field.section2d.embryo3d <- function(emb3, slice = list(), units = c("percent", "original")) {
  stopifnot(.is.interpolated(emb3))

  units <- match.arg(units)

  field <- emb3$field

  # TODO Do smth with it
  # if (identical(units, "percent")) {
  #   field$x <- (field$x - emb2$x0perc) / emb2$x1perc
  #   field$y <- (field$y - emb2$y0perc) / emb2$y1perc
  # }

  # Determine free coords
  # FIXME reorder if needed
  free.coords <- setdiff(names(field), c("f", names(slice)))
  free.coords <- free.coords[order(match(free.coords, names(field)))]

  slice.idx <- lapply(names(slice),
                      function(name) {
                        s <- slice[[name]]
                        coord <- field[[name]]
                        name <- names(s)
                        which.min(abs(coord - s))
                      })
  names(slice.idx) <- names(slice) # FIXME MB this is not needed
  stopifnot(length(slice.idx) == 1)
  tmp <- list(TRUE, TRUE)
  names(tmp) <- free.coords
  slice.idx <- c(slice.idx, tmp)
  stopifnot(length(slice.idx) == 3)

  slice.idx <- slice.idx[order(match(names(slice.idx), names(field)))] # KILLMEPLS
  names(slice.idx) <- NULL

  values <- do.call("[", c(list(field$f, drop = TRUE), slice.idx))
  list(x = field[[free.coords[1]]],
       y = field[[free.coords[2]]],
       z = values, free.coords = free.coords)
}

.plot3d.embryo3d.section2d.field <- function(x, slice, ..., add = FALSE, col = grey(0:1),
                                             grid = c("x", "y", "z")) {
  x <- field.section2d.embryo3d(x, slice = slice)
  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = x$free.coords[1],
                    ylab = x$free.coords[2],
                    zlab = "value",
                    main = "")
  if (!add) plot3d(range(x$x), range(x$y), range(x$z, na.rm = TRUE), type = "n",
                   xlab = dots$xlab, ylab = dots$ylab, zlab = dots$zlab,
                   main = dots$main)
  do.call(surface3d,
          c(list(x$x, x$y, x$z, col = .get.colors(t(x$z), col = col)), dots))

  if (!add) {
    if (!identical(grid, FALSE))
      grid3d(grid)
  }
}

.plot3d.embryo3d.section2d.nuclei <- function(x,
                                              slice,
                                              tolerance = 0.1,
                                              units = c("percent", "original"),
                                              ...,
                                              ref = FALSE, add = FALSE,
                                              grid = c("x", "y", "z")) {
  dots <- list(...)
  units <- match.arg(units)

  stripe <- subset(x, subset = slice, tolerance = tolerance)
  free.coords <- setdiff(names(x), c(names(slice), "field", "values", "x3d", "y3d", "z3d"))
  free.coords <- free.coords[order(match(free.coords, names(x)))]
  stripe$xvalues <- stripe[[free.coords[1]]]
  stripe$yvalues <- stripe[[free.coords[2]]]

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = "s",
                    xlab = free.coords[1],
                    ylab = free.coords[2],
                    zlab = "values",
                    main = "",
                    size = 1,
                    axes = TRUE,
                    bbox = TRUE,
                    col = grey(c(0, 1)))
  plot3d(stripe$xvalues, stripe$yvalues, stripe$values,
         type = dots$type,
         main = dots$main,
         size = dots$size,
         axes = dots$axes,
         xlab = dots$xlab, ylab = dots$ylab, zlab = dots$zlab,
         col = .get.colors(stripe$values, col = dots$col),
         add = add)

  if (ref) {
    surface3d(range(stripe$xvalues),
              range(stripe$yvalues),
              matrix(0, 2, 2), col = "green", alpha = 0.25)
  }

  if (!add) {
    if (!identical(grid, FALSE))
      grid3d(grid)
  }
}

plot.embryo3d <- function(x, type = c("hull", "nuclei",
                                      "field-section", "nuclei-section",
                                      "field-section3d", "nuclei-section3d"),
                          ...) {
  type <- match.arg(type)

  if (identical(type, "hull")) {
    .plot.embryo3d.hull(x, ...)
  } else if (identical(type, "nuclei")) {
    .plot.embryo3d.nuclei(x, ...)
  } else if (identical(type, "field-section")) {
    .plot1d.embryo3d.section.field(x, ...)
  } else if (identical(type, "nuclei-section")) {
    .plot1d.embryo3d.section.nuclei(x, ...)
  } else if (identical(type, "nuclei-section3d")) {
    .plot3d.embryo3d.section2d.nuclei(x, ...)
  } else if (identical(type, "field-section3d")) {
    .plot3d.embryo3d.section2d.field(x, ...)
  } else {
    stop("Unknown `type'")
  }
}

plot.BioSSA2d.reconstruction <- function(x, ...) {
  step <- attr(attr(x, "series")$field, "step")
  plot(attr(x, "rec"), aspect = step[2] / step[1], ...)
}

panel.lmline0 <- function(x, y, ..., identifier = "lmline") {
  if (length(x) > 1)
    panel.abline(lm(as.numeric(y) ~ as.numeric(x) + 0), ...,
                 identifier = identifier)
}

.plot.res.vs.res <- function(x, groups.x, groups.y, ..., lm = TRUE) {
  groups.x <- sort(unique(unlist(groups.x)))
  groups.y <- sort(unique(unlist(groups.y)))

  dots <- list(...)
  dots <- .defaults(dots,
                    xlab = sprintf("group = (%s)", paste0(groups.x, collapse = ", ")),
                    ylab = sprintf("group = (%s)", paste0(groups.y, collapse = ", ")),
                    pch = 19)

  do.call("xyplot", c(list(residuals(x, groups = groups.y)$values ~ residuals(x, groups = groups.x)$values,
                           panel = function(...) {
                             panel.xyplot(...)

                             if (lm)
                               panel.lmline0(...)
                           }),
                    dots))
}

.plot.BioSSA2d.residuals <- function(x, groups, groups.y, ...) {
  if (missing(groups.y))
    plot(noise.model(x, groups = groups, ...), ...)
  else {
    .plot.res.vs.res(x, groups.x = groups, groups.y = groups.y, ...)
  }
}

plot.BioSSA2d <- function(x, type = c("ssa-values", "ssa-vectors", "ssa-series", "ssa-wcor",
                                      "residuals"), ...) {
  type <- match.arg(type)
  if (grepl("^ssa-", type)) {
    ssa.type <- sub("ssa-", "", type)

    plot(x$ssa, type = ssa.type, ...)
  } else if (identical(type, "residuals")) {
    .plot.BioSSA2d.residuals(x, ...)
  } else {
    stop("Unknown `type'")
  }
}

.plot2d.embryo2d.nuclei <- function(x, ..., voronoi = TRUE, col = grey(0:1)) {
  dots <- list(...)

  dots <- .defaults(dots,
                    par.settings = list(regions = list(col = colorRampPalette(col))),
                    aspect = "iso",
                    xlab = "x",
                    ylab = "y",
                    pch = 19,
                    cex = 0.1,
                    points = FALSE,
                    border = TRUE)
  if (voronoi) {
    do.call("levelplot", c(list(values ~ x2d * y2d,
                                data = x,
                                panel = panel.voronoi,
                                use.tripack = TRUE),
                           dots))
  } else {
    do.call("xyplot", c(list(y2d ~ x2d, data = x), dots))
  }
}

.plot2d.embryo2d.field <- function(x, ...,
                                   col = grey(0:1),
                                   cuts = 20,
                                   zlim,
                                   at,
                                   useRaster = TRUE) {
  dots <- list(...)

  dots <- .defaults(dots,
                    par.settings = list(regions = list(col = colorRampPalette(col))),
                    aspect = "iso",
                    xlab = "x",
                    ylab = "y")

  data <- expand.grid(x = x$field$x, y = x$field$y)
  data$z <- as.vector(as.numeric(t(x$field$z)))

  if (missing(zlim)) {
    zlim <- data$z
  }

  if (missing(at)) {
    at <- pretty(zlim, n = cuts)
  }

  # Cutoff outstanding values
  data$z[data$z < min(at)] <- min(at)
  data$z[data$z > max(at)] <- max(at)

  do.call("levelplot", c(list(z ~ x * y,
                              data = data,
                              at = at,
                              useRaster = useRaster),
                         dots))
}

.plot3d.embryo2d.field <- function(x, ..., add = FALSE) {
  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    xlab = "x",
                    ylab = "y",
                    zlab = "z",
                    main = "",
                    col = grey(c(0, 1)))
  if (!add) plot3d(x$x2d, x$y2d, x$values, type = "n",
                   xlab = dots$xlab, ylab = dots$ylab, zlab = dots$zlab,
                   main = dots$main)
  surface3d(x$field$x, x$field$y, t(x$field$z),
            col = .get.colors(t(x$field$z), col = dots$col))
  grid3d(c("x", "y", "z"))
}

.plot3d.embryo2d.nuclei <- function(x, ..., add = FALSE) {
  dots <- list(...)

  # Provide convenient defaults
  dots <- .defaults(dots,
                    type = "s",
                    xlab = "x",
                    ylab = "y",
                    zlab = "z",
                    main = "",
                    size = 1,
                    axes = TRUE,
                    bbox = TRUE,
                    col = grey(c(0, 1)))
  plot3d(x$x2d, x$y2d, x$values,
         type = dots$type,
         main = dots$main,
         size = dots$size,
         axes = dots$axes,
         xlab = dots$xlab, ylab = dots$ylab, zlab = dots$zlab,
         col = .get.colors(x$values, col = dots$col),
         add = add)

  grid3d(c("x", "y", "z"))
}

field.section <- function(emb2, at, coord = "y", units = c("percent", "original")) {
  units <- match.arg(units)
  if (inherits(emb2, "embryo3d")) emb2 <- embryo2d(emb2)
  field <- emb2$field

  if (identical(units, "percent")) {
    field$x <- (field$x - emb2$x0perc) / emb2$x1perc
    field$y <- (field$y - emb2$y0perc) / emb2$y1perc
  }

  axe <- if (coord == "y") field$y else field$x
  caxe <- if (coord == "x") field$y else field$x
  idx <- which.min(abs(axe - at) %% diff(range(axe)))

  values <- if (coord == "y") field$z[idx, ] else field$z[, idx]

  list(scale = caxe, values = values)
  list(x = caxe, y = values)
}

nuclei.stripe <- function(emb2, at, units = c("percent", "original"), tolerance = 0.1, coord = "y") {
  units <- match.arg(units)

  if (inherits(emb2, "embryo3d")) emb2 <- embryo2d(emb2)

  if (identical(units, "percent")) {
    emb2$x2d <- (emb2$x2d - emb2$x0perc) / emb2$x1perc
    emb2$topology[1] <- (emb2$topology[1] - emb2$x0perc) / emb2$x1perc
    emb2$y2d <- (emb2$y2d - emb2$y0perc) / emb2$y1perc
    emb2$topology[2] <- (emb2$topology[2] - emb2$y0perc) / emb2$y1perc
  }

  axe <- emb2[[paste0(coord, "2d")]]

  along.coord.index <- if (coord == "x") 1 else 2
  topology <- emb2$topology[along.coord.index]
  idx <- (abs(axe - at) %% min(topology, 10e10)) < tolerance # TODO Get rid off this

  out <- list(x2d = emb2$x2d[idx],
              y2d = emb2$y2d[idx],
              values = emb2$values[idx])

  attr(out, "topology") <- attr(emb2, "topology")
  attr(out, "topology")[-along.coord.index] <- FALSE

  class(out) <- c("embryo2d", class(out))
  out
}

.plot1d.embryo2d.section.field <- function(x,
                                           at = 50,
                                           coord = c("y", "x"),
                                           units = c("percent", "original"),
                                           ...,
                                           ref = FALSE) {
  dots <- list(...)
  coord <- match.arg(coord)
  units <- match.arg(units)

  section <- field.section(x, at, coord = coord, units = units)
  anticoord <- if (coord == "y") "x" else "y"

  dots <- .defaults(dots,
                    type = "l",
                    col = "black",
                    ylab = "Expression",
                    xlab = sprintf("Spatial coordinate: %s%s", anticoord, if (identical(units, "percent")) ", %" else ""))

  res <- do.call("xyplot", c(list(y ~ x, data = section), dots))

  if (ref) {
    res <- res + layer(panel.abline(h = 0, reference = TRUE))
  }

  res
}

.plot1d.embryo2d.section.nuclei <- function(x,
                                            at = 50,
                                            coord = c("y", "x"),
                                            tolerance = 0.1,
                                            units = c("percent", "original"),
                                            ...,
                                            ref = FALSE) {
  dots <- list(...)
  coord <- match.arg(coord)
  units <- match.arg(units)

  stripe <- nuclei.stripe(x, at, coord = coord, tolerance = tolerance, units = units)
  anticoord <- if (coord == "y") "x" else "y"

  dots <- .defaults(dots,
                    type = "p",
                    col = "green",
                    pch = 18,
                    ylab = "Expression",
                    xlab = sprintf("Spatial coordinate: %s%s", anticoord, if (identical(units, "percent")) ", %" else ""))

  stripe$xvalues <- stripe[[if (coord == "x") "y2d" else "x2d"]]
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

rgl.move <- function(x = c(1, 0, 0),
                     y = c(0, 1, 0),
                     z = c(0, 0, 1),
                     offset = c(0, 0, 0)) {
  tensor <- rbind(cbind(x, y, z, offset), c(0, 0, 0, 1))
  par3d(userMatrix = tensor)
}
