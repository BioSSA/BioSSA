# TODO Make one intepolate function

load.python <- function() {
  stopifnot(require(rPython))

  py.file <- system.file("extdata/python", "interp.py", package = "BioSSA")
  python.load(py.file)

  return(TRUE)
}

linear.interpolate.old <- function(x, points, values, scale = FALSE) {
  stopifnot(load.python())

  points <- as.matrix(points)
  points <- base::scale(points, center = FALSE, scale = scale)
  if (!is.null(attr(points, "scaled:scale"))) {
    scale <- attr(points, "scaled:scale")
  }
  x <- as.matrix(x)
  x <- base::scale(x, center = FALSE, scale = scale)
  values <- as.vector(values)

  if (length(values) < nrow(points)) {
    # Recycling
    values <- values + numeric(nrow(points))
  }
  stopifnot(length(values) == nrow(points))

  d <- ncol(points)

  x <- as.double(x)
  points <- as.double(points)
  values <- as.double(values)
  storage.mode(x) <- storage.mode(points) <- storage.mode(values) <- "double"
  storage.mode(d) <- "integer"

  res <- python.call("interpolaten", points, values, x, d)
  res <- lapply(res, function(x) if (is.null(x)) NA_real_ else as.double(x))
  res <- unlist(res)
  res <- as.double(res)

  res[is.nan(res)] <- NA_real_

  res
}

linear.interpolate <- function(...) {
  cgal.result <- LinearInterpolator:::linear.interpolate(...)
  print("Interpolated")
  # python.result <- linear.interpolate.old(...)
  #
  # eq.res <- all.equal(cgal.result, python.result)
  # print(eq.res)
  # mask <- is.na(cgal.result) | is.na(python.result)
  #
  # print(all.equal(cgal.result[!mask], python.result[!mask]))
  #
  cgal.result
}

interp.with.NAs <- function(x, y, z, ..., na.process = c("NA", "omit")) {
  mask <- is.na(z)
  if (!any(mask)) return(interp(x, y, z, ...))

  na.process <- match.arg(na.process)

  if (identical(na.process, "NA")) {
    M <- 10 * sum(abs(z[!mask]))
    zup <- zdown <- z
    zup[mask] <- M
    zdown[mask] <- -M
    iup <- interp(x, y, zup, ...)
    idown <- interp(x, y, zdown, ...)

    rmask <- iup$z != idown$z
    iup$z[rmask] <- NA

    return(iup)
  } else if(identical(na.process, "omit")) {
    x <- x[!mask]
    y <- y[!mask]
    z <- z[!mask]

    return(interp(x, y, z, ...))
  }
}
