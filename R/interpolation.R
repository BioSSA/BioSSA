load.python <- function() {
  stopifnot(require(rPython))

  py.file <- system.file("extdata/python", "interp.py", package = "BioSSA")
  python.load(py.file)

  return(TRUE)
}

linear.interpolate <- function(x, points, values) {
  stopifnot(load.python())

  x <- as.matrix(x)
  points <- as.matrix(points)
  values <- as.vector(values)

  if (length(values) < nrow(points)) {
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
