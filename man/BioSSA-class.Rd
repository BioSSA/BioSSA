\name{BioSSA-class}
\alias{BioSSA-class}
\alias{BioSSA2d-class}
\alias{BioSSA3d-class}
\alias{BioSSA2d3d-class}

\title{S3 classes for BioSSA decompositions}

\description{
  Class for holding BioSSA decomposition
}

\details{
  `BioSSA2d' contains 2d data, their interpolation into a regular grid and
  the deomposition of the interpolated data (of their trajectory matrix).

  `BioSSA3d' contains 3d data and the list of its planar projections in the form of
  `BioSSA2s' class (at the moment, it consists of only one cylindrical projection).
}

\seealso{
  \code{\link[BioSSA:BioSSA]{BioSSA}}.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bss <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
}
