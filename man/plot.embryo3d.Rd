\name{plot.embryo3d}
\alias{plot.embryo3d}
\title{Plot `embryo3d' object}

\description{
  \code{plot} method for class `embryo3d'
}

\usage{
\method{plot}{embryo3d}(x, type = c("hull", "nuclei"), \dots)
}

\arguments{
  \item{x}{object of class `embryo3d'}
  \item{type}{plot type}
  \item{\dots}{further arguments passed to plot functions}
}

\details{
  For \code{type = 'hull'}, the convex hull of nuclei is depicted,
  the color is determined by intensities. At edges the color is obtained by linear interpolation.

  For \code{type = 'nuclei'} the nuclei are deposted as color spheres.
  Color is determined by intensities.
  \strong{Warning:} Drawing can take a long time if the number of nuclei is large.

  Parameter \code{col} is vector of two colors corresponding to minimal and maximal intensivities.
  Default value is \code{c('blue', 'red')}.

  Other parameters are passed to functions
  \code{\link[rgl:shade3d]{shade3d}} and \code{\link[rgl:plot3d]{plot3d}}
  correspondingly.
}

\examples{
file <- system.file("extdata/pce", "v5_s12166-23au07-22.pce", package = "BioSSA")
pce <- read.pce.gene(file, gene = "eve")
plot(pce, type = "hull")
plot(pce, type = "nuclei")
}
