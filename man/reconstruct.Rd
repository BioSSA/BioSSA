\name{reconstruct}
\alias{reconstruct}
\alias{reconstruct.BioSSA2d}
\alias{reconstruct.BioSSA3d}
\title{Component reconstruction for BioSSA decomposition object}

\description{
  Reconstruct the data given the BioSSA decomposition and the desired
  grouping of the elementary components.
}

\usage{
\method{reconstruct}{BioSSA2d}(x, groups, ...)
\method{reconstruct}{BioSSA3d}(x, groups, ...)
}

\arguments{
  \item{x}{BioSSA decomposition object}
  \item{groups}{list of numeric vectors, indices of elementary components
    used for reconstruction, the entries of the list can be named, see
    'Value' for more information}
  \item{\dots}{further arguments passed to corresponding SSA reconstruction method}
}

\details{
  For `embryo2d' object, the function run \code{reconstruct} method for stored
  `2d.ssa' object and then performs back interpolation to nuclei positions.

  For `embryo3d' object, the function run reconstruction for all `embryo2d'
  that are stored within `embrio3d', then perform desweeping and joining reconstructed
  fragments of whole 3d embryo.
}

\value{
  List of reconstructed embryos (`embrio2d' or `embrio3d' correspondingly).
  Elements of the list have the same
  names as elements of \code{groups}. If the group is unnamed, then
  corresponding component will obtain name `Fn', where `n' is its index
  in \code{groups} list.
}

\seealso{
  \code{\link[BioSSA:BioSSA]{BioSSA}},
  \code{\link[Rssa:reconstruct]{reconstruct}}.
}
