\name{decompose}
\alias{decompose}
\alias{decompose.BioSSA2d}
\alias{decompose.BioSSA3d}
\title{Performs the SSA decomposition of stored `2d.ssa. object(s)}

\description{
  Perform the SSA decomposition of stored `2d.ssa' object(s).
}

\usage{
\method{decompose}{BioSSA2d}(x, ...)
\method{decompose}{BioSSA3d}(x, ...)
}

\arguments{
  \item{x}{object of class `BioSSA2d' or `BioSSA3d'}
  \item{\dots}{further arguments passed to the internal \code{decompose} method call}
}

\details{
  For the `BioSSA2d' object \code{decompose} calls the \code{\link[Rssa:decompose]{decompose}}
  method from `Rssa' package applied to the `2d.ssa' object stored in input `BioSSA2d' object.

  For the `BioSSA3d` object \code{decompose} calls the \code{decompose} method applied to
  all stored `BioSSA2d' objects.
}

\seealso{
  \code{\link[Rssa:decompose]{decompose}} from package `Rssa'.
}

