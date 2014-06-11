\name{residuals.BioSSA}
\alias{residuals.BioSSA2d}
\alias{residuals.BioSSA3d}
\alias{residuals.BioSSA2d.reconstruction}
\alias{residuals.BioSSA3d.reconstruction}

\title{
  Residuals of BioSSA decomposition
}

\description{
  Returns BioSSA residuals in the specified model.
}

\usage{
\method{residuals}{BioSSA2d}(object, groups, model = "additive", offset = 0, \dots)
\method{residuals}{BioSSA3d}(object, groups, model = "additive", offset = 0, \dots)
\method{residuals}{BioSSA2d.reconstruction}(object, model = "additive", offset = 0, \dots)
\method{residuals}{BioSSA3d.reconstruction}(object, model = "additive", offset = 0, \dots)
}

\arguments{
  \item{object}{object of class `BioSSA2d', `BioSSA3d', `BioSSA2d.reconstruction'
    or `BioSSA3d.reconstruction'}
  \item{groups}{numeric vector, indices of elementary components}
  \item{model}{model name ('additive', 'multiplicative' or 'poisson') or multiplicity power value}
  \item{offset}{numeric value, trend offset}
  \item{\dots}{additional arguments, passed to \code{reconstruct}}
}

\value{
  Object of class `embryo2d' or `embryo3d' for 2d or 3d method correspondingly.
}

\seealso{
  \code{\link[BioSSA:noise.model]{noise.model}}
}
