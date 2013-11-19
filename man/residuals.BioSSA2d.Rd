\name{residuals.BioSSA2d}
\alias{residuals.BioSSA2d}

\title{
  Residuals of BioSSA2d decomposition
}

\description{
  Функция возвражает остатки. Опциоонально можно передать модель и сдвиг, тогда еще и поделит.
}

\usage{
\method{residuals}{BioSSA2d}(object, groups, \dots, model = "additive", offset = 0)
}

\arguments{
  \item{object}{object of class `BioSSA2d'}
  \item{groups}{numeric vector, indices of elementary components}
  \item{\dots}{additional arguments, passed to \code{reconstruct}}
  \item{model}{model name ('additive', 'multiplicative' or 'poisson') or multiplicity power value}
  \item{offset}{numeric value, trend offset}
}

\value{
  Object of class `embryo2d'.
}
