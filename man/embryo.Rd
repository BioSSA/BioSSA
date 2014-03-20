\name{embryo}
\alias{embryo2d}
\alias{embryo3d}
\alias{embryo2d.formula}
\alias{embryo3d.formula}
\title{Create a new `embryo' object}
\description{
  Creates a new object of class `embrio2d' or `embrio3d' correspondingly,
  input and output data for \code{BioSSA} procedure.

  ASH: Probably, we should join these methods into one and determine dimension by formula length.
}

\usage{
\method{embryo2d}{formula}(x, data = NULL, topology = c(Inf, Inf), ...)
\method{embryo3d}{formula}(x, data = NULL, ...)
}

\arguments{
  \item{x}{formula which describes 2d or 3d spatial embryo data}
  \item{\dots}{further arguments passed to \code{\link[stats:model.frame]{model.frame}} function}
  \item{data}{a data frame (or more precisely, anything that is a valid \code{envir} argument in \code{eval},
     e.g., a list or an environment) containing values for any variables in the formula.
     If not found in \code{data}, or if \code{data} is unspecified, the variables are looked for
     in the environment of the formula}
  \item{topology}{field topology, a vector of two positive (may be \code{Inf}) values.
    \code{Inf} means no cycling by corresponding coordinate and finite number denotes circle length}
}

\details{
  \code{topology} argument is used if 2d data given by cylindrical sweeping of 3d data.
}

\value{
  Object of class `embryo2d' or `embryo3d' correspondingly.

  `embrio2d' is list with following fields: two coordinates 'x2d', 'y2d' and gene activity 'value'.
  List has 'topology' attribute, which stores \code{topology}.
  If object returned from \code{reconstruct} or \code{residuals} methods of `BioSSA2d' class,
  it also contains field 'field', which stores interpolated values and interpolation grid parameters.

  `embrio3d' is list with following fields: three coordinates 'x3d', 'y3d', 'z3d' and gene activity 'value'.
}

\examples{
file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

emb2d <- embryo2d(cad ~ AP + DV, data = df)
plot(emb2d, "nuclei-2d")
}
