\name{BioSSA}
\alias{BioSSA}
\alias{BioSSA.formula}
\alias{BioSSA.embryo2d}
\alias{BioSSA.embryo3d}
\title{Construct a new BioSSA decomposition object}

\description{
  Construct a new BioSSA decomposition object from passed embryo object or formula.
}

\usage{
\method{BioSSA}{formula}(x, data = NULL, ...)
\method{BioSSA}{embryo2d}(x,
       ...,
       step = 0.5,
       L = c(50, 50),
       xlim = c(-Inf, Inf),
       ylim = c(-Inf, Inf),
       xperc = c(0, 100),
       yperc = c(0, 100),
       units = c("percent", "original"))
\method{BioSSA}{embryo3d}(x, ..., sweep = sweep.function)
}

\arguments{
  \item{x}{formula or object of class `embryo2d' or `embryo3d'}
  \item{\dots}{further arguments passed to \code{ssa} with parameter
    \code{kind = '2d-ssa'}; \code{ssa} is the main function of the `Rssa' package,
    which performs Singular Spectum Analysis}
  \item{units}{units for \code{L}, \code{step}, \code{xlim} and \code{ylim} arguments. `percent'
    means percent of whole embryo length}
  \item{L}{vector of two elements, window length (in `units')}
  \item{step}{grid step for interpolation (in `units')}
  \item{xlim, ylim}{numeric vectors of two elements, cutoff bounds by AP and DV axes (in `units')}
  \item{xperc, yperc}{numeric vectors of two elements,
    denote the size (width and height) of presented embryo part in percent
    (usefull if data present only a part of the whole embryo)}
  \item{data}{a data frame (or more precisely, anything that is a valid \code{envir} argument in \code{eval},
     e.g., a list or an environment) containing values for any variables in the formula.
     If not found in \code{data}, or if \code{data} is unspecified, the variables are looked for
     in the environment of the formula}
  \item{sweep}{sweeping function for planarize 3d embryo to 2d ones.
    Just now is only default function (\code{sweep.function}) is supported, which provides
    cylindrical sweeping (MB, projection; please, somebody, fix my lexics\dots). Interface for
    user-specified planarization functions in developing now (In future, it will be function which
    takes object of class `embryo3d' and returns list of `embryo2d' objects)}
}

\details{
  For the `formula': formula is evaluated with using \code{data}.
  The specified data are transformed to `embryo2d' or `embryo3d' object (it's determined by formula length).
  the obtained object is passed (along with arguments specified by \dots) to the \code{BioSSA} method for
  corresponding class.

  For the `embryo3d' object: 2d projection is performed (using \code{sweep} function, topology is kept,
  that is, along Y (DV) axis the circulat topology is considered. Thus, `embtyo3d' contains
  `embryo2d' object. Then the \code{BioSSA} method for `embryo2d' class is called.

  For the `embryo2d' object: the linear interpolation (see \code{\link[akima:interp]{interp}})
  is performed to the regular grid with step \code{step}.
  For regular data, the conventional 2D-SSA procedure (\code{ssa} from `Rssa') is called.
}

\value{
  `embryo2d' or `embryo3d' object, according to the type of \code{x},
  stores the result of decomposition of the trajectory matrix.
  To obtain the decomposition of the data, \code{\link[BioSSA:reconstruct]{reconstruct}} should be called.
}

\seealso{
  \link[BioSSA:BioSSA-class]{BioSSA-class},
  \code{\link[BioSSA:reconstruct]{reconstruct}},
  \code{\link[Rssa:ssa]{ssa}},
  \code{\link[Rssa:plot.ssa]{plot.ssa}}.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bss <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
summary(bss)
plot(bss)
plot(bss, type = "ssa-vectors")
plot(bss, type = "ssa-vectors", vectors = "factor")
}
