\name{plot.BioSSA2d}
\alias{plot.BioSSA2d}
\title{Plot method for class `BioSSA2d'}

\description{
  Plot the decomposition of 2d embryo.
  Plot graphs from the `Rssa' package for `ssa' object.
  Also, it can draw two plots of residuals:
  plot of the corresponding `noise.model' and
  scatterplot of residuals for different trends reconstructed by different sets of elementary components.
}

\usage{
\method{plot}{BioSSA2d}(x, type = c("ssa-values", "ssa-vectors", "ssa-series", "ssa-wcor",
                                      "residuals"), \dots)
}

\arguments{
  \item{x}{object of class `BioSSA2d'}
  \item{type}{plot type}
  \item{\dots}{additional arguments passed to plot functions}
}

\details{
  For the \code{type} parameter 'ssa-*' are just corresponding `Rssa' plots.

  For \code{type = 'residuals'} one can use following additional arguments:
  \describe{
    \item{groups}{numeric vector, elementary component indices for trend extraction}
    \item{groups.y}{like previous, used for `residuals-vs-residuals' plot}
    \item{lm}{logical, whether plot linear regression line (with zero offset) on `residuals-vs-residuals' plot}
  }

  If \code{groups.y} is missing then noise model will be estimated and plotted
  (additional arguments will be passed to both \code{noise.model} and \code{plot} methods).
  Otherwise, `residuals-vs-residuals' plot will be plotted.
}

\value{
  A Trellis object.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bss <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
plot(bss, type = "ssa-values")
plot(bss, type = "ssa-vectors")
plot(bss, type = "ssa-series", groups = 1:3)

plot(bss, type = "residuals", groups = 1:3)
plot(bss, type = "residuals", groups = 1:3, groups.y = 1:2)
}


