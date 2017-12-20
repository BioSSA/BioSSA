\name{noise.model}
\alias{noise.model}
\alias{noise.model.default}
\alias{noise.model.BioSSA2d}
\alias{plot.noise.model}
\alias{print.noise.model}
\alias{summary.noise.model}
\title{Estimate model parameters for noise}

\description{
  Functions for storing, plotting and model construction for noise
}

\usage{
\method{noise.model}{default}(x, trend,
        offset = 0,
        model = "estimate",
        reg.type = c("winsor", "trim"),
        reg.level = 0,
        averaging.type = c("sliding-window", "quantile-break", "equal-break", "none"),
        breaks = 2, window = 51,
        FUN = median,
        FUN.trend = FUN,
        \dots)
\method{noise.model}{BioSSA2d}(x, groups, \dots)
\method{plot}{noise.model}(x,
        absolute = FALSE, #TODO Mb, remove it????
        relative = TRUE,
        draw.residuals = TRUE,
        draw.means.fitted = FALSE,
        print.alpha = TRUE,
        ref = TRUE,
        symmetric = !absolute,
        \dots,
        dots.residuals = list(),
        dots.means.fitted = list(),
        digits = max(3, getOption("digits") - 3))
\method{summary}{noise.model}(object, digits = max(3, getOption("digits") - 3), \dots)
\method{print}{noise.model}(x, digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{x, object}{the input object. This might be `BioSSA2d' object for \emph{BioSSA2d} method,
    or just a numeric vector of residuals for \emph{default} implementation}
  \item{trend}{numeric vector, trend for noise model estimation; this parameter can be used only
    for \emph{default} method}
  \item{offset}{numeric value, trend offset for noise model estimation}
  \item{model}{model name ('additive', 'multiplicative' or 'poisson') or multiplicity power value.
    Use 'estimate' for perform model estimation}
  \item{reg.type}{regularization type for residuals and trend values averaging}
  \item{reg.level}{numeric value, quantile level for regularization}
  \item{averaging.type}{character, averaging method for power estimation and model visualization}
  \item{breaks}{the number of intervals to trend values breaking. Used for 'equal-break'
    and 'quantile-break' averaging types}
  \item{window}{length of sliding window. Used for 'sliding-window' averaging type}
  \item{FUN}{averaging function for logarithms of residuals absolute values
    (like \code{\link[base:mean]{mean}} or \code{\link[stats:median]{median}})
    should take numeric vector and return one numeric value. Used for all averaging types except 'none'}
  \item{FUN.trend}{averaging function for logarithms of trend values, like previous}
  \item{\dots}{additional arguments, passed to inner function calls}
  \item{groups}{numeric vector, component indices in BioSSA2d decomposition for trend extraction}
  \item{absolute}{logical, whether plot absolute values of residuals instead of residuals itself}
  \item{relative}{logical, whether plot relative residuals instead of absolute ones}
  \item{draw.residuals}{logical, whether plot residuals}
  \item{draw.means.fitted}{logical, whether plot regression line drawn by averaged residuals}
  \item{print.alpha}{logical, whether output multiplicity power estimation}
  \item{ref}{logical, whether plot zero level line}
  \item{symmetric}{logical, whether y-scale should be symmetric by default}
  \item{dots.residuals}{list of additional arguments passed to \code{\link[lattice:xyplot]{xyplot}} function
    used for plotting of the residuals}
  \item{dots.means.fitted}{same as previous, but for regression line plot (not used, supposed to be excluded)}
  \item{digits}{integer, how many significant digits are to be used for numbers formatting}
}

\details{
  \code{noise.model} applied to `BioSSA2d' object calculates the residual and trend
  and then call the \emph{default} method passing the further arguments.

  \code{plot} method plots residuals and smoothed residuals in dependence on the trend.

  Data consisting of trend and residual values are ordered by the trend values,
  then \code{FUN} function is applied to residuals and the result is depicted
  again \code{FUN.trend} applied to the ordered trend values.
}

\value{
   Object of class `noise.model' for \code{noise.model}. Trellis object for \code{plot}.
   Original object (invisibly) for \code{print} and \code{summary}.

   Object of `noise.model' is a list with following fields:
   \describe{
     \item{alpha}{is numeric value, which is estimated if \code{model = 'estimate'},
       or is equal to the user specified value (multiplicity power)}
     \item{sigma}{is standard deviation of the relative noise, (see `Model description')}
     \item{sd}{is square root of mean square deviation of relative residuals}
     \item{residuals}{is vector of residuals}
     \item{trend}{is vector of trend values}
     \item{residuals.means}{is vector of smoothed (averaged) residuals absolute values}
     \item{residuals.means.fitted}{is vector of fitted values (by \code{lm})
       of smoothed residuals absolute values}
     \item{trend.means}{is vector of smoothed trend values}
     \item{offset}{is offset value}
     \item{averaging.type}{is character name of used averaging type}
     \item{call}{is object of class `call', constructor call (MB, drop it? This is command, which
       has created object, `call' means `function call')}
   }
}

\section{Model description}{
  Generalized multiplicative noise model is considered:
  \deqn{res_i = \sigma (trend_i + offset)^\alpha \cdot \xi_i,}
  where \eqn{\xi_i} have standard normal distribution.

  Value of \eqn{\alpha}{alpha} is estimated as follows:
  logarithms of absolute values of residuals and ofsetted trend are considered,
  then they both ordered by the trend values and averaging procedure performed.
  There are following averaging methods:
  'none' means nothing averaging,
  'sliding-window' means sliding window averaging (default approach)
  with window length denoted by \code{window} argument,
  'equal-break' and 'quantile-break' mean splitting all trend values
  into buckets ('quantile-break' means buckets with equal quantity of elements in each bucket and
  'equal-break' means equal size buckets), correspondingly residuals splitting and averaging
  logarithmed residuals and trend values in each bucket.

  Then linear regression on averaged logarithms is provided. Slope is \eqn{\alpha}{alpha} estimation
  and intercept is estimation of logarithm of the standard deviation of relative noise, i.e.
  \eqn{\sigma = e^{Intercept}}{sigma = exp(Intercept)}.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bs <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
nm <- noise.model(bs, 1:3, averaging.type = "none")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, averaging.type = "sliding")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, averaging.type = "equal")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, averaging.type = "quantile")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, model = "poisson")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, model = "additive")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, model = "multiplicative")
plot(nm)
summary(nm)

nm <- noise.model(bs, 1:3, model = -1.2)
plot(nm)
summary(nm)



nm.none <- noise.model(bs, 1:3, model = "estimate", averaging.type = "none")
nm <- noise.model(bs, 1:3, model = nm.none$alpha, averaging.type = "sliding")

plot(nm)
}

