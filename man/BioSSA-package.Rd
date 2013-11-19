\name{BioSSA-package}
\alias{BioSSA-package}
\docType{package}
\title{A collection of methods for singular spectrum analysis of drosophila embryo (2d and 3d) gene expression}

\description{
  This package provides routines for SSA analysis of gene expression.
}

\details{
This package is designed for analysis of
drosophila embryo gene expression given in 2d or 3d form.
Data should contain intensities of gene expression and coordinates, 2d or 3d correspondingly.

The main function is \code{\link[BioSSA:BioSSA]{BioSSA}},
which performs singular spectrum analysis for further decomposition on pattern and noise.
The data should be prepared as a data frame and the user can set the names of columns contining the data.
After \code{BioSSA} is applied, the reconstruction of the image components can be obtained by the
\code{\link[BioSSA:reconstruct]{reconstruct}} function.
Analysis of residuals is performed by \code{\link[BioSSA:noise.model]{noise.model}} methods.

S3-classes \link[BioSSA:BioSSA-class]{`BioSSA'} contain decompositions of 2d and 3d embryos.
S3-class \link[BioSSA:noise.model]{`noise.model'} constais residuals analysis results.

Visual presentation of the decomposition and residuals' analysis results can be performed by
\code{plot}-methods of corresponding classes.
}

\author{
Alex Shlemov, Nina Golyandina

Maintainer: Alex Shlemov <shlemovalex@gmail.com>
}

\references{
  Golyandina N. et al (2012):
  \emph{Measuring Gene Expression Noise in Early Drosophila Embryos: Nucleus-to-nucleus Variability.}
  \href{http://www.sciencedirect.com/science/article/pii/S1877050912001615}{Procedia Computer Science Volume 9,
    2012, Pages 373--382.}

  Golyandina N., Korobeynikov A., Shlemov A. and Usevich K. (2013):
  \emph{Multivariate and 2D Extensions of Singular Spectrum Analysis with the Rssa Package}
  \url{http://arxiv.org/abs/1309.5050}

  Shlemov A. and Golyandina N. (2014):
  \emph{Shaped extensions of singular spectrum analysis}
  \url{http://arxiv.org/abs/1401.4980}
}

\keyword{package}

\seealso{
  \code{\link[BioSSA:BioSSA]{BioSSA}},
  \code{\link[BioSSA:noise.model]{noise.model}}.
}

\examples{
bad <- 1
good <- 3
wave <- 3
xlim <- c(22, 88)
ylim <- c(32, 68)
atx <- 50
aty <- 50
tolx <- 5
toly <- 5
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bss <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
bss <- BioSSA(cad ~ AP + DV,  data = df,
              L = L,
              step = 0.5,
              xlim = xlim, ylim = ylim)

# w-correlations for identification
plot(plot(bss, type = "ssa-wcor", groups = 1:30))

# Reconstruction of elementary components
rec.elem <- reconstruct(bss, groups = 1:6)
plot(plot(rec.elem))

# Sections for testing the reconstruction quality
rec <- reconstruct(bss, groups = list(good = 1:good, bad = 1:bad))
p.ny <- plot(attr(rec, "series"), type = "nuclei-section", at = aty, coord = "y", tol = toly)
p.fy1 <- plot(rec$bad, type = "field-section", at = aty, coord = "y")
p.fy2 <- plot(rec$good, type = "field-section", at = aty, coord = "y")
p.nx <- plot(attr(rec, "series"), type = "nuclei-section", at = atx, coord = "x", tol = tolx)
p.fx1 <- plot(rec$bad, type = "field-section", at = atx, coord = "x")
p.fx2 <- plot(rec$good, type = "field-section", at = atx, coord = "x")


# y-sections, bad and good
pls <- list()
pls[[1]] <- p.ny + p.fy1
pls[[2]] <- plot(residuals(bss, 1:bad), type = "nuclei-section",
                 at = aty, coord = "y", tol = toly,
                 ref = TRUE, col = "blue")
pls[[3]] <- p.ny + p.fy2
pls[[4]] <- plot(residuals(bss, 1:good), type = "nuclei-section",
                 at = aty, coord = "y", tol = toly,
                 ref = TRUE, col = "blue")

print(pls[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pls[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pls[[3]], split = c(1, 2, 2, 2), more = TRUE)
print(pls[[4]], split = c(2, 2, 2, 2))

# x-sections, bad and good
pls <- list()
pls[[1]] <- p.nx + p.fx1
pls[[2]] <- plot(residuals(bss, 1:bad), type = "nuclei-section",
                 at = atx, coord = "x", tol = tolx,
                 ref = TRUE, col = "blue")
pls[[3]] <- p.nx + p.fx2
pls[[4]] <- plot(residuals(bss, 1:good), type = "nuclei-section",
                 at = atx, coord = "x", tol = tolx,
                 ref = TRUE, col = "blue")

print(pls[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pls[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pls[[3]], split = c(1, 2, 2, 2), more = TRUE)
print(pls[[4]], split = c(2, 2, 2, 2))

# Dependence of noise on trend
p1 <- plot(bss, type = "residuals", model = "additive", groups = 1:good)
p2 <- plot(bss, type = "residuals", model = "poisson", groups = 1:good)
p3 <- plot(bss, type = "residuals", model = "multiplicative", groups = 1:good)

print(p1, split = c(1, 1, 3, 1), more = TRUE);
print(p2, split = c(2, 1, 3, 1), more = TRUE);
print(p3, split = c(3, 1, 3, 1));


# 3d-figure of reconstruction
plot(rec$good, type = "nuclei-3d", col = c("blue", "red"))

# 2d-figures with triangulation
plot(rec$good, type = "nuclei-2d")
plot(rec.elem[[wave]], type = "nuclei-2d")

# Model estimation
noise.model(bss, 1:good)
}
