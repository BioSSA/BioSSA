\name{plot.embryo2d}
\alias{plot.embryo2d}
\title{Plot `embryo2d' object}
\description{
  The function plots `embryo2d' object in different ways:
  (1) 2d color levelplot of intesities,
  (2) 3d plot of intensities and
  (3) 1d plots of sections along coordinate (AP or DV) directions.
}

\usage{
\method{plot}{embryo2d}(x, type = c("nuclei-2d", "field-2d",
                                    "field-section", "nuclei-section",
                                    "field-3d", "nuclei-3d"),
                                  \dots)
}

\arguments{
  \item{x}{object of class `embryo2d'}
  \item{type}{plot type (see Details)}
  \item{\dots}{further arguments passed to plotting function}
}

\details{
  A lot of drawing methods, 1d, 2d, 3d, for regular grid and for orginal nuclear centers.
  Type variants:
  \describe{
    \item{nuclei-2d}{is 2d nuclei plot. One point for one nucleus. Also one can pass the following arguments:
      \code{voronoi} --- logical whether plot Voronoi diagram for nuclei, 'TRUE' by default;
      \code{col} --- argument for \code{\link[grDevices:colorRampPalette]{colorRampPalette}}, this pallete used
      for Voronoi diagram, color intensity corresponds to nuclei intensity value,
      greyscale pallete by default}
    \item{field-2d}{is Lattice levelplot for intensity regularization. Argument \code{col} defines color palette
      (see previous item)}
    \item{field-section}{is Lattice xyplot for of section for intensities in regular grid. Additional arguments:
      \code{at}, \code{units} --- units for \code{at} ('percent' or 'original'),
      \code{coord} --- 'x' or 'y', 'y' by default}
    \item{nuclei-section}{is Lattice xyplot of section for intensities in nucleus. Additional arguments:
      \code{at}, \code{units}, \code{coord} --- see previous item, \code{tolerance} ---
      is the half-width of the strip around the section position}
    \item{field-3d}{is 3d rgl plot, nuclei are depicted by colored spheres,
      color intensity corresponds to nuclei intensity value. Additional argument \code{col} ---
      character vector of two values, denotes colors for the lowest and the highest intensities correspondingly}
    \item{nuclei-3d}{is 3d rgl surface plot, values in regular grid are depicted,
      color intensity corresponds to nuclei intensity value. Additional argument \code{col} --- see previous item}
  }

  Other parameters are passed to plot functions.
  % Остальные параметры передаются дальше в функции
  % Тут у нас есть маленькое ТУДУ для меня. Регуляризация может отсутсвовать в эмбрионе,
  % одни голые ядра (если эмбрион не был обработан BioSSA). Поэтому сейчас некоторые рисовалки не работают


  % Не очень понятно, как с этим быть, параметры регулярицации задаются в BioSSA. Либо надо выносить регуляризацию
  % из BioSSA, либо надо отслеживать случай нерегуляризованного эмбриона
}

\value{
  `trellis' object (see \link[lattice:lattice-package]{lattice}) for 2d plots,
  a vector of object IDs for 3d plots.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)
good <- 3

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bss <- BioSSA(cad ~ AP + DV,  data = df,
              L = L,
              step = 0.5,
              xlim = xlim, ylim = ylim)


# Reconstruction of elementary components
rec.elem <- reconstruct(bss, groups = 1:6)
plot(plot(rec.elem))

# Sections for testing the reconstruction quality
rec <- reconstruct(bss, groups = list(good = 1:good))
p.ny <- plot(attr(rec, "series"), type = "nuclei-section", at = 50, coord = "y", tol = 5)
p.fy2 <- plot(rec$good, type = "field-section", at = 50, coord = "y")

# y-sections
pls <- list()
pls[[1]] <- p.ny + p.fy2
pls[[2]] <- plot(residuals(bss, 1:good), type = "nuclei-section",
                 at = 50, coord = "y", tol = 5,
                 ref = TRUE, col = "blue")

print(pls[[1]], split = c(1, 1, 2, 1), more = TRUE)
print(pls[[2]], split = c(2, 1, 2, 1))


# Dependence of noise on trend
plot(bss, type = "residuals", model = "multiplicative", groups = 1:good)

# 3d-figure of reconstruction
plot(rec$good, type = "nuclei-3d", col = c("blue", "red"))

# 2d-figures with triangulation
plot(rec$good, type = "nuclei-2d")
}
