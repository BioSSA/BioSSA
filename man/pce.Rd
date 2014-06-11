\name{pce}
\alias{read.pce.gene}
\alias{read.pce}
\alias{extract.gene.pce}
\title{Read gene data from PCE file}

\description{
  Read gene expression data from PointCloud File
}

\usage{
read.pce(file)
read.pce.gene(file, gene, \dots)
extract.gene.pce(pce, gene, measurement = c("default", "apical", "basal", "nuclear", "cellular"))
}

\arguments{
  \item{file}{character, PCE file name}
  \item{gene}{character, gene name for extraction}
  \item{\dots}{further arguments passed to \code{extract.gene.pce} (e.g. \code{measurement})}
  \item{measurement}{measurement kind}
  \item{pce}{parsed PCE content returned by \code{read.pce}}
}

\details{
  \code{read.pce.gene} is combination of \code{read.pce} and \code{extract.gene.pce}.
}

\value{
  Object of class `embryo3d' for \code{read.pce.gene} and \code{extract.gene.pce}.
  Data frame with additional attributes for \code{read.pce}.
}

\seealso{
  \code{\link[BioSSA:embryo3d]{embryo3d}},
  \code{\link[BioSSA:plot.embryo3d]{plot.embryo3d}}
}

\references{
  PCE format description: \url{http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=pceFormat}
}

\examples{
file <- system.file("extdata/pce", "v5_s12166-23au07-22.pce", package = "BioSSA")
pce <- read.pce.gene(file, gene = "eve")
plot(pce, type = "hull")
plot(pce, type = "nuclei")
}
