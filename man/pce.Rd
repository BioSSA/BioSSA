\name{pce}
\alias{read.pce.gene}
\alias{read.vpc.gene}
\alias{read.pce}
\alias{extract.gene.pce}
\alias{read.vpc}
\alias{extract.gene.vpc}
\title{Read gene data from PCE file}

\description{
  Read gene expression data from PointCloud file
}

\usage{
read.vpc(file)
read.pce(file)
extract.gene.vpc(vpc, gene, cohort)
extract.gene.pce(pce, gene, measurement = c("default", "apical", "basal", "nuclear", "cellular"))
read.vpc.gene(file, gene, \dots)
read.pce.gene(file, gene, \dots)
}

\arguments{
  \item{file}{character, PCE file name}
  \item{gene}{character, gene name for extraction}
  \item{vpc}{parsed VPC file content returned by \code{read.vpc}}
  \item{pce}{parsed PCE file content returned by \code{read.pce}}
  \item{cohort}{cohort number in VPC file}
  \item{measurement}{measurement kind}
  \item{\dots}{further arguments passed to \code{extract.gene.vpc} or \code{extract.gene.pce} (i.e., \code{cohort} or \code{measurement})}
}

\details{
  \code{read.pce.gene} is combination of \code{read.pce} and \code{extract.gene.pce}.
  \code{read.vpc.gene} is combination of \code{read.vpc} and \code{extract.gene.vpc}.
}

\value{
  Object of class `embryo3d' for \code{extract.gene.vpc}, \code{extract.gene.pce}, \code{read.vpc.gene}, and \code{read.pce.gene}.
  Data frame with additional attributes for \code{read.vpc} and \code{read.pce}.
}

\seealso{
  \code{\link[BioSSA:embryo3d]{embryo3d}},
  \code{\link[BioSSA:plot.embryo3d]{plot.embryo3d}}
}

\references{
  PCE format description: \url{http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=pceFormat}
  VPC format description: \url{http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=vpcFormat}
}

\examples{
file <- system.file("extdata/pce", "v5_s12166-23au07-22.pce", package = "BioSSA")
eve <- read.pce.gene(file, gene = "eve")
plot(eve, type = "hull")
plot(eve, type = "nuclei")


file <- system.file("extdata/vpc", "D_pse_wt_atlas_rc4.vpc", package = "BioSSA")
vpc <- read.vpc(file)
eve <- extract.gene.vpc(vpc, gene = "eve", cohort = 4)
plot(eve, type = "hull")
plot(eve, type = "nuclei")
}
