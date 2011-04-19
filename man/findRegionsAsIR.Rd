\name{findRegionsAsIR}
\Rdversion{1.1}
\alias{findRegionsAsIR}
\title{findRegionsAsIR - finding regions of high coverage using Lindell-Aumann algorithm.
}
\description{
The function is running Lindell-Aumann algorithm to find regions of irreducible expression on the coverage data in the \code{NucleotideDistr} object. The function may be used to find the location and boundaries of significant expression of exons and small RNA. 
}
\usage{
findRegionsAsIR(nd, mi, minsup=5, exp)
}

\arguments{
  \item{nd}{An object of \code{NucleotideDistr} class that has coverage values for a given region
}
  \item{mi}{The threshold of coverage that makes the region significant
}
  \item{minsup}{Minimal support of the numeric association rule - namely, in this case, the mininmal length of the discovered region
}
  \item{exp}{ Sample (experiment) number
}
}
\value{
   IRanges object with irreducible regions with high coverage. 
}
\examples{
  if (xmapConnected())  
  {
    rs <- newSeqReads(1,1,20000,1)
    rs <- addExperimentsToReadset(rs,1:3)
    nd.cov <- getCoverageFromRS(rs,1:3)
    nd.regs <- findRegionsAsND(nd.cov, 10)
  }
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
