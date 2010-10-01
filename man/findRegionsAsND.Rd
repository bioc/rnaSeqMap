\name{findRegionsAsND}
\Rdversion{1.1}
\alias{findRegionsAsND}
\alias{regionmining}
\title{findRegionsAsND - finding regions of high coverage using Lindell-Aumann algorithm.
}
\description{
The function is running Lindell-Aumann algorithm to find regions of irreducible expression on the coverage data in the \code{NucleotideDistr} object. The function may be used to find the location and boundaries of significant expression of exons and small RNA. 
}
\usage{
findRegionsAsND(nd, mi, minsup=5)
}

\arguments{
  \item{nd}{An object of \code{NucleotideDistr} class that has coverage values for a given region
}
  \item{mi}{The threshold of coverage that makes the region significant
}
  \item{minsup}{Minimal support of the numeric association rule - namely, in this case, the mininmal length of the discovered region
}
}

\value{
   NucleotideDistr object that includes a matrix with zeros where no region was found and the value of mi for all the nucleotides included in the region. The type fo the object is "REG".
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

