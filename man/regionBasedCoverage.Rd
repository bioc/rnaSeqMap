\name{regionBasedCoverage}
\Rdversion{1.1}
\alias{regionBasedCoverage}
\title{regionBasedCoverage - transformation of the region coverage by the Lindell-Aumann regions
}
\description{
The function builds a \code{NucleotideDistr} object from another object of coverage, using sequential call of Lindell-Aumann algorithm on the same data with a sequence of mi-levels. Each nucleotide is assigned the maximum mi-value of a region that covers it. 

The output \code{NucleotideDistr} object has the distribution without peaks and small drops of coverage, but the thade-off is that the level of coverage are discrete: seq\*maxexp.
}
\usage{
regionBasedCoverage(nd, seqq=1:10, maxexp=20, minsup=5)
}

\arguments{
  \item{nd}{An object of \code{NucleotideDistr} class that has coverage values for a given region
}
  \item{seqq}{Vector of numbers used to divide the range of coverage for subsequent mi-levels
}
  \item{maxexp}{The maximal mi-level for coverage
}
  \item{minsup}{Minimal support of the numeric association rule - namely, in this case, the mininmal length of the discovered region
}
}

\value{
   NucleotideDistr object that includes a matrix with zeros where no region was found and a maximum of mi-levels used for the sequential region searched. The distributions are similar to coverage, but have removed outliers of coverage peaks and short drops of coverage.
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
  rs <- newSeqReads(1,1,20000,1)
  rs <- addExperimentsToReadset(rs,1:3)
  nd.cov <- getCoverageFromRS(rs,1:3)
  nd.regs <- regionBasedCoverage(nd.cov, 1:10, 100)
  #runs the Lindell-Aumann algorithm at 100, 90, ... and picks maximal mi-level, where the nucleotide has a region found
  }
}
