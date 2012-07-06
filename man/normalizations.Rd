\name{normalizations}
\alias{standarizationNormalize}
\alias{min_maxNormalize}
\alias{densityNormalize}
\alias{globalCountsNormalize}

\title{Normalization Methods 
}
\description{Various normalization methods.
}
\usage{

standarizationNormalize(nd) 
min_maxNormalize(nd)
densityNormalize(nd)
globalCountsNormalize(nd, sums)

}

\arguments{
  \item{nd}{nucleotide distribution object 
}
  \item{sums}{sum of reads in a sequencing sample 
}

Normalizations of a single coverage profile for multiple 
samples contained in the \code{NucleotideDistr} object. 
Full description will follow in a paper. 

}

\author{
    Anna Lesniewska,Michal Okoniewski
}
\examples{
 if (xmapConnected())
  {
  s <- newSeqReads('chr2', 220238268, 220254744, -1)
  f <- c("test1.bam", "test2.bam", "test3.bam", "test4.bam", "test5.bam")
  ff <- sapply(f, function(x) system.file("extdata", x, package = "rnaSeqMap"))
  rs <- getBamData(rs, 1:5, files = ff)    	
  nd <- getCoverageFromRS(rs, 1:5)
  min_maxNormalize(nd)
  }
}
