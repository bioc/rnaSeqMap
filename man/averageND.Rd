\name{averageND}
\Rdversion{1.1}
\alias{averageND }
\alias{sumND }
\alias{combineND }
\alias{log2ND }
\title{averageND, sumND, combineNS, log2ND - operations on distributions 
}
\description{ Set of functions to operate on \code{NucleotideDistr} objects.

averageND calculates the mean for samples, sumND adds up selected samples' distributions, combineND adds two objects with the same size of distribution matrix, log2ND transforms all numeric data in the object into log space. 

}
\usage{
averageND(nd, exps);
sumND(nd, exps);
combineND(nd1, nd2);
log2ND(nd);
}
\arguments{
  \item{nd, nd1, nd2}{ \code{NucleotideDistr} objects
}
  \item{exps}{a pair of numbers of samples in the experiment
}
}

\value{
   \code{NucleotideDistr} object of the same type as input objects
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    rs <- newSeqReads(1,1,20000,1)
    nd.cov <- getCoverageFromRS(rs,1:3)
    nd.avg <- averageND(nd.cov,c(1,3))
    nd.sum <- averageND(nd.cov,c(1,3))
    nd.sum <- combineND(nd.cov,nd.cov)
    nd.log <- log2ND(nd.cov)
  }
}
