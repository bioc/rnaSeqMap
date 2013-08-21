\name{getFCFromND}
\Rdversion{1.1}
\alias{getFCFromND }
\title{getFCFromND - calculating fold change of coverages 
}
\description{This function calculates the fold change of two sample coverages from a \code{NucleotideDistr} objects. The coverages are assumed to be after logarithmic transformation, so the function basically subtracts the value and generates new \code{NucleotideDistr} object with a single vector of fold changes. 
}
\usage{
getFCFromND(nd, exps)
}

\arguments{
  \item{nd}{ \code{NucleotideDistr} object with coverages
}
  \item{exps}{a pair of numbers of samples in the experiment
}
}

\value{
   \code{NucleotideDistr} object of type "FC" with a single vector of fold changes
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  if (xmapConnected())  
#  {
#    rs <- newSeqReads(1,1,20000,1)
#    rs <- addExperimentsToReadset(rs,1:3)
#    nd.cov <- getCoverageFromRS(rs,1:3)
#    nd.fc <- getFCFromND(nd.cov,c(1,3))
#  }
}
