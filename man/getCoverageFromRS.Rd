\name{getCoverageFromRS}
\Rdversion{1.1}
\alias{getCoverageFromRS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCoverageFromRS - conversion to coverage object
}
\description{
Calculates the coverage function for the reads encapsulated in the SeqReads object. 
}
\usage{
getCoverageFromRS(rs, exps)
}

\arguments{
  \item{rs}{SeqReads object to modify
}
  \item{exps}{Vector of numbers of experimental samples in xXMAP
}
}

\value{
   NucleotideDistr object with coverage matrix in assayData slot and type "COV".
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  if (xmapConnected())  
#  {
#    rs <- newSeqReads(1,1,20000,1)
#    rs <- addExperimentsToReadset(rs,1:6)
#    nd.cov <- getCoverageFromRS(rs,1:3)
#  }
}
