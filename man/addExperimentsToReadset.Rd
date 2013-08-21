\name{addExperimentsToReadset}
\Rdversion{1.1}
\alias{addExperimentsToReadset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{addExperimentsToReadset - getting sample data from the database.  
}
\description{
Add data from experimental samples in the xXMAP database to the readset. Connection to the database required. 
}
\usage{
addExperimentsToReadset(rs, exps)
}

\arguments{
  \item{rs}{SeqReads object to modify
}
  \item{exps}{Vector of numbers of experimental samples in xXMAP
}
}

\value{
   SeqReads object with samples added from the database.
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  if (xmapConnected())  
#  {
#    rs <- newSeqReads(1,1,20000,1)
#    rs <- addExperimentsToReadset(rs,1:3)
#  }
}
