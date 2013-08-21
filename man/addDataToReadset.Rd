\name{addDataToReadset}
\Rdversion{1.1}
\alias{addDataToReadset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{addDataToReadset - adding one more sample in the SeqRead on R level
}
\description{
Add another reads matrix to the readset. No control of region consistency, 
the matrix needs just 2 columns: starts and ends. 
}
\usage{
addDataToReadset(rs, datain, spl)
}

\arguments{
  \item{rs}{
}
  \item{datain}{
}
  \item{spl}{Number or name of the experimental sample
}
}

\value{
   SeqReads object with one more sample added. 
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  rs <- newSeqReads(1,1,20000,1)
#  my.data1 <- rbind(c(1,50), c(3,53), c(11,60))
#  rs <- addDataToReadset(rs, my.data1, 1)
}
