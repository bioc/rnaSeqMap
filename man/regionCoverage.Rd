\name{regionCoverage}
\alias{regionCoverage}
\alias{gcoverage}

\title{ regionCoverage
}
\description{
Finds all the reads for a genomic range
}
\usage{
regionCoverage(chr, start, end, strand, ex, db = "FALSE" )
}

\arguments{
  \item{chr}{Chromosome
}
  \item{start}{Start of the region on a chromosome
}
  \item{end}{End of the region on a chromosome
}
  \item{strand}{Genome strand: 1 or -1
}
  \item{ex}{experiment
}
  \item{db}{ Use database (SQL) implementation of the algorithm
}

}

\value{coverage vector, independent from \code{NucleotideDistr}
}
\examples{
  if (xmapConnected())  
  {
   tmp <- regionCoverage( 1, 10000, 20000, 1,3)
  }
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
