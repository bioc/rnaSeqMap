\name{readsInRange}
\alias{readsInRange}

\title{ readsInRange
}
\description{
Finds all the reads for a genomic range
}
\usage{
readsInRange(chr, start, end, strand, ex)
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

}

\value{table of reads, as in the database
}
\examples{
#  if (xmapConnected())  
 # {
#   tmp <- readsInRange( 1, 10000, 20000, 1,3)
#  }
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
