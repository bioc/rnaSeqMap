\name{spaceInChromosome}
\alias{spaceInChromosome}

\title{ spaceInChromosome
}
\description{
Finds all the intergenic spaces in the given chromosome region
}
\usage{
spaceInChromosome(chr, start, end, strand)
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

}

\value{ table of the intergenic spaces in a given regions, produced with stored procedure
}
\examples{
  if (xmapConnected())  
  {
    spaceInChromosome(1, 1, 80000, 1)
  }
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
