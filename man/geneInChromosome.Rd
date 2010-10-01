\name{geneInChromosome}
\alias{geneInChromosome}

\title{ geneInChromosome
}
\description{
Finds all the genes in the given chromosome regions
}
\usage{
geneInChromosome(chr, start, end, strand)
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

\value{ table of the genes in a given regions, produced with stored procedure
}
\examples{
  if (xmapConnected())  
  {
    geneInChromosome(1, 1, 80000, 1)
  }
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
