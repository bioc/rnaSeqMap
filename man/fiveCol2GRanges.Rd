\name{fiveCol2GRanges}
\Rdversion{1.1}
\alias{fiveCol2GRanges}
\title{fiveCol2GRanges  
}
\description{ Set of functions to operate on \code{NucleotideDistr} objects.

averageND calculates the mean for samples, sumND adds up selected samples' distributions, combineND adds two objects with the same size of distribution matrix, log2ND transforms all numeric data in the object into log space. 

}
\usage{
fiveCol2GRanges(t)
}
\arguments{
  \item{t}{ A matrix or data frame including genomic regions in 5 columns: ID, chr/contig name, start, end, strand
}
}

\value{
  GenomicRanges object with the same values 
}
\author{
   Michal Okoniewski
}
