\name{NucleotideDistr-class}
\alias{NucleotideDistr-class}
\alias{newNuctleotideDistr}
\alias{distribs}
\alias{getDistr}
\docType{class}
\title{Numeric distributions by nucleotide - class}

\description{
An S4 class that inherits from \code{eSet} and holds all the numeric distributions of functions defined over the genome. 
The values may include coverage, splicing, fold change, etc. for a region defined by 
genomic coordinates. 
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{chr}: numeric matrix containing the read counts.
  
  \code{start}: data.frame containing the library size and group labels.

  \code{end}: data.frame containing the library size and group labels.

  \code{strand}: data.frame containing the library size and group labels.

  \code{start}: data.frame containing the library size and group labels.

}

\seealso{SeqReads, NDtransforms}

\section{Methods}{
 
  \code{distribs} gives the matrix of distributions from \code{assayData}
   
  \code{getDistr} gives a single distributions from \code{assayData} as a vector

  \code{newNuctleotideDistr (distribs, chr, start, end, strand, type="UNKNOWN", phenoData=NULL, featureData=NULL)} constructor from a matrix of data and chromosome coordinates. 
}

\author{Anna Lesniewska, Michal Okoniewski}

\keyword{classes}
