\name{plotGeneCoverage}
\alias{plotGeneCoverage}
\alias{plotRegionCoverage}
\alias{plotExonCoverage}
\alias{plotCoverageHistogram}
\alias{plotGeneExonCoverage}
\alias{plotSI}
\alias{ghistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Genomic plots with rnaSeqMap
}
\description{Various plots of genomic coverage for experiments.
}
\usage{

plotGeneCoverage(gene_id, ex)
plotRegionCoverage(chr, start, end, strand, ex) 
plotExonCoverage (exon_id,ex)
plotCoverageHistogram (chr,start,end,strand,ex, skip)
plotGeneExonCoverage(gene_id, ex)
plotSI(chr,start,end,strand, exp1, exp2 )

}

\arguments{
  \item{ex}{ vectors of experiment numbers to plot
}
  \item{exp1, exp2}{ experiment numbers for splicing index
}
  \item{gene_id}{ Ensembl gene ID
}
  \item{exon_id}{ Ensembl exon ID
}
  \item{chr}{ Chromosome
}
  \item{start}{ Start position of region on the chromosome
}
  \item{end}{ Start position of region on the chromosome
}
  \item{strand}{ Strand
}
  \item{skip}{ size of the bucket in histogram
}
}

\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    plotGeneCoverage( "ENSG00000144567", 1:3) # plotting FAM134A for experiments 1,2,3
    plotRegionCoverage( 2, 220040947, 220050201, 1, 1:3 ) # the same, using coordinates
  }
}
