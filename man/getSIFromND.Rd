\name{getSIFromND}
\Rdversion{1.1}
\alias{getSIFromND }
\alias{splicingind}
\title{getSIFromND - calculating splicing index of two coverages 
}
\description{This function calculates the splicing index value of two sample coverages from a \code{NucleotideDistr} object. It is assumed that the region in the \code{NucleotideDistr} is a single gene. 
Splicing index is calculated in similar way to the implementation for exon Affy microarrays (see Gardina et al, Genome Biology, 2007 for details), but it is run for each nucleotide in the region and instead of gene-level average expression values, it uses sums of reads for both samples. 
}
 
\usage{
getSIFromND(nd, exps)
}

\arguments{
  \item{nd}{ \code{NucleotideDistr} object with coverages
}
  \item{exps}{a pair of numbers of samples in the experiment
}
}

\value{
   \code{NucleotideDistr} object of type "FC" with a single vector of splicing index values 
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  if (xmapConnected())  
#  {
#    rs <- newSeqReads(1,1,20000,1)
#    nd.cov <- getCoverageFromRS(rs,1:3)
#    nd.fc <- getSIFromND(nd.cov,c(1,3))
#  }
}
