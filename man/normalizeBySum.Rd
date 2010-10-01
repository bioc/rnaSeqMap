\name{normalizeBySum}
\Rdversion{1.1}
\alias{normalizeBySum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalization of NucleotideDistr by global number of reads
}
\description{
\code{normalizeBySum} function normalizes the coverage values in \code{NucleotideDistr} by dividing all the numbers for all samples by the sum of reads for each sample. 
The number of reads from each sample may be taken from the database by the function getSumsExp, which is a wrapper for an appropriate SQL procedure. Alternatively, it is passed directly as a vector of numeric values of the same length as the number of samples analyzed.
Such simple normalization allows comparisons of the coverage values for samples with different number of reads
 }
\usage{
normalizeBySum(nd,r=NULL)
}

\arguments{
  \item{nd}{\code{NucleotideDistr} object with raw read counts
}
  \item{r}{Vector of numbers. If there is no such parameter, a database procedure summarizing reads is run
}
}

\value{
   NucleotideDistr object
}
\seealso{getSumsExp}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    rs <- newSeqReads(1,10000,20000,1)
    nd.cov <- getCoverageFromRS(rs,1:3)
    nd.norm <- normalizeBySum(nd.cov)
    nd.norm <- normalizeBySum(nd.cov, r=c(100, 200, 1000))
  }
}
