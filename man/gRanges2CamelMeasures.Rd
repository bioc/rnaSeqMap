\name{gRanges2CamelMeasures}
\alias{gRanges2CamelMeasures}
\alias{allCamelMeasuresForRegion}

\title{ Genomic plots based upon NucleotideDistr objects
}
\description{Various plots of genomic coverage for data from \code{NucleotideDistr} objects
}
\usage{

gRanges2CamelMeasures(gR, cvd,  sample.idx1, sample.idx2, sums=NULL, progress=NULL)
allCamelMeasuresForRegion(ch, st, en, str, cvd, sample.idx1, sample.idx2, sums=NULL)
}

\arguments{
  \item{ch}{ chromosome name
}
  \item{st}{ genomic start
}
  \item{en}{ genomic end
}
  \item{str}{ strand
}
  \item{cvd}{ name of the file with BAM description - covdesc
}
  \item{gR}{ GenomicRanges object to use as a set of genomic regions to query
}
  \item{sample.idx1,sample.idx2}{ sample indices
}
  \item{sums}{ the vector of sums for normalization
}
  \item{progress}{ every how many regions print a dot for progress indicator
}

}

\author{
   Michal Okoniewski
}
\examples{
# 
  }
