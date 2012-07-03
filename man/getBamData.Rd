\name{getBamData}
\alias{getBamData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getBamData - getting sample data from BAM file.  
}

\description{
Add data from experimental samples stored in BAM file. 
}

\usage{
getBamData(rs, exps=NULL,  rs, exps=NULL, cvd=NULL, covdesc.file="covdesc")
}
\arguments{
  \item{rs}{SeqReads object to modify
}
  \item{exps}{Vector of numbers of experimental samples
}
  \item{cvd}{Covdesc-like data frame - BAM files are read from row names 
}
  \item{covdesc}{Alternatively, the experiment description file
}
}

\value{
   SeqReads object with samples added from the BAM files. 
   List of BAM files comes from the covdesc. The covdesc 
   content becomes phenoData of the object. 
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    rs <- newSeqReads(1,1,20000,1)
    rs <- getBamData(rs,1:3)
  }
}

