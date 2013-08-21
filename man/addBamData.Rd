\name{addBamData}
\alias{addBamData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{addBamData - getting sample data from BAM file.  
}

\description{
Add data from experimental samples stored in BAM file. 
}
\usage{
addBamData(rs, file, exp, phenoDesc=NULL)
}

\arguments{
  \item{rs}{SeqReads object to modify
}
  \item{file}{BAM file to read
}
  \item{exp}{Numbers of sample slot in the object
}
  \item{phenoDesc}{A vector to add to phenoData
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
#  if (xmapConnected())  
#  {
#    rs <- newSeqReads(1,1,20000,1)
#    rs <- addBamData(rs,1:3)
#  }
}

