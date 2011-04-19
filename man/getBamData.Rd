\name{getBamData}
\alias{getBamData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getBamData - getting sample data from BAM file.  
}

\description{
Add data from experimental samples stored in BAM file. 
}

\usage{
getBamData(rs, exps=NULL, files=NULL, unstranded=FALSE, covdesc="covdesc")
}
\arguments{
  \item{rs}{SeqReads object to modify
}
  \item{exps}{Vector of numbers of experimental samples
}
  \item{files}{Vector of BAM files to read
}
  \item{unstranded}{Flag which type of data are using (with distinguishing strand or not)
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

