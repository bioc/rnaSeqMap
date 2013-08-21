\name{parseGff3}
\alias{parseGff3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{parseGff3 - parsing gff3 file format 
}

\description{
Parses gff3 file into genes, transcripts and exons. 
}

\usage{
parseGff3(filegff, fileg="genes.txt", filet="transcripts.txt", filee="exons.txt", nofiles=FALSE)
}
\arguments{
  \item{filegff}{ Input file in GFF3 format
}
  \item{fileg}{Filename for output: genes 
}
  \item{filet}{Filename for output: transcripts
}
  \item{filee}{Filename for output: exons
}
  \item{nofiles}{Flag: just optput list, no files
}
}

\value{
     List with elements "genes", "transcripts", "exons" with appropriate tables. 
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
#  if (xmapConnected())  
#  {
#     parseGff3("Athaliana.gff3")
#  }
}

