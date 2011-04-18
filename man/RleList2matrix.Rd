\name{RleList2matrix}
\alias{RleList2matrix}
\title{RleList2matrix}
\description{ Function transforms list of Rle objects to matrix.
}
\usage{
RleList2matrix(list);
}
\arguments{
  \item{list}{list of Rle objects. 
}
  
}
\value{
 Produces the full, unpacked coverage matrix from a list of Rle objects. 
 Used to re-format the coverage data. 
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    rs <- newSeqReads(1,1,20000,1)
    rs <- getBamData(rs,1:3)
    nd.cov <- getCoverageFromRS(rs,1:3)
    RleList2matrix(nd.cov@data)  
  }
}

