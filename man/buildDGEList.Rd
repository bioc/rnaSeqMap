\name{buildDGEList}
\alias{buildDGEList}

\title{ buildDGEList - create DGEList (edgeR)
}
\description{Creates \code{DGEList} from the data in the database using the list 
of genes supplied - for further analysis with edgeR
}
\usage{
buildDGEList(genes,exps,conds=NULL)
}

\arguments{
  \item{genes}{ vector of Ensembl gene IDs
}
  \item{exps}{ vector of experiments
}
  \item{conds}{ Vector of experimental condition descriptions for the samples
}

}

\value{\code{DGEList} object filled with the data of gene-level counts of reads
   
}
\seealso{buildDESeq}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    data(sample_data_rnaSeqMap)
    gg <- names(rs.list)
    cds <- buildDGEList(gg,1:6, c("a","b","b","a","a","b"))
  }
}
