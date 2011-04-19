\name{buildDESeq}
\alias{buildDESeq}

\title{ buildDESeq - create CountDataSet 
}
\description{Creates \code{CountDataSet} from the data in the database using the list 
of genes supplied - for further analysis with DESeq
}
\usage{
buildDESeq(genes,exps,conds=NULL)
}

\arguments{
  \item{genes}{ vector of Ensembl gene IDs
}
  \item{exps}{ vector of experiments
}
  \item{conds}{ Vector of experimental condition descriptions for the samples
}

}

\value{\code{CountDataSet} object filled with the data of gene-level counts of reads
   
}
\seealso{buildDGEList}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
    data(sample_data_rnaSeqMap)
    gg <- names(rs.list)
    cds <- buildDESeq(gg,1:6, c("a","b","b","a","a","b"))
  }
}
