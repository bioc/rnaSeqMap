\name{bam2sig}
\alias{bam2sig}

\title{ bam2sig - encapsulated pipeline of finding significant expression
}
\description{
Reads BAM files according to annotation and produces output table processed with DESeq negative binomial test. 
}
\usage{
bam2sig(annotlib, covdesc="covdesc", species=NULL, level="gene")
}
\arguments{
  \item{annotlib}{Character table or data frame with colums: chr, start, end, strand, name 
}
  \item{covdesc}{Name of the file that includes BAM files (experiment description file) 
}
  \item{species}{Species name - needed for .chr.convert function - to match BAM and annotation chromosome names 
}
  \item{level}{The level of annotation for calculating the counts: gene, transcript of exon
}
}
\value{ Output table with significant expression results, as from DESeq
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (1==0)  
  {
   all.g <- all.genes(as.vector=F)
   ss <- sample(1:20000, 10)
   genes <- as.data.frame(all.g[ss,])

   deseqRes <- bam2sig("cassava.db")
   deseqRes[1:10,]
  }
}
