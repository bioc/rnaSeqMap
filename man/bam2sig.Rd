\name{bam2sig}
\alias{bam2sig}

\title{ bam2sig - encapsulated pipeline of finding significant expression
}
\description{
Reads BAM files according to annotation and produces output table processed with DESeq negative binomial test. 
}
\usage{
bam2sig(annot, covdesc="covdesc", species="homo_sapiens")
}
\arguments{
  \item{annot}{Character table or data frame with colums: chr, start, end, strand, name 
}
  \item{covdesc}{Name of the file that includes BAM files (experiment description file) 
}
  \item{species}{Species name - needed for .chr.convert function - to match BAM and annotation chromosome names 
}
}
\value{ Output table with significant expression results, as from DESeq
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
\examples{
  if (xmapConnected())  
  {
   all.g <- all.genes(as.vector=F)
   ss <- sample(1:20000, 10)
   genes <- as.data.frame(all.g[ss,])

   genes <- cbind(as.vector(genes[,"stable_id"]), as.vector(genes[,"space"]), as.vector(genes[,"start"]), as.vector(genes[,"end"]) ,as.vector(genes[,"strand"]))
  colnames(genes) <- c("name", "chr", "start", "end", "strand")

   deseqRes <- bam2sig()
   deseqRes[1:10,]
  }
}
