\name{SeqReads}
\Rdversion{1.1}
\alias{SeqReads}
\alias{SeqReads-class}
\alias{newSeqReads}
\alias{newSeqReadsFromGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SeqReads - a container for RNAseq reads
}
\description{
SeqReads objects keep the reads information in the form of a list, containing one matrix of reads per experiment. Matrices of dimension n x 2 should come from a mapping to the regions defined by genome coordinates (chromosome, start, end, strand) in the SeqReads object. 

The object may be filled in from the database or from list with read data. It is recommended to create one SeqReads object per gene or intergenic region. The object are used then ot create object of class NucleotideDistr

}
\usage{
newSeqReads(chr, start, end, strand, datain=NULL)
newSeqReadsFromGene(g)
}

\arguments{
  \item{chr}{Chromosome
}
  \item{start}{Start of the region on a chromosome
}
  \item{end}{End of the region on a chromosome
}
  \item{strand}{Genome strand: 1 or -1
}
  \item{datain}{If supplied, it must be a list of matrices of reads start and stop
}
  \item{g}{Ensembl identifier of a gene
}
}

\value{
   Object of a class SeqReads
}
\author{
   Michal Okoniewski, Anna Lesniewska
}
