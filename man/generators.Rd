\name{generators}
\alias{generatorAddSquare}
\alias{generatorAdd}
\alias{generatorMultiply}
\alias{generatorTrunc}
\alias{generatorSynth}
\alias{generatorPeak}

\title{ Generators for synt data and 
}
\description{Various generators for experiments.
}
\usage{

generatorAddSquare(nd, deg, length.prop=0.5)
generatorAdd(nd, deg, length.prop=0.5)
generatorMultiply(nd, deg, length.prop=0.5)
generatorTrunc(nd,deg)
generatorSynth(nd, deg, length.prop=0.5)
generatorPeak(nd, deg, sr=10, mult=10)

}

\arguments{
  \item{nd}{nucleotide distribution object 
}
  \item{deg}{ degeneration level for the output profile
}
  \item{length.prop}{ a fraction of the genome region to be degenerated - (0,1)
}
  \item{sr}{ distance from the 5' end for the peak
}
  \item{mult}{ multiplier - how many times the peak is supposed to be higher than the 
  maximum of the distribution
}
Generators of synthetic and semi-synthetic coverage profiles, for RNA-seq 
measures testing. 
}

\author{
    Anna Lesniewska,Michal Okoniewski
}
\examples{
   if (1==0)
   {
   rs <- newSeqReads('chr2', 220238268, 220254744, -1)
   f <- c("test1.bam", "test2.bam", "test3.bam", "test4.bam", "test5.bam")
   ff <- sapply(f, function(x) system.file("extdata", x, package = "rnaSeqMap"))
   rs <- getBamData(rs, 1:5)    	 
   nd <- getCoverageFromRS(rs, 1:5) 
   generatorTrunc(nd,0.5)
   }
}
