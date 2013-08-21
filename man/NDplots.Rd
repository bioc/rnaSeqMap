\name{NDplots}
\alias{distrCOVPlot}
\alias{distrSIPlot}
\alias{distrCOVPlotg}

\title{ Genomic plots based upon NucleotideDistr objects
}
\description{Various plots of genomic coverage for data from \code{NucleotideDistr} objects
}
\usage{
distrCOVPlot(nd,exps)
distrSIPlot(nd, ex1, ex2, mi,minsup=5)
}

\arguments{
  \item{nd}{ \code{NucleotideDistr} object
}
  \item{exps}{ vectors of experiment numbers to plot
}
  \item{ex1,ex2}{ experiment numbers to plot
}
  \item{mi}{ threshold in the region mining algorithm
}
  \item{minsup}{ minimal support - minimal length of the irreducible region found
}

}

\author{
   Michal Okoniewski, Anna Lesniewska
}
%\examples{
%# data(sample_data_rnaSeqMap)
%# rs <- rs.list[[1]]  
% # if (xmapConnected())  
% # {
%#   nd.cov <- getCoverageFromRS(rs,1:6)
%#   distrSIPlot(nd.cov,  1,3, mi=5, minsup=10) 
% # }
%}
