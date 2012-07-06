\name{measures}
\alias{ks_test}
\alias{pp_plot}
\alias{diff_area}
\alias{diff_derivative_area}
\alias{qq_plot}
\alias{qq_derivative_plot}
\alias{pp_derivative_plot}
\alias{hump_diff1}
\alias{hump_diff2}

\title{ Measures 
}
\description{Various measures to find differential expression.
}
\usage{

ks_test(dd)
diff_area(dd, cconst)
diff_derivative_area(dd, cconst)
qq_plot(dd)
qq_derivative_plot(dd)
pp_plot(dd)
pp_derivative_plot(dd)
hump_diff1(dd)
hump_diff2 (dd)

}

\arguments{
  \item{dd}{a matrix with 2 columns for samples and rows for nucleotides, containing 
  coverage data (like from BED files)}
  \item{cconst}{NULL default}

 
 The measures give various assesment of the difference between two sequencing samples 
 shapes. Full description will follow in the paper. 

}

\author{
    Anna Lesniewska,Michal Okoniewski
}

\examples{
 if (xmapConnected())  
  {
   ks_test(dd)
  }
}
