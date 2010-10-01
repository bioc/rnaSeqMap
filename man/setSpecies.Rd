\name{setSpecies}
\alias{setSpecies}

\title{ setSpecies
}
\description{
Sets the species name for chromosomes X, Y and MT translation into consecutive numbers. If you use \code{xmap.connect}, no need to call \code{setSpecies}. Both set the internal variable of \code{xmapcore}.
}
\usage{
setSpecies(name=NULL)
}

\arguments{
  \item{name}{Species name
}
}
\examples{
  setSpecies("mus_musculus")
}

\author{
   Michal Okoniewski, Anna Lesniewska
}
