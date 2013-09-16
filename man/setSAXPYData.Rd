\name{setSAXPYData}
\alias{setSAXPYData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data accessor function for rnaSeqMap objects containing  'data' field
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function sets the 'data' field at i position. The new value is the old one multiplied by a iParam.
}
\usage{
setSAXPYData(iND1,iParam,i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iND1}{
 rnaSeqMap object containing 'data' field
}
\item{iParam}{
Scalining parameter
}
\item{i}{Index of the 'data' field to be modified}
}

\value{
NULL
}
\author{
   Michal Okoniewski, Anna Lesniewska, Marek Wiewiorka
}
