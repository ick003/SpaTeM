% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getBonds}
\alias{getBonds}
\title{Function to get the bonds based on neighbourhood}
\usage{
getBonds(Location, NN, th = 2)
}
\arguments{
\item{Location}{matrix of coordinates.}

\item{NN}{maximum number of neighbours.}
}
\description{
Function to get the bonds based on neighbourhood
}
\examples{
Obs.loc = as.matrix(expand.grid(1:4,1:4))
getBonds(Obs.loc, NN = 4)
}
