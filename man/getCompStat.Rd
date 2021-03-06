% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getCompStat}
\alias{getCompStat}
\title{Function to summarize the components' statistics (size and colors).}
\usage{
getCompStat(Comp, CompID, Cols)
}
\arguments{
\item{Bds}{matrix of bonds.}

\item{nvert}{number of vertices.}
}
\description{
Function to summarize the components' statistics (size and colors).
}
\examples{
Obs.loc = expand.grid(1:10,1:10)
Bds = getBonds(Obs.loc)
Bds = Bds[which(Bds[,1]>0),]
Betas = c(0.8, 0.8, 0.2)
Cols = sample(1:3, 100, replace = TRUE)
Bd = selectBonds(Bds, Cols, Betas)
Bd = Bd[which(Bd[,1]>0),]
par(mfrow = c(1,1))
col = grey.colors(3)
image(matrix(Cols, ncol = 10), col = col)
segments((Obs.loc[Bd[,1],1]-1)/9, (Obs.loc[Bd[,1],2]-1)/9, (Obs.loc[Bd[,2],1]-1)/9, (Obs.loc[Bd[,2],2]-1)/9)
segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 33, col = "black")
segments((Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, (Obs.loc[,1]-1)/9, (Obs.loc[,2]-1)/9, lwd = 32, col = "white")
Comp = setComp(Bd, 100)
col.comp = rainbow(100)
text((Obs.loc[,1]-1)/div,(Obs.loc[,2]-1)/div, labels = 1:100, col = col.comp[Comp])
}
