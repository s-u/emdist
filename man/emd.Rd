\name{emd}
\alias{emd}
\title{
Earth Mover's Distance
}
\description{
\code{emd} computes Earth Mover's Distance (related to 1st Mallows and
Wasserstein distances) between distributions. In particular, the current
code works on a two dimensional grid of matching points.
}
\usage{
emd(A, B, xdist = 1, ydist = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{matrix A}
  \item{B}{matrix B}
  \item{xdist}{
    distance between columns (scalar) or a vector of positions of the columns
  }
  \item{ydist}{
    distance between rows (scalar) of a vector of positions of the rows
  }
}
%\details{
%}
\value{
  EMD of the distributions \code{A} and \code{B} (it is assumed that
  each have the sum of 1).
}
%\references{
%}
%\author{
%}
%\note{
%}
\examples{
A <- matrix(1:6 / sum(1:6), 2)
B <- matrix(c(0, 0, 0, 0, 0, 1), 2)
emd(A, B)
# if we bring the rows closer, the distance will be reduced
# since mass from the first row has to move down
emd(A, B,, 0.1)
# the positions can overlap - this is a degenerate case of that
emd(A, B, rep(0, 3), rep(0, 2))
# and just a sanity check
emd(A, A) + emd(B, B)
}
\keyword{multivariate}
