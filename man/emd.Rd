\name{emd}
\alias{emd}
\alias{emd2d}
\alias{emdw}
\title{
Earth Mover's Distance
}
\description{
\code{emd} computes Earth Mover's Distance (related to 1st Mallows and
Wasserstein distances) between distributions. \code{emd} and \code{emdw}
use (weight,location) notation whereas \code{emd2d} compares two
distributions over a grid represented as matrices.
}
\usage{
emd(A, B, ...)
emdw(A, wA, B, wB, ...)
emd2d(A, B, xdist = 1, ydist = 1)
}
\arguments{
  \item{A}{matrix A}
  \item{B}{matrix B}
  \item{wA}{weights for locations specified by A}
  \item{wB}{weights for locations specified by B}
  \item{xdist}{
    distance between columns (scalar) or a vector of positions of the columns
  }
  \item{ydist}{
    distance between rows (scalar) of a vector of positions of the rows
  }
  \item{...}{additional parameters for future use, currently unused}
}
\details{
  \code{emd2d} interprets the two matrices \code{A} and \code{B} as
  a distibution over a two-dimensional grid. The distance between the
  grid points in each direction is defined by \code{xdist} and
  \code{ydist}. Both matrices must have the same dimensionality.

  \code{emd3d} uses first column of each matrix as the weigths and the
  remaining columns as location coordinates in a up to four-dimensional
  space. \code{A} and \code{B} must have the same number of columns.
  
  \code{emdw} separates the weights from the location matrices but is
  otherwise identical to \code{emd}.

  Note that in all cases the weights are transformed to integers by a
  multiplication with 1e6 since the underlying code is working on bin
  counts. The implication is that weights cannot exceed approximately
  the value of 2000 (which should have no real implications as the
  weight should sum up to one!).
}
\value{
  EMD of the distributions \code{A} and \code{B} (it is assumed that
  each have the sum of 1).
}
%\references{
%}
\author{
  R code by Simon Urbanek, EMD code by Shankar Krishnan,
  solver by David Applegate.
}
%\note{
%}
\examples{
A <- matrix(1:6 / sum(1:6), 2)
B <- matrix(c(0, 0, 0, 0, 0, 1), 2)
emd2d(A, B)
# if we bring the rows closer, the distance will be reduced
# since mass from the first row has to move down
emd2d(A, B,, 0.1)
# the positions can overlap - this is a degenerate case of that
emd2d(A, B, rep(0, 3), rep(0, 2))
# and just a sanity check
emd2d(A, A) + emd2d(B, B)
# and the weight/location code should, hopefully have the same results
A. <- matrix(c(1:6 / sum(1:6), 1,2,1,2,1,2, 1,1,2,2,3,3), 6)
B. <- matrix(c(1, 2, 3), 1)
stopifnot(emd(A., B.) == emd2d(A, B))
stopifnot(emd(A., B.) == emdw(A.[,-1], A.[,1], B.[,-1,drop=FALSE], B.[,1]))
}
\keyword{multivariate}
