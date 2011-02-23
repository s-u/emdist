\name{emd}
\alias{emd}
\alias{emd3d}
\alias{emd3dw}
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
emd3d(A, B)
emd3dw(A, wA, B, wB, ...)
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
  \item{wA}{weights for locations specified by A}
  \item{wB}{weights for locations specified by B}
  \item{...}{unused}
}
\details{
  \code{emd} interprets the two matrices \code{A} and \code{B} as
  a distibution over a two-dimensional grid. The distance between the
  gridn points in each direction is defined by \code{xdist} and
  \code{ydist}. Both matrices must have the same dimensionality.

  \code{emd3d} uses first column of each matrix as the weigths and the
  remaining columns as location coordinates in a up to three-dimensional
  space. \code{A} and \code{B} must have the same number of columns.
  
  \code{emd3dw} separates the weights from the location matrices but is
  otherwise identical to \code{emd3d}.

  Note that in all cases the weights are transformed to integers by a
  multiplication with 1e6 since the underlying code is working on bin
  counts. The implication is that weights cannot exceed approximately
  the value of 2000.
}
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
# and the 3d code should, hopefully have the same results
A. <- matrix(c(1:6 / sum(1:6), 1,1,2,2,3,3, 1,2,1,2,1,2), 6)
B. <- matrix(c(1, 2, 3), 1)
emd3d(A., B.)
}
\keyword{multivariate}
