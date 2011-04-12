# this was the original API but now we just map it to the matrix API
emd2d <- function(A, B, xdist = 1, ydist = 1, ...) {
  if (!is.matrix(A) || !is.matrix(B) || !identical(dim(A), dim(B))) stop("A and B must be matrices of the same dimensions")
  m = dim(A)[1]
  n = dim(B)[2]
  A[is.na(A)] = 0
  B[is.na(B)] = 0
  A = matrix(c(A, rep((1:m) * ydist, n), rep((1:n) * xdist, each = m)), m * n)
  B = matrix(c(B, A[,2], A[,3]), m * n)
  emd(A[A[,1] != 0,,drop=FALSE], B[B[,1] != 0,,drop=FALSE])
}

emdr <- function(A, B, extrapolate=NA, flows=FALSE, ...) .Call("emd_r", A, B, extrapolate, flows, PACKAGE="emd")

emd <- function(A, B, ...) emdr(A, B, ...)

emdw <- function(A, wA, B, wB, ...) emd(cbind(wA, A), cbind(wB, B), ...)
