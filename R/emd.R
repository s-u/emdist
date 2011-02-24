emd2d <- function(A, B, xdist = 1, ydist = 1) .Call("emd_2d", A, B, xdist, ydist, PACKAGE="emd")

emd <- function(A, B, ...) .Call("emd_4d", A, B, PACKAGE="emd")

emdw <- function(A, wA, B, wB, ...) emd3d(cbind(wA, A), cbind(wB, B), ...)
