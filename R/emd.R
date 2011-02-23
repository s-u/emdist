emd <- function(A, B, xdist = 1, ydist = 1) .Call("emd_2d", A, B, xdist, ydist, PACKAGE="emd")

emd3d <- function(A, B) .Call("emd_3d", A, B, PACKAGE="emd")

emd3dw <- function(A, wA, B, wB, ...) emd3d(cbind(wA, A), cbind(wB, B), ...)
