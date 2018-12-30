


CalculateC <- function(a, b, epsAlpha = 1e-07, AlphaHandler = warning) {
  if (is.null(epsAlpha)) 
    return(chol(t(a) %*% a - t(b) %*% b))
  matrixC <- try(chol(t(a) %*% a - t(b) %*% b), silent = TRUE)
  if (!inherits(matrixC, "try-error")) {
    attr(matrixC, "alpha") <- 1
    return(matrixC)
  }
  alpha <- min(1, FindAlpha(a, b)/(1 + epsAlpha))  #
  if (alpha < 1) 
    AlphaHandler(paste("alpha = ", alpha))
  
  matrixC <- chol(t(a) %*% a - alpha^2 * t(b) %*% b)
  attr(matrixC, "alpha") <- alpha
  matrixC
}



#' Calculation of alpha
#' 
#' Function to find the largest alpha that makes equation 10 in the paper solvable.
#'
#' @param a matrix E in paper 
#' @param b matrix Eg in paper
#'
#' @return alpha
#' 
#' @note FindAlphaSimple performs the calculations by a simple/direct method. 
#' FindAlpha is made to handle problematic special cases.
#'  
#' @export
#'
FindAlpha <- function(a, b) {
  aCol <- IndependentCol(a)
  bCol <- IndependentCol(b)
  cCol <- aCol | bCol
  if (any(!cCol)) {
    stop("Problematic collinearity detected. Not handled in this version.")
  }
  ata <- t(a) %*% a
  btb <- t(b) %*% b
  if (!any(aCol != cCol)) {
    # ata invertible
    m <- solve(ata, btb)
    return(sqrt(MinEigen(m, inverse = TRUE)))
  }
  if (!any(bCol != cCol)) {
    # btb invertible
    m <- solve(btb, ata)
    return(sqrt(MinEigen(m, inverse = FALSE)))
  }
  stop("Could not calculate alpha.")
}


#' @rdname FindAlpha
#' @export
FindAlphaSimple <- function(a, b) {
  m <- t(a) %*% a %*% solve(t(b) %*% b)
  sqrt(min(eigen(m, only.values = TRUE)$values))
}

IndependentCol <- function(x, tol = 1e-07) {
  a <- rep(FALSE, NCOL(x))
  qrX <- qr(x, tol = tol)
  a[qrX$pivot[seq_len(qrX$rank)]] <- TRUE
  a
}

# Use of tol inspierd by ginv i MASS
MinEigen <- function(x, tol = sqrt(.Machine$double.eps), inverse = TRUE, ...) {
  a <- Re(eigen(x, symmetric = FALSE, only.values = TRUE, ...)$values)
  if (inverse) {
    ok <- abs(a) > max(tol * max(abs(a)), 0)
    a <- 1/a[ok]
  }
  min(a)
}

