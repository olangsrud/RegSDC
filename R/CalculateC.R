


#' Calculation of C by solving equation 10 in the paper
#' 
#' The limit calculated by \code{\link{FindAlpha}} is used when alpha =1 cannot be chosen (warning produced). 
#' In output, alpha is attribute.
#'
#' @param a matrix E in paper 
#' @param b matrix Eg in paper
#' @param epsAlpha Precision constant for alpha calculation 
#' @param AlphaHandler Function (warning or stop) to be used when alpha<1 
#' @param alpha Possible with alpha as input instead of computing
#' 
#' @rdname CalculateC
#'
#' @return Calculated C with attributes alpha and viaQR (when CalculateC) 
#' 
#' @details  When epsAlpha=NULL calculations are performed directly (alpha=1) and alpha is not attribute. 
#' 
#' @export
#'
#' @examples
#' x <- 1:10
#' y <- matrix(rnorm(30) + 1:30, 10, 3)
#' a <- residuals(lm(y ~ x))
#' b <- residuals(lm(2 * y + matrix(rnorm(30), 10, 3) ~ x))
#' 
#' a1 <- a
#' b1 <- b
#' a1[, 3] <- a[, 1] + a[, 2]
#' b1[, 3] <- b[, 1] + b[, 2]
#' 
#' alpha <- FindAlpha(a, b)
#' FindAlphaSimple(a, b)  # Same result as above
#' CalculateC(a, b)
#' CalculateCdirect(a, b)  # Same result as above without viaQR attribute 
#' CalculateCdirect(a, b, alpha = alpha/(1 + 1e-07))  # Same result as above since epsAlpha = 1e-07
#' CalculateCdirect(a, b, alpha = alpha/2)  # OK
#' # CalculateCdirect(a,b, alpha = 2*alpha) # Not OK
#' 
#' FindAlpha(a, b1)
#' # FindAlphaSimple(a,b1) # Not working since b1 is collinear
#' CalculateC(a, b1, returnAlpha = TRUE)  # Almost same alpha as above (epsAlpha cause difference)
#' 
#' FindAlpha(b, a)
#' CalculateC(b, a, returnAlpha = TRUE)  # 1 returned (not same as above)
#' CalculateC(b, a)
#' 
#' FindAlpha(b1, a)   # alpha smaller than epsAlpha is set to 0 in CalculateC
#' CalculateC(b1, a)  # When alpha = 0 C is calculated by GenQR insetad of chol
CalculateCdirect <- function(a, b, epsAlpha = 1e-07, AlphaHandler = warning, alpha = NULL) {
  if (is.null(alpha)) {
    if (is.null(epsAlpha)) 
      return(chol(t(a) %*% a - t(b) %*% b))
    matrixC <- try(chol(t(a) %*% a - t(b) %*% b), silent = TRUE)
    if (!inherits(matrixC, "try-error")) {
      attr(matrixC, "alpha") <- 1
      return(matrixC)
    }
    alpha <- min(1, FindAlpha(a, b)/(1 + epsAlpha))
    if(alpha < epsAlpha)
      alpha <- 0
    if (alpha < 1) 
      AlphaHandler(paste("alpha = ", alpha))
  }
  if (alpha == 0) 
    matrixC <- GenQR(a)$R else matrixC <- chol(t(a) %*% a - alpha^2 * t(b) %*% b)
    attr(matrixC, "alpha") <- alpha
    matrixC
}

#' @rdname CalculateC
#' @param ... Arguments to CalculateCdirect
#' @param viaQR When TRUE QR is involved. This may be needed to handle colinear data. When NULL viaQR is set to TRUE if ordinary computations fail.
#' @param returnAlpha When TRUE alpha (1 or value below 1) is returned instead of C. Attribute viaQR is included.
#' @export
CalculateC <- function(a, b, ..., viaQR = NULL, returnAlpha = FALSE) {
  if (!is.null(viaQR)) {
    if (!viaQR) {
      matrixC <- CalculateCdirect(a, b, ...)
      if (returnAlpha) {
        matrixC <- attr(matrixC, "alpha")
      }
      attr(matrixC, "viaQR") <- viaQR
      return(matrixC)
    } else {
      abQR <- GenQR(rbind(a, b))
      inda <- seq_len(NROW(a))
      cQ <- CalculateCdirect(abQR$Q[inda, ], abQR$Q[-inda, ], ...)
      if (returnAlpha) {
        matrixC <- attr(cQ, "alpha")
      } else {
        matrixC <- cQ %*% abQR$R
        attr(matrixC, "alpha") <- attr(cQ, "alpha")
      }
      attr(matrixC, "viaQR") <- viaQR
      return(matrixC)
    }
  }
  z <- try(CalculateC(a, b, ..., viaQR = FALSE, returnAlpha = returnAlpha), silent = TRUE)
  if (!inherits(z, "try-error")) 
    return(z)
  z <- try(CalculateC(a, b, ..., viaQR = TRUE, returnAlpha = returnAlpha), silent = TRUE)
  if (!inherits(z, "try-error")) 
    return(z)
  if ("alpha" %in% names(list(...))) 
    stop(paste("Could not calculate C with alpha =", list(...)$alpha, "as input"))
  CalculateC(a, b, ..., viaQR = FALSE, alpha = 0, returnAlpha = returnAlpha)
}


#' Calculation of alpha
#' 
#' Function to find the largest alpha that makes equation 10 in the paper solvable.
#'
#' @param a matrix E in paper 
#' @param b matrix Eg in paper
#' @param tryViaQR When TRUE QR transformation used (to handle collinearity) 
#' when ordinary calculations fail.  
#'
#' @return alpha
#' 
#' @note FindAlphaSimple performs the calculations by a simple/direct method. 
#' FindAlpha is made to handle problematic special cases.
#'  
#' @export
#' @seealso See examples in the documentation of \code{\link{CalculateC}}
FindAlpha <- function(a, b, tryViaQR = TRUE) {
  ata <- t(a) %*% a
  btb <- t(b) %*% b
  m <- try(solve(ata, btb), silent = TRUE)
  if (!inherits(m, "try-error")) 
    return(sqrt(MinEigen(m, inverse = TRUE)))
  m <- try(solve(btb, ata), silent = TRUE)
  if (!inherits(m, "try-error")) 
    return(sqrt(MinEigen(m, inverse = FALSE)))
  if(tryViaQR){
    abQR = GenQR(rbind(a,b))
    inda = seq_len(NROW(a))
    alpha = FindAlpha(abQR$Q[inda,],abQR$Q[-inda,], tryViaQR = FALSE)
    attr(alpha, "viaQR") <- TRUE
    return(alpha)
  }
  warning("Could not calculate alpha. 0 returned")
  0
}


#' @rdname FindAlpha
#' @export
FindAlphaSimple <- function(a, b) {
  m <- t(a) %*% a %*% solve(t(b) %*% b)
  sqrt(min(eigen(m, only.values = TRUE)$values))
}


# Use of tol inspierd by ginv i MASS
MinEigen <- function(x, tol = sqrt(.Machine$double.eps), inverse = TRUE, ...) {
  a <- Re(eigen(x, symmetric = FALSE, only.values = TRUE, ...)$values)
  if (inverse) {
    ok <- abs(a) > max(tol * max(abs(a)), 0)
    a <- 1/a[ok]
  }
  max(0, min(a))
}

