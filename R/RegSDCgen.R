#' Regression-based SDC Tools - General data generation
#'
#' IPSO by QR or SVD, scores from arbitrary data, and ROMM
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param x Matrix of non-confidential variables
#' @param doSVD SVD when TRUE and QR when FALSE
#' @param yNew Matrix of y-data for new scores (simulated when NULL)
#' @param lambda ROMM parameter
#' @param makeunique Parameter to be used in GenQR 
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' @param returnParts Alternative output two matrices: yHat (fitted) and yRes (generated residuals).
#' 
#' @details doSVD has effect on decomposition of y and yNew. Input matrices are subjected to \code{\link{EnsureMatrix}}.
#' 
#' @return Generated version of y
#' @keywords internal
#' @importFrom stats rnorm
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' exY <- matrix(rnorm(15), 5, 3)
#' RegSDCgen(exY)
#' RegSDCgen(exY, yNew = exY + 0.001 * matrix(rnorm(15), 5, 3))  # Close to exY
#' RegSDCgen(exY, lambda = 0.001)  # Close to exY
RegSDCgen <- function(y, x = NULL, doSVD = FALSE, yNew = NULL, lambda = Inf, makeunique = TRUE, 
                      ensureIntercept = TRUE, returnParts = FALSE) {
  y <- EnsureMatrix(y)
  x <- EnsureMatrix(x, nrow(y))
  if(ensureIntercept)
    x <- EnsureIntercept(x)
  xQ <- GenQR(x, findR = FALSE)
  
  if (NROW(xQ) == NCOL(xQ)) {
    if (returnParts) 
      return(list(yHat = y, yRes = 0 * y)) 
    else 
      return(y)
  }
  
  yHat <- xQ %*% (t(xQ) %*% y)
  # makeunique <- !is.null(yNew)
  eQR <- GenQR(y - yHat, doSVD = doSVD, makeunique = makeunique)
  n <- NROW(y)
  m <- NCOL(eQR$Q)
  if (is.null(yNew)) {
    yNew <- matrix(rnorm(n * m), n, m)
    doSVD <- FALSE
  } else
    yNew <- EnsureMatrix(yNew, nrow(y)) # ncol is not checked and could be different 
  if (is.finite(lambda)) {
    yNew <- eQR$Q + lambda * yNew
    doSVD <- FALSE
  }
  eSim <- yNew - xQ %*% (t(xQ) %*% yNew)
  eSimQ <- GenQR(eSim, doSVD = doSVD, findR = FALSE, makeunique = makeunique)
  if (NCOL(eSimQ) > m) 
    eSimQ <- eSimQ[, seq_len(m), drop = FALSE]
  if(!is.null(rownames(y))){
    rownames(yHat) <- rownames(y)
    rownames(eSimQ) <- rownames(y)
  }
  if (returnParts)
    return(list(yHat = yHat, yRes = eSimQ %*% eQR$R))
  yHat + eSimQ %*% eQR$R
}



#' Regression-based SDC Tools - Ordinary synthetic data (IPSO)
#'
#' Implementation of equation 4 in the paper.
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param x Matrix of non-confidential variables
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' 
#' @details Input matrices are subjected to \code{\link{EnsureMatrix}}.
#' 
#' @return Generated version of y
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' x <- matrix(1:5, 5, 1)
#' y <- matrix(rnorm(15) + 1:15, 5, 3)
#' ySynth <- RegSDCipso(y, x)
#' 
#' # Identical regression results
#' summary(lm(y[, 1] ~ x))
#' summary(lm(ySynth[, 1] ~ x))
#' 
#' # Identical covariance matrices
#' cov(y) - cov(ySynth)
#' cov(residuals(lm(y ~ x))) - cov(residuals(lm(ySynth ~ x)))
RegSDCipso <- function(y, x = NULL, ensureIntercept = TRUE) {
  RegSDCgen(y = y, x = x, ensureIntercept = ensureIntercept)
}


#' Regression-based SDC Tools - Scores from new data
#' 
#' Implementation of equation 12 in the paper.
#'
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param yNew Matrix of y-data for new scores 
#' @param x Matrix of non-confidential variables
#' @param doSVD SVD when TRUE and QR when FALSE
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' 
#' @details doSVD has effect on decomposition of y and yNew. Input matrices are subjected to \code{\link{EnsureMatrix}}.  
#' 
#' @return Generated version of y
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' x <- matrix(1:5, 5, 1)
#' y <- matrix(rnorm(15) + 1:15, 5, 3)
#' 
#' # Same as IPSO (RegSDCipso)
#' RegSDCnew(y, matrix(rnorm(15), 5, 3), x)
#' 
#' # Close to y
#' RegSDCnew(y, y + 0.001 * matrix(rnorm(15), 5, 3), x)
RegSDCnew <- function(y, yNew, x = NULL, doSVD = FALSE, ensureIntercept = TRUE) {
  RegSDCgen(y = y, x = x, doSVD =doSVD, yNew = yNew, ensureIntercept = ensureIntercept)
}


#' Regression-based SDC Tools - Random orthogonal matrix masking (ROMM)
#'
#' Implementation based on equations 11, 12 and 17 in the paper.
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param lambda ROMM parameter
#' @param x Matrix of non-confidential variables
#' @param doSVD SVD when TRUE and QR when FALSE
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' 
#' @details doSVD has effect on decomposition of y. 
#' The exact behaviour of the method depends on the choice of the decomposition method because of 
#' the sequentially phenomenon mentioned in the paper. 
#' The similarity to the original data will tend to be highest for the first component. 
#' Input matrices are subjected to \code{\link{EnsureMatrix}}.
#' 
#' @return Generated version of y
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' x <- matrix(1:5, 5, 1)
#' y <- matrix(rnorm(15) + 1:15, 5, 3)
#' 
#' # Same as IPSO (RegSDCipso)
#' RegSDCromm(y, Inf, x)
#' 
#' # Close to IPSO
#' RegSDCromm(y, 100, x)
#' 
#' # Close to y
#' RegSDCromm(y, 0.001, x)
RegSDCromm <- function(y, lambda = Inf, x = NULL, doSVD = FALSE, ensureIntercept = TRUE) {
  RegSDCgen(y = y, x = x, doSVD =doSVD, lambda = lambda, ensureIntercept = ensureIntercept)
}

