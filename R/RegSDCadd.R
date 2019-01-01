#' Regression-based SDC Tools - Synthetic addition
#' 
#' Residuals from arbitrary data with a synthetic addition
#'
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param yStart Arbitrary data whose residuals will be used
#' @param x Matrix of non-confidential variables (including a constant term (column of ones)) 
#' @param epsAlpha Precision constant for alpha calculation 
#' @param AlphaHandler Function (warning or stop) to be used when alpha<1 
#' @param alphaAttr When TRUE alpha is attribute in output 
#' @param makeunique Parameter to be used in GenQR
#' 
#' @details  Use epsAlpha=NULL to avoid calculation of alpha. Use of alpha (<1) will produce a warning. 
#' 
#' @return Generated version of y
#' @keywords internal
#' @export
RegSDCaddGen <- function(y, yStart, x = matrix(1, NROW(y), 1), epsAlpha = 1e-07, AlphaHandler = warning, alphaAttr = TRUE, makeunique = TRUE) {
  xQ <- GenQR(x, doSVD = FALSE, findR = FALSE)
  yHat <- xQ %*% (t(xQ) %*% y)
  yHatStart <- xQ %*% (t(xQ) %*% yStart)
  A <- y - yHat
  B <- yStart - yHatStart
  
  R <- CalculateC(A, B, epsAlpha = epsAlpha, AlphaHandler = AlphaHandler, alpha = NULL)
  alpha <- attr(R, "alpha")
  
  n <- NROW(y)
  m <- NROW(R)
  ySim <- matrix(rnorm(n * m), n, m)
  
  newQ <- GenQR(cbind(xQ, yStart, ySim), findR = FALSE, makeunique = makeunique)
  
  nColXQyStart <- NCOL(xQ) + NCOL(yStart)
  if (NCOL(newQ) < (nColXQyStart + m)) 
    nColXQyStart <- qr(cbind(xQ, yStart), tol = 1e-07)$rank
  
  newQ <- newQ[, -seq_len(nColXQyStart), drop = FALSE]
  
  if (NCOL(newQ) < m) 
    stop("Not enough dimensions")
  
  yOut <- yHat + alpha * B + newQ %*% R
  if (alphaAttr) 
    attr(yOut, "alpha") <- alpha
  yOut
}



#' Regression-based SDC Tools - Synthetic addition with residual correlation control
#' 
#' Implementation of equation 6 (arbitrary residual data) and equation 7 (residual correlations) in the paper.
#' The alpha limit is calculated (equation 9). The limit is used when alpha =1 cannot be chosen (warning produced). 
#' In output, alpha is attribute.
#'
#'
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param resCorr  Required residual correlations (possibly recycled)
#' @param x Matrix of non-confidential variables (including a constant term (column of ones)) 
#' @param yStart Arbitrary data whose residuals will be used. Will be calculated from resCorr when NULL.
#' 
#' @details  Use epsAlpha=NULL to avoid calculation of alpha. Use of alpha (<1) will produce a warning. 
#' 
#' @return Generated version of y with alpha as attribute
#' @export
#' @examples
#' x1 <- 1:10
#' x <- cbind(x0 = 1, x1 = x1)
#' y <- matrix(rnorm(30) + 1:30, 10, 3)
#' yOut <- RegSDCadd(y, c(0.1, 0.2, 0.3), x)
#' 
#' # Correlations between residuals as required
#' diag(cor(residuals(lm(y ~ x1)), residuals(lm(yOut ~ x1))))
#' 
#' # Identical covariance matrices
#' cov(y) - cov(yOut)
#' cov(residuals(lm(y ~ x1))) - cov(residuals(lm(yOut ~ x1)))
#' 
#' # Identical regression results
#' summary(lm(y[, 1] ~ x1))
#' summary(lm(yOut[, 1] ~ x1))
#' 
#' # alpha as attribute
#' attr(yOut, "alpha")
#' 
#' # With yStart as input and alpha limit in use (warning produced)
#' yOut <- RegSDCadd(y, NULL, x, 2 * y + matrix(rnorm(30), 10, 3))
#' attr(yOut, "alpha")
#' 
#' # Same correlation for all variables
#' RegSDCadd(y, 0.2, x)
#' # But in this case RegSDCcomp is equivalent and faster
#' RegSDCcomp(y, 0.2, x)
#' 
#' 
#' # Make nearly collinear data
#' y[, 3] <- y[, 1] + y[, 2] + 0.001 * y[, 3]
#' # Not possible to achieve correlations. Small alpha with warning.
#' RegSDCadd(y, c(0.1, 0.2, 0.3), x)
#' # Exact collinear data
#' y[, 3] <- y[, 1] + y[, 2]
#' # Zero alpha with warning
#' RegSDCadd(y, c(0.1, 0.2, 0.3), x)
RegSDCadd <- function(y, resCorr = NULL, x = matrix(1, NROW(y), 1), yStart = NULL) {
  if (is.null(resCorr)) {
    if (is.null(yStart)) 
      stop("resCorr or yStart must be supplied")
  } else {
    if (!is.null(yStart)) 
      stop("resCorr or yStart must be NULL")
    resCorr <- rep_len(resCorr, NCOL(y))
    yStart <- outer(rep(1, NROW(y)), resCorr) * y
  }
  RegSDCaddGen(y, yStart, x)
}




