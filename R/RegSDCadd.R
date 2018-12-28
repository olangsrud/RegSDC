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
#' 
#' @details  Use epsAlpha=NULL to avoid calculation of alpha. Use of alpha (<1) will produce a warning. 
#' 
#' @return Generated version of y
#' @keywords internal
#' @export
RegSDCaddGen <- function(y, yStart, x = matrix(1, NROW(y), 1), epsAlpha = 1e-07, AlphaHandler = warning, alphaAttr = TRUE) {
  xQ <- GenQR(x, doSVD = FALSE, findR = FALSE)
  yHat <- xQ %*% (t(xQ) %*% y)
  yHatStart <- xQ %*% (t(xQ) %*% yStart)
  A <- y - yHat
  B <- yStart - yHatStart
  if (is.null(epsAlpha)) 
    alpha <- 1 else {
      Mf <- t(A) %*% A %*% solve(t(B) %*% B)
      alpha <- min(1, sqrt(min(eigen(Mf, only.values = TRUE)$values))/(1 + epsAlpha))  #
      if (alpha < 1) 
        AlphaHandler(paste("alpha = ", alpha))
    }
  R <- chol(t(A) %*% A - alpha^2 * t(B) %*% B)  # collinear y not treated
  n <- NROW(y)
  m <- NCOL(y)  # collinear y not treated
  ySim <- matrix(rnorm(n * m), n, m)
  
  newQ <- GenQR(cbind(xQ, yStart, ySim), findR = FALSE)[, -seq_len(NCOL(xQ) + m), drop = FALSE]
  
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




