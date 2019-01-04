#' Regression-based SDC Tools - Component score correlation control
#' 
#' Implementation of equation 8 in the paper.
#' 
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param compCorr Required component score  correlations (possibly recycled)
#' @param x Matrix of non-confidential variables
#' @param doSVD SVD when TRUE and QR when FALSE
#' @param makeunique Parameter to be used in GenQR
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' 
#' @details NA component score correlation means independent random. Input matrices are subjected to \code{\link{EnsureMatrix}}.
#' 
#'
#' @return Generated version of y
#' @export
#'
#' @examples
#' x <- matrix(1:10, 10, 1)
#' y <- matrix(rnorm(30) + 1:30, 10, 3)
#' 
#' # Same as IPSO (RegSDCipso)
#' RegSDCcomp(y, NA, x)
#' 
#' # Using QR and SVD
#' yQR <- RegSDCcomp(y, c(0.1, 0.2, NA), x)
#' ySVD <- RegSDCcomp(y, c(0.1, 0.2, NA), x, doSVD = TRUE)
#' 
#' # Calculation of residuals
#' r <- residuals(lm(y ~ x))
#' rQR <- residuals(lm(yQR ~ x))
#' rSVD <- residuals(lm(ySVD ~ x))
#' 
#' # Correlations for two first components as required
#' diag(cor(GenQR(r)$Q, GenQR(rQR)$Q))
#' diag(cor(GenQR(r, doSVD = TRUE)$Q, GenQR(rSVD, doSVD = TRUE)$Q))
#' 
#' # Identical covariance matrices
#' cov(yQR) - cov(ySVD)
#' cov(rQR) - cov(rSVD)
#' 
#' # Identical regression results
#' summary(lm(y[, 1] ~ x))
#' summary(lm(yQR[, 1] ~ x))
#' summary(lm(ySVD[, 1] ~ x))
RegSDCcomp <- function(y, compCorr = NA, x = NULL, doSVD = FALSE, makeunique = TRUE, ensureIntercept = TRUE) {
  y <- EnsureMatrix(y)
  x <- EnsureMatrix(x, nrow(y))
  if(ensureIntercept)
    x <- EnsureIntercept(x)
  xQ <- GenQR(x, doSVD = doSVD, findR = FALSE)
  yHat <- xQ %*% (t(xQ) %*% y)
  eQR <- GenQR(y - yHat, doSVD = doSVD, makeunique = makeunique)
  n <- NROW(y)
  m <- NCOL(eQR$Q)
  compCorr <- rep_len(compCorr, m)
  
  colRotate <- is.na(compCorr)
  colFixed <- compCorr == 1
  colFixed[is.na(colFixed)] <- FALSE
  colCorr <- (!colRotate) & (!colFixed)
  
  m_ <- m - sum(colFixed)
  randData <- matrix(rnorm(n * m_), n, m_)
  
  nFixed <- sum(colFixed)
  nCorr <- sum(colCorr)
  nRotate <- sum(colRotate)
  nxQ <- NCOL(xQ)
  
  newQ <- GenQR(cbind(xQ, eQR$Q[, colFixed, drop = FALSE], eQR$Q[, colCorr, drop = FALSE], randData), findR = FALSE, makeunique = makeunique)
  if ((nxQ + nFixed + nCorr + m_) > NCOL(newQ)) 
    stop("Not enough dimensions")
  newQ <- newQ[, matlabColon(nxQ + nFixed + nCorr + 1, nxQ + nFixed + nCorr + m_), drop = FALSE]
  
  cCorr <- compCorr[colCorr]
  
  eSimQ <- eQR$Q
  
  eSimQ[, colRotate] <- newQ[, matlabColon(1, nRotate), drop = FALSE]
  
  eSimQ[, colCorr] <- t(cCorr * t(eSimQ[, colCorr, drop = FALSE]) + sqrt((1 - cCorr^2)) * t(newQ[, matlabColon(nRotate + 1, m_), drop = FALSE]))
  
  yHat + eSimQ %*% eQR$R
}
