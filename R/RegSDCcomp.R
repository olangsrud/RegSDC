#' Regression-based SDC Tools - Component score correlation control
#' 
#' Implementation of equation 8 in the paper.
#' 
#' @encoding UTF8
#'
#' @param y Matrix of confidential variables
#' @param compCorr Required component score  correlations (possibly recycled)
#' @param x Matrix of non-confidential variables (including a constant term (column of ones)) 
#' @param doSVD SVD when TRUE and QR when FALSE
#' @param makeunique Parameter to be used in GenQR
#' 
#' @details NA component score  correlation means independent random
#'
#' @return Generated version of y
#' @export
#'
#' @examples
#' x1 <- 1:10
#' x <- cbind(x0 = 1, x1 = x1)
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
#' r <- residuals(lm(y ~ x1))
#' rQR <- residuals(lm(yQR ~ x1))
#' rSVD <- residuals(lm(ySVD ~ x1))
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
#' summary(lm(y[, 1] ~ x1))
#' summary(lm(yQR[, 1] ~ x1))
#' summary(lm(ySVD[, 1] ~ x1))
RegSDCcomp <- function(y, compCorr = NA, x = matrix(1, NROW(y), 1), doSVD = FALSE, makeunique = TRUE) {
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
