

ColRepMatrix <- function(x, nRep) {
  if (nRep == 1) 
    return(x)
  if (nRep == 0) 
    return(x[, integer(0), drop = FALSE])
  cn <- colnames(x)
  rn <- rownames(x)
  x <- matrix(x, nrow(x), ncol(x) * nRep)
  if (!is.null(cn)) 
    colnames(x) <- rep_len(cn, ncol(x))
  if (!is.null(rn)) 
    rownames(x) <- rn
  x
}


# Own function in SSBtools
SeqInc <- function (from, to)  
{
  if (from > to){
    if(from-to> 1L)
      stop("Length of sequence (1+to-from) must be non-negative")
    integer(0)
  }
  else from:to
}






#' Yhat from X and Z
#' 
#' Implementation of equation 21 in the paper.
#' 
#' Generalized inverse is computed by \code{\link{ginv}}.
#' In practise, the computations can be speeded up using reduced versions of X and Z. See.
#'
#' @param z Z as a matrix
#' @param x X as a matrix
#' @param digits When non-NULL, output values close to whole numbers will be rounded using 
#'        \code{digits} as input to \code{\link{RoundWhole}}.
#'
#' @return Yhat as a matrix
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' # Same data as in the paper
#' z <- RegSDCdata("sec7z")
#' x <- RegSDCdata("sec7x")
#' Z2Yhat(z, x)
#' 
#' # With y known, yHat can be computed in other ways
#' y <- RegSDCdata("sec7y")  # Now z is t(x) %*% y 
#' fitted(lm(y ~ x - 1))
#' IpsoExtra(y, x, FALSE, resScale = 0)
Z2Yhat <- function(z, x, digits = 9) {
  yHat <- crossprod(ginv(as.matrix(x)), z)
  if(!is.null(rownames(x)))
    rownames(yHat) <- rownames(x)
  if (!is.null(digits)) 
    if (!is.na(digits)) 
      yHat <- RoundWhole(yHat, digits = digits)
  yHat
}




#' Round values that are close two whole numbers
#'
#' @param x vector or matrix
#' @param digits parameter to \code{\link{round}}
#' @param onlyZeros Only round values close to zero 
#'
#' @return
#' @export
#' @keywords internal 
#'
#' @examples
#' x <- c(0.0002, 1.00003, 3.00014)
#' RoundWhole(x)     # No values rounded
#' RoundWhole(x, 4)  # One value rounded
#' RoundWhole(x, 3)  # All values rounded
#' RoundWhole(x, 3, TRUE)  # One value rounded
RoundWhole <- function(x, digits = 9, onlyZeros = FALSE) {
  round_x <- round(x)
  round_x_digits <- round(x, digits = digits)
  if (onlyZeros) 
    toWhole <- round_x_digits == 0 
  else 
    toWhole <- round_x == round_x_digits
  x[toWhole] <- round_x[toWhole]
  x
}




#' Extended variant of RegSDCipso
#' 
#' Possible to generate several y's and to re-scale residuals.
#' 
#' @param y Matrix of confidential variables
#' @param x Matrix of non-confidential variables
#' @param ensureIntercept Whether to ensure/include a constant term. Non-NULL x is subjected to \code{\link{EnsureIntercept}}
#' @param returnParts Alternative output two matrices: yHat (fitted) and yRes (generated residuals).
#' @param nRep Integer, when >1, several y's will be generated. Extra columns in output.
#' @param resScale Residuals will be scaled by resScale
#' @param digits Digits used to detect perfect fit (caused by fitted values as input). 
#'      This checking will be done only when rmse is in input. When perfect fit, rmse will be used instead of resScale.
#' @param rmse Desired root mean square error (residual standard error). Will be used when resScale is 
#'          NULL or cannot be used (see parameter digits). This parameter is possible only when single y variable. 
#'
#' @return Generated version of y
#' @export
#' @keywords internal 
#'
#' @examples
#' x <- matrix(1:5, 5, 1)
#' y <- matrix(10 * (sample(7:39, 15) + 4 * (1:15)), 5, 3)
#' colnames(y) <- paste("y", 1:3, sep = "")
#' y1 <- y[, 1, drop = FALSE]
#' 
#' IpsoExtra(y, x)  # Same as RegSDCipso(y, x)
#' 
#' IpsoExtra(y, x, resScale = 0)  # Fitted values (whole numbers in this case)
#' IpsoExtra(y, x, nRep = 2, resScale = 1e-05)  # Downscaled residuals 
#' 
#' ySynth <- IpsoExtra(y1, x, nRep = 2, rmse = 0.25)  # Downscaled residuals 
#' summary(lm(ySynth ~ x))  # Identical regression results with Residual standard error: 0.25
#' 
#' IpsoExtra(fitted(lm(y1 ~ x)), x, nRep = 2, resScale = 0.1)  # resScale no effect since perfect fit
#' IpsoExtra(fitted(lm(y1 ~ x)), x, nRep = 2, resScale = 0.1, rmse = 2)  # with warning
#' 
#' # Using data in the paper
#' IpsoExtra(RegSDCdata("sec7y"), RegSDCdata("sec7x"))  # Similar to Y*
#' IpsoExtra(RegSDCdata("sec7y"), RegSDCdata("sec7x"), rmse = 1)
IpsoExtra <- function(y, x = NULL, ensureIntercept = TRUE, returnParts = FALSE, nRep = 1, resScale = NULL, digits = 9, rmse = NULL) {
  y <- EnsureMatrix(y)
  x <- EnsureMatrix(x, nrow(y))
  if (ensureIntercept) 
    x <- EnsureIntercept(x)
  xQ <- GenQR(x, findR = FALSE)
  
  if (!is.null(rmse)) 
    if (NCOL(y) > 1) 
      stop("rmse parameter only when single y")
  
  if (NROW(xQ) == NCOL(xQ)) {
    if (!resScale | !rmse) 
      warning("resScale/rmse ignored when Q from X is square.")
    if (nRep != 1) 
      y <- ColRepMatrix(y, nRep)
    if (returnParts) 
      return(list(yHat = y, yRes = 0 * y)) else return(y)
  }
  
  yHat <- xQ %*% (t(xQ) %*% y)
  
  n <- NROW(y)
  
  eQRR <- NULL
  if (!is.null(digits) & !is.null(resScale)) {
    if (!is.null(rmse)) 
      if (max(abs(round(y - yHat, digits = digits))) == 0) {
        warning("rmse used instead of resScal since perfect fit.")
        resScale <- NULL
        eQRR <- matrix(1, 1, 1)  # Changed below
        m <- 1L
      }
  }
  
  if (is.null(eQRR)) {
    eQRR <- GenQR(y - yHat, makeunique = TRUE)$R
    m <- NROW(eQRR)
  }
  if (!is.null(resScale)) {
    eQRR <- resScale * eQRR
  } else {
    if (!is.null(rmse)) 
      eQRR[] <- sqrt(n - NCOL(xQ)) * rmse
  }
  
  if (nRep != 1) {
    yHat <- ColRepMatrix(yHat, nRep)
    yRes <- 0 * yHat
  }
  for (i in seq_len(nRep)) {
    yNew <- matrix(rnorm(n * m), n, m)
    eSim <- yNew - xQ %*% (t(xQ) %*% yNew)
    eSimQ <- GenQR(eSim, findR = FALSE, makeunique = TRUE)
    if (nRep == 1) 
      yRes <- eSimQ %*% eQRR 
    else 
      yRes[, SeqInc(1 + m * (i - 1), m * i)] <- eSimQ %*% eQRR
  }
  if(!is.null(rownames(y))){
    rownames(yHat) <- rownames(y)
    rownames(yRes) <- rownames(y)
  }
  if (returnParts) 
    return(list(yHat = yHat, yRes = yRes))
  yHat + yRes
}