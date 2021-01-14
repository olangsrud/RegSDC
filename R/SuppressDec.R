#' Suppressed tabular data: Inner cell frequencies as decimal numbers
#' 
#' Assume that frequencies to be published, \code{z}, can be computed from inner 
#' frequencies, \code{y}, via \code{ z = t(x) \%*\% y}, 
#' where \code{x} is a dummy matrix. 
#' Assuming correct suppression, this function will generate safe inner cell frequencies as decimal numbers.
#' 
#' This function makes use of \code{\link{ReduceX}} and \code{\link{RegSDCipso}}.
#' 
#'
#' @param x Dummy matrix where the dimensions matches z and/or y input. Sparse matrix (Matrix package) is possible.
#' @param z Frequencies to be published. All, only the safe ones or with suppressed as NA.
#' @param y Inner cell frequencies.
#' @param suppressed Logical vector defining the suppressed elements of z.
#' @param digits Output close to whole numbers will be rounded using \code{digits} as input to \code{\link{RoundWhole}}.
#' @param nRep Integer, when >1, several y's will be generated. Extra columns in output.
#' @param yDeduct Values to be subtracted from y and added back after the calculations. 
#'           Can be used to perform the modulo method described in the paper (see examples).
#' @param resScale Residuals will be scaled by resScale
#' @param rmse Desired root mean square error (residual standard error). Will be used when resScale is NULL or cannot be used.
#' 
#' @note Capital letters, X, Y and Z, are used in the paper.
#'
#' @return The inner cell frequencies as decimal numbers
#' @importFrom SSBtools RoundWhole
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' # Same data as in the paper
#' z <- RegSDCdata("sec7z")
#' x <- RegSDCdata("sec7x")
#' y <- RegSDCdata("sec7y")  # Now z is t(x) %*% y 
#' zAll <- RegSDCdata("sec7zAll")
#' zAllSupp <- RegSDCdata("sec7zAllSupp")
#' xAll <- RegSDCdata("sec7xAll")
#' 
#' # When no suppression, output is identical to y
#' SuppressDec(xAll, zAll, y)
#' SuppressDec(xAll, zAll)  # y can be seen in z
#' 
#' # Similar to Y* in paper (but other random values)
#' SuppressDec(x, z, y)
#' 
#' # Residual standard error forced to be 1
#' SuppressDec(x, z, y, rmse = 1)
#' 
#' # Seven ways of obtaining the same output
#' SuppressDec(x, z, rmse = 1)  # slower, y must be estimated
#' SuppressDec(x, y = y, rmse = 1)
#' SuppressDec(xAll, zAllSupp, y, rmse = 1)
#' SuppressDec(xAll, zAllSupp, rmse = 1)  # slower, y must be estimated
#' SuppressDec(xAll, zAll, y, is.na(zAllSupp), rmse = 1)
#' SuppressDec(xAll, zAll, suppressed = is.na(zAllSupp), rmse = 1)  # y seen in z
#' SuppressDec(xAll, y = y, suppressed = is.na(zAllSupp), rmse = 1)
#' 
#' # YhatMod4 and YhatMod10 in Table 2 in paper
#' SuppressDec(xAll, zAllSupp, y, yDeduct = 4 * (y%/%4), resScale = 0)
#' SuppressDec(xAll, zAllSupp, y, yDeduct = 10 * (y%/%10), rmse = 0)
#' 
#' # As data in Table 3 in paper (but other random values)
#' SuppressDec(xAll, zAllSupp, y, yDeduct = 10 * (y%/%10), resScale = 0.1)
#' 
#' # rmse instead of resScale and 5 draws
#' SuppressDec(xAll, zAllSupp, y, yDeduct = 10 * (y%/%10), rmse = 1, nRep = 5)
SuppressDec <- function(x, z = NULL, y = NULL, suppressed = NULL, digits = 9, nRep = 1, yDeduct = NULL, resScale = NULL, rmse = NULL) {
  origY <- !is.null(y)
  if (!is.null(z)) 
    z <- EnsureMatrix(z, NCOL(x))
  if (!is.null(y)) 
    y <- EnsureMatrix(y, NROW(x))
  
  
  if (!is.null(suppressed)) {
    suppr <- suppressed
  } else {
    
    if (is.null(z)) {
      suppr <- rep(FALSE, NCOL(x))
    } else {
      suppr <- is.na(z[, 1])
    }
  }
  
  if (is.null(y)) 
    if (!is.null(suppressed)) 
      if (!any(is.na(z))) 
        if (any(suppr)) {
          a <- ReduceX(x = x, z = z, digits = digits)
          y <- a$y
          origY <- !any(!a$yKnown)
        }
  
  nonSuppr <- which(!suppr)
  if (is.null(y)) {
    if (is.null(z)) {
      stop("z needed in input when y is NULL")
    }
  } else {
    z <- NULL
  }
  
  a <- ReduceX(x = x[, nonSuppr, drop = FALSE], z = z[nonSuppr, , drop = FALSE], y = y, digits = digits)
  
  
  if (is.null(y)) 
    origY <- !any(!a$yKnown)
  
  if (!origY) {
    if (is.null(rmse)) {
      stop("Without original y values, rmse must be supplied.")
    } else {
      if (!is.null(resScale)) {
        warning("Without original y values, rmse is used instead of resScale.")
      }
    }
  }
  
  
  if (any(!a$yKnown)){ 
    if(is.null(yDeduct)){
      rw <- RoundWhole(IpsoExtra(y = a$y[which(!a$yKnown), , drop = FALSE], x = a$x, nRep = nRep, 
                               ensureIntercept = FALSE, rmse = rmse, resScale = resScale), digits = digits)
    } else {
      yDeduct <- EnsureMatrix(yDeduct)[which(!a$yKnown), , drop = FALSE]
      rw <- RoundWhole(IpsoExtra(y = a$y[which(!a$yKnown), , drop = FALSE] - yDeduct, x = a$x, nRep = nRep, 
                                 ensureIntercept = FALSE, rmse = rmse, resScale = resScale), digits = digits)
      if (nRep != 1)
        yDeduct <- ColRepMatrix(yDeduct, nRep)
      rw <- rw + yDeduct
    }
  }
  
  if (nRep != 1) 
    a$y <- ColRepMatrix(a$y, nRep)
  if (any(!a$yKnown)) 
    a$y[which(!a$yKnown), ] <- rw
  a$y
}





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






#' Suppressed tabular data: Yhat from X and Z
#' 
#' Implementation of equation 21 in the paper.
#' 
#' Generalized inverse is computed by \code{\link{ginv}}.
#' In practise, the computations can be speeded up using reduced versions of X and Z. See \code{\link{ReduceX}}.
#'
#' @param z Z as a matrix
#' @param x X as a matrix
#' @param digits When non-NULL, output values close to whole numbers will be rounded using 
#'        \code{digits} as input to \code{\link{RoundWhole}}.
#'
#' @return Yhat as a matrix
#' 
#' @seealso \code{\link{IpsoExtra}}
#' 
#' @importFrom MASS ginv
#' @export
#' @author Øyvind Langsrud
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





#' Suppressed tabular data: Reduce dummy matrix, X (and estimate Y)
#' 
#' In section 7 in the paper \code{ Z = t(X) \%*\% Y} where \code{X} is a dummy matrix. 
#' Some elements of Y can be found directly as elements in Z. Corresponding rows of X will be removed. 
#' After removing rows, some columns will only have zeros and these will also be removed.
#' 
#' To estimate Y, this function finds some values directly from Z and other values by running \code{\link{Z2Yhat}} on reduced versions of X and Z.  
#' 
#'
#' @param z Z as a matrix
#' @param x X as a matrix
#' @param y Y as a matrix
#' @param digits When non-NULL and when NULL y input, output y estimates close to whole numbers will be rounded using 
#'        \code{digits} as input to \code{\link{RoundWhole}}.
#'
#' @return A list of four elements:
#'         \item{\code{x}}{Reduced \code{x}}
#'         \item{\code{z}}{Corresponding reduced \code{z} or NULL when no \code{z} in input}
#'         \item{\code{yKnown}}{Logical vector specifying elements of y that can be found directly as elements in z}
#'         \item{\code{y}}{As \code{y} in input (not reduced) or estimated \code{y} when NULL y in input}
#'         
#' @keywords internal
#' @importFrom  methods as
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' # Same data as in the paper
#' z <- RegSDCdata("sec7z")
#' x <- RegSDCdata("sec7x")
#' y <- RegSDCdata("sec7y")  # Now z is t(x) %*% y 
#' 
#' a <- ReduceX(x, z, y)
#' b <- ReduceX(x, z)
#' d <- ReduceX(x, z = NULL, y)  # No z in output
#' 
#' # Identical output for x and z
#' identical(a$x, b$x)
#' identical(a$x, d$x)
#' identical(a$z, b$z)
#' 
#' # Same y in output as input
#' identical(a$y, y)
#' identical(d$y, y)
#' 
#' # Estimate of y (yHat) when NULL y input
#' b$y
#' 
#' # These elements of y can be found directly in in z
#' y[a$yKnown, , drop = FALSE]
#' # They can be found by searching for unit colSums
#' colSums(x)[colSums(x) == 1]
#' 
#' # These trivial data rows can be omitted when processing data
#' x[!a$yKnown, ]
#' # Now several columns can be omitted since zero colSums
#' colSums0 <- colSums(x[!a$yKnown, ]) == 0
#' # The resulting matrix is output from the function
#' identical(x[!a$yKnown, !colSums0], a$x)
#' 
#' # Output z can be computed from this output x
#' identical(t(a$x) %*% y[!a$yKnown, , drop = FALSE], a$z)
ReduceX <- function(x, z = NULL, y = NULL, digits = 9) {
  yNULL <- is.null(y)
  zNULL <- is.null(z)
  if(yNULL & zNULL)
    stop("z or y must be supplied")
  colSums_1 <- which(colSums(x) == 1)
  x1 <- x[, colSums_1, drop = FALSE]
  x1dgT <- as(x1, "dgTMatrix")
  nonDub <- x1dgT@j[x1dgT@x != 0][!duplicated(x1dgT@i[x1dgT@x != 0])] + 1L
  x1 <- x1[, nonDub, drop = FALSE]
  if (!zNULL) 
    zA <- z[colSums_1[nonDub], , drop = FALSE]
  zA1 <- matrix(1, NCOL(x1), 1)
  yKnown1 <- round(x1 %*% zA1)
  yKnown1_0 <- which(yKnown1 == 0)
  if (yNULL) {
    yHat <- x1 %*% zA
  } else {
    yHat <- y
    yHat[yKnown1_0, ] <- 0
  }
  if (!zNULL) 
    z <- z - crossprod(x, yHat)
  x <- x[yKnown1_0, , drop = FALSE]
  colSums_ok <- which(colSums(x) != 0)
  if (!zNULL) 
    z <- z[colSums_ok, , drop = FALSE]
  x <- x[, colSums_ok, drop = FALSE]
  if (yNULL) {
    if (length(yKnown1_0))
      yHat[yKnown1_0, ] <- Z2Yhat(z, x, digits = NA)
    if (!is.null(digits))
      if (!is.na(digits)) 
        yHat <- RoundWhole(yHat, digits = digits)
  } else {
    yHat <- y
  }
  list(x = x, z = z, yKnown = yKnown1 != 0, y = yHat)
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
#'          NULL or cannot be used (see parameter digits). This parameter forces the rmse value for one y variable (the first). 
#'
#' @return Generated version of y
#' @export
#' @author Øyvind Langsrud
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
  
  #if (!is.null(rmse)) 
  #  if (NCOL(y) > 1) 
  #    stop("rmse parameter only when single y")
  
  if (NROW(xQ) == NCOL(xQ)) {
    if (!is.null(resScale) | !is.null(rmse)) 
      warning("resScale/rmse ignored when Q from X is square.")
    if (nRep != 1) 
      y <- ColRepMatrix(y, nRep)
    if (returnParts) 
      return(list(yHat = y, yRes = 0 * y)) else return(y)
  }
  
  yHat <- xQ %*% (t(xQ) %*% y)
  
  n <- NROW(y)
  ncoly <- NCOL(y)
  
  eQRR <- NULL
  if (!is.null(digits) & !is.null(resScale)) {
    if (!is.null(rmse)) 
      if (max(abs(round(y - yHat, digits = digits))) == 0) {
        if (ncoly > 1){
          warning("rmse with identical residual vectors used instead of resScal since perfect fit.")
        } else {
          warning("rmse used instead of resScal since perfect fit.")
        }
        resScale <- NULL
        eQRR <- matrix(1, 1, ncoly)  # Changed below
        m <- 1L
      }
  }
  
  if (is.null(eQRR)) {
    eQRR <- GenQR(y - yHat, makeunique = TRUE)$R
    m <- NROW(eQRR)
  }
  
  
  if (!is.null(rmse)) 
    if (ncoly > 1){ 
      rmseVar <- match(TRUE,!is.na(rmse))
      rmse <- rmse[rmseVar] 
      resScale <- rmse * sqrt((n - NCOL(xQ))/sum(eQRR[, rmseVar]^2)) 
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
    else {
      yRes[, SeqInc(1 + ncoly * (i - 1), ncoly * i)] <- eSimQ %*% eQRR
    }
  }
  if(!is.null(rownames(y))){
    rownames(yHat) <- rownames(y)
    rownames(yRes) <- rownames(y)
  }
  if (returnParts) 
    return(list(yHat = yHat, yRes = yRes))
  yHat + yRes
}