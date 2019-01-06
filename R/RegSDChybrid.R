
  
GroupModelMatrix <- function(group, x = NULL) {
  if (is.null(group)) 
    return(x)
  xGroup <- droplevels(as.factor(group))
  
  if (is.null(x)) {
    if (length(levels(xGroup)) <= 1) 
      m <- matrix(1, length(xGroup), 1) else m <- model.matrix(~xGroup)[, ]
      return(m)
  }
  if (length(group) != NROW(x)) 
    stop("length(group) != NROW(x)")
  if (length(levels(xGroup)) <= 1 | NCOL(x) < 1) 
    return(x)
  model.matrix(~xGroup:x - 1)[, ]
}


AdjustQR <- function(y, x = NULL, xOrth = GenQR(x, findR = FALSE), ...) {
  n <- ncol(xOrth)
  z <- GenQR(cbind(xOrth, y), ...)
  if (is.null(n)) 
    return(z)
  if (n == 0) 
    return(z)
  if (is.list(z)) {
    z$Q <- z$Q[, -seq_len(n), drop = FALSE]
    z$R <- z$R[-seq_len(n), -seq_len(n), drop = FALSE]
    return(z)
  }
  z[, -seq_len(n), drop = FALSE]
}






#' Regression-based SDC Tools - Generalized microaggregation
#' 
#' Implementation of the methodology in section 6 in the paper
#'
#'
#' @encoding UTF8
#' @param y Matrix of confidential variables
#' @param clusters Vector of cluster coding
#' @param xLocal Matrix of x-variables to be crossed with clusters
#' @param xGlobal Matrix of x-variables NOT to be crossed with clusters
#' @param clusterPieces Vector of coding of cluster pieces
#' @param xClusterPieces Matrix of x-variables to be crossed with cluster pieces
#' @param groupedClusters Vector of coding of grouped clusters
#' @param xGroupedClusters Matrix of x-variables to be crossed with grouped clusters
#' @param alternative One of "" (default), "a", "b" or "c" 
#' @param alpha Possible to specify parameter used internally by alternative "c"
#' @param ySim  Possible to specify the internally simulated data manually
#' @param returnParts Alternative output six matrices:  
#'        y1 and y2 (fitted),    e3s and e4s (new residuals),   e3 and e4 (original residuals) 
#' @param epsAlpha Precision constant for alpha calculation
#' @param makeunique  Parameter to be used in GenQR
#'
#' 
#' @return Generated version of y
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' y <- matrix(rnorm(30) + 1:30, 10, 3)
#' x <- matrix(1:10, 10, 1)
#' 
#' # Same as RegSDCipso(y)
#' yOut <- RegSDChybrid(y)
#' 
#' # With a single cluster both are same as RegSDCipso(y, x)
#' yOut <- RegSDChybrid(y, xLocal = x)
#' yOut <- RegSDChybrid(y, xGlobal = x)
#' 
#' # Define two clusters
#' clust <- rep(1:2, each = 5)
#' 
#' # MHa and MHb in paper
#' yMHa <- RegSDChybrid(y, clusters = clust, xLocal = x)
#' yMHb <- RegSDChybrid(y, clusterPieces = clust, xLocal = x)
#' 
#' # An extended variant of MHb as mentioned in paper paragraph below definition of MHa/MHb
#' yMHbExt <- RegSDChybrid(y, clusterPieces = clust, xClusterPieces = x)
#' 
#' # Identical means within clusters
#' aggregate(y, list(clust = clust), mean)
#' aggregate(yMHa, list(clust = clust), mean)
#' aggregate(yMHb, list(clust = clust), mean)
#' aggregate(yMHbExt, list(clust = clust), mean)
#' 
#' # Identical global regression results
#' summary(lm(y[, 1] ~ x))
#' summary(lm(yMHa[, 1] ~ x))
#' summary(lm(yMHb[, 1] ~ x))
#' summary(lm(yMHbExt[, 1] ~ x))
#' 
#' # MHa: Identical local regression results
#' summary(lm(y[, 1] ~ x, subset = clust == 1))
#' summary(lm(yMHa[, 1] ~ x, subset = clust == 1))
#' 
#' # MHb: Different results
#' summary(lm(yMHb[, 1] ~ x, subset = clust == 1))
#' 
#' # MHbExt: Same estimates and different std. errors
#' summary(lm(yMHbExt[, 1] ~ x, subset = clust == 1))
#' 
#' # More examples will be added
RegSDChybrid <- function(y, clusters = NULL, xLocal = NULL, xGlobal = NULL, clusterPieces = NULL, 
                         xClusterPieces = NULL, groupedClusters = NULL, xGroupedClusters = NULL, 
                         alternative = NULL, alpha = NULL, ySim = NULL, returnParts = FALSE, 
                         epsAlpha = 1e-07, makeunique = TRUE) {
  y <- EnsureMatrix(y)
  n <- nrow(y)
  m <- ncol(y)
  if (is.null(clusters)) 
    clusters <- rep(1, n)
  xLocal <- EnsureMatrix(xLocal, n)
  xGlobal <- EnsureMatrix(xGlobal, n)
  xClusterPieces <- EnsureMatrix(xClusterPieces, n)
  xClusterPieces <- EnsureIntercept(xClusterPieces)
  xGroupedClusters <- EnsureMatrix(xGroupedClusters, n)
  if (is.null(alternative)) 
    alternative <- ""
  
  if (!is.null(groupedClusters)) {
    xGlobal <- cbind(GroupModelMatrix(groupedClusters, xGroupedClusters), xGlobal)
  }
  
  xGlobalOrth <- matrix(NaN, n, NCOL(xGlobal))
  yHat1 <- matrix(NaN, n, m)
  ySimHat1 <- matrix(NaN, n, m)
  yHatG <- matrix(NaN, n, m)
  ySimHatG <- matrix(NaN, n, m)
  eSimHat3 <- matrix(NaN, n, m)
  eHat3 <- matrix(NaN, n, m)
  eHatStar4 <- matrix(NaN, n, m)
  Q4global <- matrix(NaN, n, m)
  
  
  # --- Item 1 ---
  if (is.null(ySim)) 
    ySim <- matrix(rnorm(n * m), n, m)
  
  
  for (clust in unique(clusters)) {
    
    k <- clusters == clust
    
    yLocal <- y[k, , drop = FALSE]
    ySimLocal <- ySim[k, , drop = FALSE]
    
    # Local XA taking into account possible modification e) 
    x1 <- cbind(GroupModelMatrix(clusterPieces[k], xClusterPieces[k, , drop = FALSE]), xLocal[k, , drop = FALSE])
    
    
    # Use QR for calculation of regression fits
    xQ1 <- GenQR(x1, findR = FALSE)
    
    # --- Item 2 ---
    yHat1local <- xQ1 %*% (t(xQ1) %*% yLocal)
    ySimHat1local <- xQ1 %*% (t(xQ1) %*% ySimLocal)
    yHat1[k, ] <- yHat1local
    ySimHat1[k, ] <- ySimHat1local
    
    # --- Item 3 --- (with Item 2)
    xG <- xGlobal[k, , drop = FALSE]
    
    xGqr <- AdjustQR(xG, xOrth = xQ1)
    
    xGlobalOrth[k, ] <- xGqr$Q %*% xGqr$R
    
    # Use QR for calculation of regression fits
    xQG <- xGqr$Q
    
    # --- Item 4 ---
    yHatGlobal <- xQG %*% (t(xQG) %*% yLocal)
    ySimHatGlobal <- xQG %*% (t(xQG) %*% ySimLocal)
    ySimHatG[k, ] <- ySimHatGlobal
    yHatG[k, ] <- yHatGlobal
    
    # --- Item 5 ---
    eHat4local <- yLocal - yHat1local - yHatGlobal
    eSimHat4local <- ySimLocal - ySimHat1local - ySimHatGlobal
    
    # --- Item 6 ---
    R4 <- GenQR(eHat4local, makeunique = makeunique)$R
    Q4 <- GenQR(eSimHat4local, findR = FALSE, makeunique = makeunique)
    if (NCOL(Q4) < NROW(R4)) {
      stop("Something went wrong. Not enough dimensions in simulated data.")
    }
    eHatStar4[k, ] <- Q4[, seq_len(NROW(R4)), drop = FALSE] %*% R4
    # Q4 needed later if modification a) b) or c)
    Q4global[k, seq_len(NCOL(Q4))] <- Q4
  }
  # Use QR for calculation of regression fits
  xQG <- GenQR(xGlobalOrth, findR = FALSE)
  
  # --- Item 7 ---
  yHat2 <- xQG %*% (t(xQG) %*% y)
  ySimHat2 <- xQG %*% (t(xQG) %*% ySim)
  
  # --- Item 8 ---
  eHat3 <- yHatG - yHat2
  eSimHat3 <- ySimHatG - ySimHat2
  
  # --- Item 9 ---
  R3 <- GenQR(eHat3, makeunique = makeunique)$R
  Q3 <- GenQR(eSimHat3, findR = FALSE, makeunique = makeunique)
  if (NCOL(Q3) < NROW(R3)) {
    stop("Something went wrong. Not enough dimensions in simulated data.")
  }
  eHatStar3 <- Q3[, seq_len(NROW(R3)), drop = FALSE] %*% R3
  
  if (alternative == "a" | alternative == "b" | alternative == "c" | !is.null(alpha)) {
    if (!is.null(alpha)) 
      alternative <- "0"  # To avoid alternative=='c'
    if (alternative == "a") 
      alpha <- 0
    if (alternative == "b") 
      alpha <- 1
    if (alternative == "c") {
      alpha <- 1
      newAlpha <- 1  # new alpha will be computed
    }
    
    eHat <- y - yHat1 - yHat2
    for (clust in unique(clusters)) {
      k <- clusters == clust
      if (alternative == "a") 
        R4 <- GenQR(eHat[k, , drop = FALSE], makeunique = makeunique)$R else {
          if (alternative == "c") {
            cAlpha <- FindAlpha(eHat[k, , drop = FALSE], eHatStar3[k, , drop = FALSE])
            newAlpha <- min(cAlpha/(1 + epsAlpha), newAlpha)
            cat("\n ", sprintf("%15s: clusterAlpha = %7.3f, Alpha = %4.3f", clust, cAlpha, newAlpha))
          } else {
            R4 <- CalculateC(a = eHat[k, , drop = FALSE], b = eHatStar3[k, , drop = FALSE], epsAlpha = epsAlpha, alpha = alpha)
          }
          # Run function two times when c and alpha is NULL. First loop used to calculate alpha.
        }
      if (alternative != "c") {
        Q4 <- Q4global[k, seq_len(NROW(R4)), drop = FALSE]
        if (anyNA(Q4)) 
          stop("Not enough dimensions in simulated data.")
        eHatStar4[k, ] <- Q4 %*% R4
      }
    }
    if (alternative == "c") 
      return(RegSDChybrid(y = y, xGlobal = xGlobal, xLocal = xLocal, clusters = clusters, 
                          clusterPieces = clusterPieces, xClusterPieces = xClusterPieces, # groupedClusters already taken care of
                          alpha = newAlpha, ySim = ySim, returnParts = returnParts))
    eHatStar3 <- alpha * eHatStar3
  }
  if (returnParts) 
    return(list(y1 = yHat1, y2 = yHat2, e3s = eHatStar3, e4s = eHatStar4, e3 = eHat3, e4 = y - yHat1 - yHat2 - eHat3))
  return(yHat1 + yHat2 + eHatStar3 + eHatStar4)
}



