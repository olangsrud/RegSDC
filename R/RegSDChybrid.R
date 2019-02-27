
  
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
#' @details Input matrices are subjected to \code{\link{EnsureMatrix}}.
#' Necessary constant terms (intercept) are automatically included. 
#' That is, a column of ones is not needed in the input matrices.
#'
#' 
#' @return Generated version of y
#' @importFrom stats model.matrix
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' #################################################
#' # Generate example data for introductory examples
#' ################################################# 
#' y <- matrix(rnorm(30) + 1:30, 10, 3)
#' x <- matrix(1:10, 10, 1)  # x <- 1:10 is equivalent
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
#' ###################################################
#' #  Generate example data for more advanced examples
#' ###################################################
#' x <- matrix((1:90) * (1 + runif(90)), 30, 3)
#' x1 <- x[, 1]
#' x2 <- x[, 2]
#' x3 <- x[, 3]
#' y <- matrix(rnorm(90), 30, 3) + x
#' clust <- paste("c", rep(1:3, each = 10), sep = "")
#' 
#' ######## Run main algorithm
#' z0 <- RegSDChybrid(y, clusters = clust, xLocal = x3, xGlobal = cbind(x1, x2))
#' 
#' # Corresponding models by lm
#' lmy <- lm(y ~ clust + x1 + x2 + x3:clust)
#' lm0 <- lm(z0 ~ clust + x1 + x2 + x3:clust)
#' 
#' # Preserved regression coef (x3 within clusters)
#' coef(lmy) - coef(lm0)
#' 
#' # Preservation of x3 coef locally can also be seen by local regression
#' coef(lm(y ~ x3, subset = clust == "c2")) - coef(lm(z0 ~ x3, subset = clust == "c2"))
#' 
#' # Covariance matrix preserved
#' cov(resid(lmy)) - cov(resid(lm0))
#' 
#' # But not preserved within clusters
#' cov(resid(lmy)[clust == "c2", ]) - cov(resid(lm0)[clust == "c2", ])
#' 
#' ######## Modification (a)
#' za <- RegSDChybrid(y, clusters = clust, xLocal = x3, xGlobal = cbind(x1, x2), alternative = "a")
#' lma <- lm(za ~ clust + x1 + x2 + x3:clust)
#' 
#' # Now covariance matrices preserved within clusters
#' cov(resid(lmy)[clust == "c2", ]) - cov(resid(lma)[clust == "c2", ])
#' 
#' # If we estimate coef for x1 and x2 within clusters, 
#' # they become identical and identical to global estimates
#' coef(lma)
#' coef(lm(za ~ clust + x1:clust + x2:clust + x3:clust))
#' 
#' ######## Modification (c) with automatic calculation of alpha 
#' # The result depends on the randomly generated data
#' # When the result is that alpha=1, modification (b) is equivalent
#' zc <- RegSDChybrid(y, clusters = clust, xLocal = x3, xGlobal = cbind(x1, x2), alternative = "c")
#' lmc <- lm(zc ~ clust + x1 + x2 + x3:clust)
#' 
#' # Preserved regression coef as above
#' coef(lmy) - coef(lmc)
#' 
#' # Again covariance matrices preserved within clusters
#' cov(resid(lmy)[clust == "c2", ]) - cov(resid(lmc)[clust == "c2", ])
#' 
#' # If we estimate coef for x1 and x2 within clusters, 
#' # results are different from modification (a) above
#' coef(lmc)
#' coef(lm(zc ~ clust + x1:clust + x2:clust + x3:clust))
#' 
#' 
#' ####################################################
#' # Make groups of clusters (d) and cluster pieces (e)
#' ####################################################
#' clustGr <- paste("gr", ceiling(rep(1:3, each = 10)/2 + 0.1), sep = "")
#' clustP <- c("a", "a", rep("b", 28))
#' 
#' ######## Modifications (c), (d) and (e)
#' zGrP <- RegSDChybrid(y, clusters = clust, clusterPieces = clustP, groupedClusters = clustGr,
#'                      xLocal = x3, xGroupedClusters = x2, xGlobal = x1, alternative = "c")
#' 
#' # Corresponding models by lm
#' lmGrP <- lm(zGrP ~ clust:clustP + x1 + x2:clustGr + x3:clust - 1)
#' lmY <- lm(y ~ clust:clustP + x1 + x2:clustGr + x3:clust - 1)
#' 
#' # Preserved regression coef
#' coef(lmY) - coef(lmGrP)
#' 
#' # Identical means within cluster pieces
#' aggregate(y, list(clust = clust, clustP = clustP), mean)
#' aggregate(zGrP, list(clust = clust, clustP = clustP), mean)
#' 
#' # Covariance matrix preserved
#' cov(resid(lmY)) - cov(resid(lmGrP))
#' 
#' # Covariance matrices preserved within clusters
#' cov(resid(lmY)[clust == "c2", ]) - cov(resid(lmGrP)[clust == "c2", ])
#' 
#' # Covariance matrices not preserved within cluster pieces
#' cov(resid(lmY)[clustP == "a", ]) - cov(resid(lmGrP)[clustP == "a", ])
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
    #eHat4local <- yLocal - yHat1local - yHatGlobal
    #eSimHat4local <- ySimLocal - ySimHat1local - ySimHatGlobal
    eHat4local <- Cdiff(yLocal , yHat1local + yHatGlobal)
    eSimHat4local <- Cdiff(ySimLocal , ySimHat1local + ySimHatGlobal)
    
    
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
  #eHat3 <- yHatG - yHat2
  #eSimHat3 <- ySimHatG - ySimHat2
  eHat3 <- Cdiff(yHatG , yHat2)
  eSimHat3 <- Cdiff(ySimHatG , ySimHat2)
  
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
            if(newAlpha < epsAlpha)
              newAlpha <- 0
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
    if (alternative == "c"){ 
      cat("\n")
      return(RegSDChybrid(y = y, xGlobal = xGlobal, xLocal = xLocal, clusters = clusters, 
                          clusterPieces = clusterPieces, xClusterPieces = xClusterPieces, # groupedClusters already taken care of
                          alpha = newAlpha, ySim = ySim, returnParts = returnParts))
      }
    eHatStar3 <- alpha * eHatStar3
  }
  if (returnParts) 
    return(list(y1 = yHat1, y2 = yHat2, e3s = eHatStar3, e4s = eHatStar4, e3 = eHat3, e4 = y - yHat1 - yHat2 - eHat3))
  return(yHat1 + yHat2 + eHatStar3 + eHatStar4)
}



