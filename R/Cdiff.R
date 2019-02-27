#' Matrix difference (a-b) including checking for equal columns
#' 
#' Each column is checked by \code{\link{all.equal}}
#'
#' @param a numerical matrix
#' @param b numerical matrix
#' @param tolerance parameter to \code{\link{all.equal}}
#'
#' @return (a-b) where equal columns are set to zero
#' @keywords internal
#' @export
#'
#' @examples
#' a <- matrix(rnorm(6), 3, 2)
#' b <- matrix(rnorm(6), 3, 2)
#' a - b
#' Cdiff(a, b)
#' b[, 1] <- a[, 1] + 1e-10 * b[, 1]
#' a - b
#' Cdiff(a, b)
#' a[, 2] <- b[, 2]
#' a - b
#' Cdiff(a, b)
Cdiff = function (a,b, tolerance = sqrt(.Machine$double.eps)){
  d = a-b
  n = ncol(d)
  equalCol = rep(FALSE, n)
  for(i in seq_len(n)){
    equalCol[i] = isTRUE(all.equal(a[,i],b[,i], tolerance =tolerance))
  }
  d[, equalCol] = 0
  d
}





