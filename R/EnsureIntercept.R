#' Ensure that input is matrix (by as.matrix) and check number of rows (and columns)
#'
#' @param x NULL or input to as.matrix
#' @param nRow Expected number of rows
#' @param nCol Expected number of columns
#'
#' @return Input as a matrix 
#' @export
#' @author Øyvind Langsrud
#' @keywords internal
#'
#' @examples
#' x <- matrix(c(5, 8, 4, 2, 7, 6), 3, 2)
#' EnsureMatrix(x)
#' EnsureMatrix(x, 3)
#' EnsureMatrix(1:4)
#' EnsureMatrix(1:4, 4)
#' EnsureMatrix(NULL, 4)
#' try(EnsureMatrix(x, 4))
#' try(EnsureMatrix(1:3, 4))
#' EnsureMatrix(x, 3, 2)
#' try(EnsureMatrix(x, 3, 3))
#' try(EnsureMatrix(NULL, 3, 3))
EnsureMatrix <- function(x, nRow = NULL, nCol = NULL) {
  if (is.null(x)) {
    if(!is.null(nCol))
      if(nCol)
        stop("nCol not allowed when x is NULL")
    return(matrix(0, nRow, 0))
  }
  x <- as.matrix(x)
  if (!is.null(nRow)) 
    if (nrow(x) != nRow) 
      stop(paste("nrow is", nrow(x), "when", nRow, "expected"))
  if (!is.null(nCol)) 
    if (ncol(x) != nCol) 
      stop(paste("ncol is", ncol(x), "when", nCol, "expected"))
  x
}


#' Ensure constant term in matrix
#' 
#' A column of ones may be added
#'
#' @param x Input matrix 
#'
#' @return The input matrix possibly with a column of ones added
#' @export
#' @author Øyvind Langsrud
#' @keywords internal
#'
#' @examples
#' x <- matrix(c(5, 8, 4, 2, 7, 6), 3, 2)
#' EnsureIntercept(x)
#' EnsureIntercept(cbind(x, 2))
#' EnsureIntercept(cbind(x, 0))
#' EnsureIntercept(matrix(0, 4, 0))
EnsureIntercept <- function(x) {
  colMax <- apply(x, 2, max)
  colMin <- apply(x, 2, min)
  if (!any((colMax - colMin == 0) & (colMax != 0))) 
    return(cbind(1, x))
  x
}



