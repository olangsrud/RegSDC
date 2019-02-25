#' Function that returns a dataset
#' 
#' 
#' @param dataset Name of data set within the RegSDC package
#'
#' @return data frame
#' 
#' @details 
#' \strong{sec7data:} Data in section 7 of the paper as a data frame
#' 
#' \strong{sec7y:} Y in section 7 of the paper as a matrix
#' 
#' \strong{sec7x:} X in section 7 of the paper as a matrix
#' 
#' \strong{sec7z:} Z in section 7 of the paper as a matrix
#' 
#' \strong{sec7xAll:} Xall in section 7 of the paper as a matrix
#' 
#' \strong{sec7zAll:} Zall in section 7 of the paper as a matrix
#' 
#' @importFrom SSBtools Hierarchies2ModelMatrix
#' @import Matrix
#' @export
#' @author Øyvind Langsrud
#'
#' @examples
#' RegSDCdata("sec7data")
#' RegSDCdata("sec7y")
#' RegSDCdata("sec7x")
#' RegSDCdata("sec7z")
#' RegSDCdata("sec7xAll")
#' RegSDCdata("sec7zAll")
RegSDCdata <- function(dataset) {
  if (dataset == "sec7data" | dataset == "sec7y" | dataset == "sec7x" | dataset == "sec7xAll"
      | dataset == "sec7z" | dataset == "sec7zAll" ) {
    y <- matrix(c(3, 1, 12, 18, 11, 9, 22, 19, 32, 13, 2, 16, 30, 8, 2, 3), ncol=1)
    cols <- paste("col", col(matrix(1:12, 4, 4)), sep = "")
    rows <- paste("row", row(matrix(1:12, 4, 4)), sep = "")
    rownames(y) <- paste(rows,cols,sep="_")
    colnames(y) <- "freq"
    if(dataset == "sec7y")
      return(y)
    inner <- data.frame(rows, cols, y)
    h2m <- Hierarchies2ModelMatrix(inner[, 1:2], list(rows = "Total", cols = "Total"), crossTable = TRUE)
    m <- h2m$modelMatrix
    rownames(m) <- rownames(y)
    if(dataset == "sec7xAll")
      return(as.matrix(m))
    if(dataset == "sec7zAll")
      return(as.matrix(t(m) %*% y))
    z <- as.vector(as.matrix(t(m) %*% y))
    suppressed <- z
    suppressed[z %in% c(1, 2, 3, 9, 11, 13, 18)] <- NA
    if(dataset == "sec7x")
      return(as.matrix(m)[,!is.na(suppressed)])
    if(dataset == "sec7z")
      return(as.matrix(t(as.matrix(m)[,!is.na(suppressed)]) %*% y))
    y2 <- z
    y2[(h2m$crossTable[, 1] %in% "Total") | (h2m$crossTable[, 2] %in% "Total")] <- NA
    return(cbind(h2m$crossTable, y = y2, z = z, suppressed = suppressed))
  }
  stop(paste("No data with dataset =", dataset))
}
