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
#' \strong{sec7zAllSupp:} As Zall with suppressed values set to NA
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
#' RegSDCdata("sec7zAllSupp")
RegSDCdata <- function(dataset) {
  if (dataset %in% c("sec7data", "sec7y", "sec7x", "sec7xAll", "sec7z", "sec7zAll", "sec7zAllSupp") ) {
    y <- matrix(c(3, 1, 12, 18, 11, 9, 22, 19, 32, 13, 2, 16, 30, 8, 2, 3), ncol=1)
    cols <- paste("col", col(matrix(1, 4, 4)), sep = "")
    rows <- paste("row", row(matrix(1, 4, 4)), sep = "")
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
    if(dataset == "sec7zAll" | dataset == "sec7zAllSupp"){
      z <- as.matrix(t(m) %*% y)
      if(dataset == "sec7zAll")
        return(z)
      z[z[,1] %in% c(1, 2, 3, 9, 11, 13, 18),] <- NA
      return(z)
    }
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
