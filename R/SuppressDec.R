

ColRepMatrix <- function(x, nRep) {
  if (nRep == 1) 
    return(x)
  if (nRep == 0) 
    return(x[, integer(0), drop = FALSE])
  cn <- colnames(x)
  rn <- rownames(x)
  x <- matrix(x, nrow(x), ncol(x) * nRep)
  if (is.null(cn)) 
    colnames(x) <- cn
  if (is.null(rn)) 
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