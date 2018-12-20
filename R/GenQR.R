#' Generalized QR decomposition
#' 
#' Matrix X decomposed as Q and R (X=QR) where columns of Q are orthogonal.
#' Ordinary QR or SVD may be used. 
#' 
#' @encoding UTF8
#'
#' @param x Matrix to be decomposed
#' @param doSVD When TRUE SVD instead of QR
#' @param findR When FALSE only Q returned
#' @param makeunique When TRUE force uniqueness by positive diagonal elements (QR) or by column sums (SVD)
#' @param tol As input to qr or, in the case of svd(), similar as input to MASS::ginv().
#'  
#' @details
#' To handle dependency a usual decomposition of X is PX=QR where P is a permutation matrix. 
#' This function returns RP^T as R. When SVD, Q=U and R=SV^T.
#'
#' @return List with Q and R or just Q
#' @export
#'
#' @examples
#'    GenQR(matrix(rnorm(15),5,3))
#'    GenQR(matrix(rnorm(15),5,3)[,c(1,2,1,3)])
#'    GenQR(matrix(rnorm(15),5,3)[,c(1,2,1,3)],TRUE)
GenQR = function(x,doSVD=FALSE,findR=TRUE,makeunique=findR,tol = 1e-07){
  if(is.null(x)) return(NULL)
  if(doSVD){  # inspierd by ginv i MASS
    xSvd <- svd(x)
    Positive <- xSvd$d > max(tol * xSvd$d[1L], 0)
    Q = xSvd$u[, Positive, drop = FALSE]
    if(findR|makeunique){ 
      if(makeunique){
        v = t(xSvd$v[, Positive, drop = FALSE])
        #print(v)
        sgn =  sign(rowSums(v))
        Q   <- Q %*% Diag(sgn)
        if(!findR) return(Q)
        R = sgn*xSvd$d[Positive]*v
      } else
        R = xSvd$d[Positive] * t(xSvd$v[, Positive, drop = FALSE])
    }
    else return(Q)
    return(list(Q=Q,R=R))
  }
  qrX = qr(x,tol=tol)
  Q = qr.Q(qrX)[,seq_len(qrX$rank),drop=FALSE]
  if(findR|makeunique) {
    R = qr.R(qrX)
    if(makeunique){
      R = R[seq_len(qrX$rank), ,drop=FALSE]
      sgn <- sign(Diag(R))
      Q   <- Q %*% Diag(sgn)
      if(!findR) return(Q)
      R <- Diag(sgn) %*% R
      R = R[ ,sort.list(qrX$pivot),drop=FALSE]
    }
    else
      R = R[seq_len(qrX$rank),sort.list(qrX$pivot),drop=FALSE]
  }
  else 
    return(Q)
  return(list(Q=Q,R=R))
}

# See note in documentation od diag
Diag = function(x){
  if(is.matrix(x)) return(diag(x))
  diag(x, nrow = length(x))
}


GenQRQR = function(x,...){
  z=GenQR(x,...)
  z$Q%*%z$R
}

if(FALSE){
  x=matrix(rnorm(15),5,3)[,c(1,2,1,3)]
  for(a in c(FALSE,TRUE))
    for(b in c(FALSE,TRUE))
      print(max(abs(LangsrudWorks:::GenQRQR(x,a,TRUE,b) -x)))
  x=matrix(rnorm(15),5,3)[,c(1,2,1,3)]
  for(a in c(FALSE,TRUE))
    print(GenQR(x,a,FALSE,TRUE)/GenQR(x,a,FALSE,FALSE))
  }
