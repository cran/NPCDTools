#' @title Entry-wise and vector-wise recovery rates
#'
#' @description
#' Function \code{RR} is used to compute the recovery rates for an estimate Q-matrix.
#' In general, it can be used to compute the agreement rate between two matrices with identical dimensionalities.
#'
#' @param Q1 The first Q-matrix.
#' @param Q2 The second Q-matrix that has the same dimensionality as Q1.
#'
#' @return The function returns
#' \describe{
#' \item{entry.wise}{The entry-wise recovery rate}
#' \item{item.wise}{The item-wise recovery rate}
#' }
#'
#' @export
#'
#'
RR = function(Q1,Q2){
  if (identical(dim(Q1),dim(Q2)) != TRUE){
    stop("The dimension for these two matrix should be same.")
  }
  J = nrow(Q1); K = ncol(Q1)
  tmp1 = sum(Q1==Q2) / (J*K) ## by element
  tmp2 = sum(apply((Q1==Q2)+0,1,prod))/J ## by row
  info = list(entry.wise = tmp1, item.wise = tmp2)
  return(info)
}

