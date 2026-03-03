#' @title Entry-wise and vector-wise agreement rate between two Q-matrices
#'
#' @description
#' \code{RR} is used to compute the agreement rate between two Q-matrices with 
#' identical dimensions.
#'
#' @param Q1 The first Q-matrix.
#' @param Q2 The second Q-matrix that has the same dimensionality as \code{Q1}.
#'
#' @return The function returns
#' \describe{
#' \item{entry.wise}{The entry-wise agreement rate}
#' \item{item.wise}{The item-wise agreement rate}
#' }
#'
#' @export
#' @seealso See the examples for the \code{\link{QR}} and \code{\link{TSQE}} functions.
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

