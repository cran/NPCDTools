#' @title Column permutation of a Q-matrix with respect to a benchmark Q-matrix.
#'
#' @description
#' Function \code{bestQperm} is used to permute the columns of a Q-matrix so that
#' the order of the columns best matches that of the benchmark Q-matrix. This function
#' is useful in a Q-matrix estimation process.
#'
#' @param Q The targeted Q-matrix.
#' @param bench.Q The benchmark Q-matrix.
#'
#' @return The function returns a Q-matrix in which the order of the columns best matches
#' that of the benchmark Q-matrix.
#' @examples
#' # See examples used for TSQE.
#' 
#' @export
#'
#'
bestQperm = function(Q,bench.Q){
  K=ncol(Q)
  combrow = gtools::permutations(K,K)
  result = list()
  empt1 = vector()
  for(r in 1:nrow(combrow)){
    temp1 = as.numeric(combrow[r,])
    temp2 = Q[,temp1]
    temp3 = RR(bench.Q,temp2)
    empt1 = rbind(empt1,temp3)
  }
  temp4 = which.max(empt1)
  temp5 = as.numeric(combrow[temp4,])
  resultQ = Q[,temp5]
  #result[1] = list(resultQ)
  #result[2] = list(temp4)
  return(resultQ)
}
