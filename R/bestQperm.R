#' @title Column permutation of the estimated Q-matrix
#'
#' @description
#' Function \code{bestQperm} is used to rearrange the columns of the estimated Q so that
#' the order of the columns best matches that of the true Q-matrix.
#'
#' @param estQ The estimated Q-matrix.
#' @param trueQ The true Q-matrix.
#'
#' @return The function returns a Q-matrix in which the order of the columns best matches
#' that of the true Q-matrix.
#'
#' @export
#'
#'
bestQperm = function(estQ,trueQ){
  K=ncol(estQ)
  combrow = gtools::permutations(K,K)
  result = list()
  empt1 = vector()
  for(r in 1:nrow(combrow)){
    temp1 = as.numeric(combrow[r,])
    temp2 = estQ[,temp1]
    temp3 = RR(trueQ,temp2)
    empt1 = rbind(empt1,temp3)
  }
  temp4 = which.max(empt1)
  temp5 = as.numeric(combrow[temp4,])
  resultQ = estQ[,temp5]
  #result[1] = list(resultQ)
  #result[2] = list(temp4)
  return(resultQ)
}
