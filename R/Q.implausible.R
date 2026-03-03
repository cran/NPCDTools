#' @title Generation of a Q-matrix that contains implausible MC items
#'
#' @description This function turns a proper and plausible Q-matrix into an implausible Q-matrix with the 
#' user-specified number of implausible distractors.
#'
#' @param Q A proper and plausible Q-matrix for MC items.
#' @param n.implausible The number of items that have implausible distractors.
#'
#' @return The function returns
#' \describe{
#' \item{Q.implausible}{The generated Q-matrix}
#' \item{item.implausible}{The ID of the items that have implausible distractors.}
#' }
#'
#' @references
#' Chiu, C.-Y., Köhn, H. F. & Wang, Y. (Online first). Plausible and proper multiple-choice items 
#' for diagnostic classification. \emph{Psychometrika}.
#' \doi{10.1017/psy.2025.10074}
#' 
#' @export
#'
#'
Q.implausible = function(Q, n.implausible){
  Q.implausible = Q
  K=ncol(Q.implausible)-2
  J=max(Q.implausible[,1])
  pattern <- diag(K)
  for (l in 2:K) {
    pattern <- rbind(pattern, t(apply(utils::combn(K, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  M <- nrow(pattern)
  
  can.implausible = NULL
  for (j in 1:J){
    temp = which(Q.implausible[,1]==j)
    n.temp = length(temp)
    key = Q.implausible[temp[1],-(1:2)]
    if (length(temp) >1 & sum(key)!=K) can.implausible = c(can.implausible, j)
  }
  
  if(is.null(can.implausible)==T){
    return(cat("The phenomenon of having implausible items will not occur.","\n"))
  }
  
  item.plausible=NULL
  for (j in 1:n.implausible){
    if (length(can.implausible)==1) jj=can.implausible else{
      jj=sample(can.implausible, 1)
      can.implausible=can.implausible[-which(can.implausible==jj)]
    }
    
    item.plausible=c(item.plausible, jj)
    index = which(Q.implausible[,1]==jj)
    subQ=Q.implausible[index,]
    
    used.option = subQ[,2]
    
    key=subQ[1,-(1:2)]
    dif = matrix(rep(key,M),M, ncol(subQ)-2,byrow=T) - pattern
    nest = which(dif%*%key==0 & rowSums(abs(dif))!=0)
    if (length(nest)==1) new.distractor = pattern[nest,] else new.distractor = pattern[sample(nest, 1),]
    
    if (length(index)==2) replaced.option = index[2] else replaced.option=sample(index[-1], 1)
    
    Q.implausible[replaced.option, -(1:2)] = new.distractor
    Q.implausible=as.matrix(Q.implausible)
  }
  output = list(Q.implausible=Q.implausible, item.plausible = item.plausible)
  return(output)
}
