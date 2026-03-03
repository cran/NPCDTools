#' @title Generation of a Q-matrix with improper MC items that are not proper based on a 'plausible and proper' Q
#'
#' @description This function turns a proper and plausible Q-matrix into an improper Q-matrix with the user-specified
#' number of improper distractors.
#'
#' @param Q A proper and plausible Q-matrix for MC items.
#' @param n.improper The number of improper distractors.
#'
#' @return The function returns
#' \describe{
#' \item{Q.improper}{The improper Q-matrix generated from Q}
#' \item{item.improper}{The ID of the items that are improper.}
#' }
#'
#' @references
#' Chiu, C.-Y., Köhn, H. F. & Wang, Y. (Online first). Plausible and proper multiple-choice items 
#' for diagnostic classification. \emph{Psychometrika}.
#' \doi{10.1017/psy.2025.10074}
#' @export
#'
#'
Q.improper=function(Q, n.improper){
  Q.improper = Q
  K = ncol(Q.improper)-2
  J = max(Q.improper[,1])
  M = 2^K
  
  pattern <- diag(K)
  for (l in 2:K) {
    pattern <- rbind(pattern, t(apply(utils::combn(K, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  
  can.improper=NULL
  for (j in 1:J){
    temp = which(Q.improper[,1]==j)
    n.temp = length(temp)
    key = Q.improper[temp[1],-(1:2)]
    if (length(temp) > 2 & sum(key)!=K) can.improper = c(can.improper, j)
  }
  
  if(is.null(can.improper)==T){
    return(cat("The phenomenon of having improper items will not occur.","\n"))
  }
  
  ## If length(can.improper) < n.improper, length(can.improper) is used to generate Q.improper.
  if (length(can.improper) < n.improper){
    n.improper = length(can.improper)
    cat("The maximal number of items that improperness could happen (",length(can.improper), ") is less than the specified n.improper.","\n",
        "The maximal number is used to generate the improper items.","\n", sep="")
  }
  
  item.improper=NULL
  for (j in 1:n.improper){
    if (length(can.improper)==1) jj = can.improper else{
      jj=sample(can.improper,1)
      can.improper=can.improper[-which(can.improper==jj)]
    }
    
    index = which(Q.improper[,1]==jj)
    subQ=Q.improper[index,]
    item.improper=c(item.improper, jj)
    used.option = subQ[,2]
    
    ## Those that are nested within the key and those that the key is nested with are eliminated
    ## If (A+B>0)*1 = A, B is a subset of A
    ## If (B-A)A = 0, A is a subset of B
    key=subQ[1,-(1:2)]
    dif = matrix(rep(key,M),M, ncol(subQ)-2,byrow=T) - pattern
    nested = which(dif%*%key==0 & rowSums(abs(dif))!=0) # the key is nested within
    #nest = which(rowSums(((matrix(rep(key,M),M, ncol(subQ)-2,byrow=T) + pattern)>0)*1 - key)==0) # nest within the key
    rest = (1:M)[-nested]
    
    for (i in 1:length(index)){
      dif2 = matrix(rep(subQ[i, -(1:2)],M),M, ncol(subQ)-2,byrow=T) - pattern
      rest = rest[rest%in%which(rowSums(abs(dif2))==0)==F]
    }
    
    ##======================
    # For a special case like this, the selected.option can only be the distractor 11011
    # Option Att1 Att2 Att3 Att4 Att5
    #      1    0    0    1    0    0
    #      3    1    1    0    1    1
    #      2    0    1    0    0    0
    #====================================================================================
    selected.option = sample(index[-1],1)
    for (i in 2:length(index)){
      if (sum(subQ[1, -(1:2)])==1 & sum(abs(subQ[i, -(1:2)] + subQ[1, -(1:2)] - rep(1,K)))==0) selected.option = index[i]
    }
    
    n=sample(rest, 1)
    Q.improper[selected.option, -(1:2)]=pattern[n,]
    new.subQ=Q.improper[index,]
    for (k in (2:length(index))[-which(index[-1]==selected.option)]){
      ## (pattern[n,]+subQ[k,-(1:2)]>0)*1 is the union
      ind2 = rowSums(abs(matrix(rep((pattern[n,]+subQ[k,-(1:2)]>0)*1, nrow(subQ)), nrow(subQ), ncol(subQ)-2,byrow=T)-new.subQ[,-(1:2)]))
      while (0%in%ind2==T | ((pattern[n,]+subQ[k,-(1:2)]>0)*1 - key)%*%key==0){
        n=sample(rest[-which(rest==n)], 1)
        Q.improper[selected.option, -(1:2)]=pattern[n,]
        new.subQ=Q.improper[index,]
        ind2 = rowSums(abs(matrix(rep((pattern[n,]+subQ[k,-(1:2)]>0)*1, nrow(subQ)), nrow(subQ), ncol(subQ)-2,byrow=T)-new.subQ[,-(1:2)]))
      }
    }
    Q.improper=as.matrix(Q.improper)
  }
  return(list(Q.improper=Q.improper, item.improper = item.improper))
}
