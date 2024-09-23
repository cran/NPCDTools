#' @title Correction rate of a Q-matrix refinement method
#'
#' @description This function computes the proportion of corrected q-entries that
#' were originally misspecified in the provisional Q-matrix. This function is used
#' only when the true Q-matrix is known.
#'
#' @param ref.Q the \eqn{J \times K} binary Q-matrix obtained from applying the refinement procedure.
#' @param mis.Q A \eqn{J \times K} binary provisional Q-matrix.
#' @param true.Q The \eqn{J \times K} binary true Q-matrix.
#'
#' @return The function returns a value between 0 and 1 indicating the proportion of corrected q-entries in \code{ref.Q}
#' that were originally missepcified in \code{mis.Q}.
#'
#' @export
#'
#'
correction.rate = function(ref.Q = ref.Q, mis.Q = mis.Q, true.Q = true.Q)
{
  if (sum(abs(mis.Q - true.Q))==0) cor.rate = "NA" else
  {
    mis=mis.Q-true.Q
    rec=ref.Q-true.Q
    loc.miss=which(mis!=0,arr.ind = TRUE)
    count=0
    for (m in 1:nrow(loc.miss)) if (rec[loc.miss[m,1],loc.miss[m,2]]==0) count=count+1
    cor.rate = count/nrow(loc.miss)
    return(cor.rate)
  }
}
