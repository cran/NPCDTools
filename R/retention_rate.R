#' @title Retention rate of a Q-matrix refinement method
#'
#' @description This function computes the proportion of correctly specified q-entries
#' in a provisional Q-matrix that remain correctly specified after a Q-matrix refinement
#' procedure is applied. This function is used only when the true Q-matrix is known.
#'
#' @param ref.Q the \eqn{J \times K} binary Q-matrix obtained from applying a refinement procedure.
#' @param mis.Q A \eqn{J \times K} binary provisional Q-matrix.
#' @param true.Q The \eqn{J \times K} binary true Q-matrix.
#'
#' @return The function returns a value between 0 and 1 indicating the proportion of
#' correctly specified q-entries in \code{mis.Q} that remain correctly specified in \code{ref.Q}
#' after a Q-matrix refinement procedure is applied to \code{mis.Q}.
#'
#' @export
#'
#'
retention.rate = function(ref.Q = ref.Q, mis.Q = mis.Q, true.Q = true.Q)
{
  mis=mis.Q-true.Q
  rec=ref.Q-true.Q
  loc.cor=which(mis==0,arr.ind = TRUE)
  count=0
  for (m in 1:nrow(loc.cor)) if (rec[loc.cor[m,1],loc.cor[m,2]]==0) count=count+1
  ret.rate = count/nrow(loc.cor)
  return(ret.rate)
}
