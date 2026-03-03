#' @title Detect the implausible and improper distractors in a Q-Matrix for multiple-choice items
#' @description
#' Function \code{distractor.check} is used to assess whether the distractors of a given Q-matrix for multiple-choice items are plausible and/or proper.
#'
#' @param Q The given Q-matrix for multiple-choice items. It has to be organized in the following manner. The Q-matrix should contain \eqn{K+2} columns
#' and \eqn{J \times (H_j)} rows, where \eqn{H_j} is the number of coded options for item \eqn{j}. Among the \eqn{H_j} rows for item \eqn{j},
#' the first row is the key, followed by the coded distractors.
#' The first column of the Q-matrix lists the ID of the items and second column indicates the corresponding options of the coded options.
#' If the Q-matrix is not organized in such a way, an argument \code{key} (see below) that indicates the options that the keys
#' are located needs to be given.
#' @param key A vector that indicates the options where the keys are located.
#'
#' @return A list with class "distractor.check" containing:
#' \describe{
#'   \item{not.plausible}{A matrix indicating items and options that are not plausible, or NULL if all items are plausible.}
#'   \item{not.proper}{A vector of item IDs that are not proper, or NULL if all items are proper.}
#'   \item{all.plausible}{Logical; TRUE if all items are plausible.}
#'   \item{all.proper}{Logical; TRUE if all items are proper.}
#' }
#'
#' @references
#' Chiu, C.-Y., Köhn, H. F. & Wang, Y. (Online first). Plausible and proper
#' multiple-choice items for diagnostic classification. \emph{Psychometrika}.
#' 
#' @examples
#' \dontrun{
#' library(NPCDTools)
#' Q1 <- Q_Ozaki
#' distractor.check(Q1)
#' 
#' Q2 <- GDINA::sim10MCDINA2$simQ
#' key <- c(1, 2, 4, 1, 1, 3, 2, 4, 1, 4)
#' distractor.check(Q2, key)
#' }
#' 
#' @export
#'
#'
#'
#'
distractor.check = function(Q = Q, key = NULL){
  
  # Input validation
  if (missing(Q)) {
    stop("Argument 'Q' is missing. Please provide a Q-matrix.")
  }
  
  if (!is.matrix(Q) && !is.data.frame(Q)) {
    stop("Q must be a matrix or data frame.")
  }
  
  if (is.data.frame(Q)) {
    Q <- as.matrix(Q)
  }
  
  if (ncol(Q) < 3) {
    stop("Q must have at least 3 columns (Item ID, Option, and at least one attribute).")
  }
  
  if (nrow(Q) < 1) {
    stop("Q must have at least 1 row.")
  }
  
  if (!is.null(key)) {
    if (!is.numeric(key) && !is.integer(key)) {
      stop("'key' must be a numeric vector.")
    }
    if (length(key) != length(unique(Q[,1]))) {
      stop("Length of 'key' must equal the number of unique items in Q.")
    }
  }
  
  # Check if first column is monotonically non-decreasing (items should be grouped together)
  if (any(diff(Q[,1]) < 0)) {
    stop("Please make sure the input Q-matrix is MC Q-matrix.")
  }
  
  # Original algorithm
  if (is.null(key)==F){
    J = length(key)
    Q.temp = NULL
    for (j in 1:J){
      q = Q[which(Q[,1]==j),]
      if (is.matrix(q)==T){
        qq = rbind(q[which(q[,2]==key[j]),], q[-which(q[,2]==key[j]),])
        Q.temp = rbind(Q.temp, qq)
      }else {
        Q.temp = rbind(Q.temp, q)
      }
    }
    Q = Q.temp
  }
  Q=as.matrix(Q)
  J=length(unique(Q[,1]))
  K=ncol(Q)-2
  
  flag.plausible = flag.proper = NULL
  for (j in 1:J){
    subQ=Q[which(Q[,1]==j),]
    if (is.matrix(subQ)==T){
      ## Check plausible. Note: If (B-A)A=0, A is nested in B
      dif = matrix(rep(subQ[1,3:ncol(Q)],nrow(subQ)-1), nrow(subQ)-1, K, byrow=T) - subQ[2:nrow(subQ),3:ncol(Q)]
      if (length(which(subQ[1,3:ncol(Q)]%*%(t(abs(dif)))==0))!=0) flag.plausible = rbind(flag.plausible, cbind(j,which(subQ[1,3:ncol(Q)]%*%(t(dif))==0)))
      ## Check proper
      if (nrow(subQ)>2){
        # index of distractors
        index = as.matrix(utils::combn(2:nrow(subQ), 2))
        for (k in 1:ncol(index)){
          union=(subQ[index[1,k],3:ncol(Q)] | subQ[index[2,k],3:ncol(Q)])*1
          # Check whether the union is in subQ or the key is nested in the union
          condition = which(rowSums(abs(subQ[,3:ncol(Q)] - matrix(rep(union, nrow(subQ)), nrow(subQ), ncol(Q)-2, byrow=T)))==0)
          if (length(condition)==0 & (union-subQ[1, 3:ncol(subQ)])%*%subQ[1, 3:ncol(subQ)]!=0){
            flag.proper = c(flag.proper, j)
            break
          }else{
            k = k+1
          }
        }
      }
    }
  }
  
  # Create result object
  result <- list(
    not.plausible = flag.plausible,
    not.proper = flag.proper,
    all.plausible = is.null(flag.plausible),
    all.proper = is.null(flag.proper)
  )
  
  class(result) <- "distractor.check"
  return(result)
}

#' @title Print the Summary of Distractor Check Results
#'
#' @description
#' Prints a summary of whether the distractors of a given Q-matrix for multiple-choice items are plausible and/or proper.
#'
#' @param x An object of class "distractor.check" returned by \code{\link{distractor.check}}.
#' @param ... Additional arguments (currently not used).
#'
#'
#' @method print distractor.check
#' @rdname print.distractor.check
#' @export
print.distractor.check <- function(x, ...) {
  cat("\n")
  
  if (x$all.plausible & x$all.proper) {
    cat("All the items are plausible and proper.\n")
  } else if (x$all.plausible & !x$all.proper) {
    cat("All the items are plausible.\n")
    cat("However, items", paste(unique(x$not.proper), collapse = ", "), "are not proper.\n")
  } else if (!x$all.plausible & x$all.proper) {
    cat("All the items are proper.\n")
    cat("However, items", paste(unique(x$not.plausible[,1]), collapse = ", "), "are not plausible.\n")
  } else {
    cat("Items", paste(unique(x$not.plausible[,1]), collapse = ", "), "are not plausible.\n")
    cat("Items", paste(unique(x$not.proper), collapse = ", "), "are not proper.\n")
  }
  
  cat("\n")
  invisible(x)
}

