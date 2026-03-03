#' @title Q-matrix refinement method
#'
#' @description The \code{QR} function refines a provisional Q-matrix by minimizing the residual sum of squares (RSS)
#'     between the observed and ideal item responses across all possible q-vectors, given the estimates of
#'     examinees' attribute profiles. 
#'
#' @details
#' This function implements the Q-matrix refinement (QR) method developed by Chiu
#' (2013). The NPC method (Chiu & Douglas, 2013) is first used to classify examinees and the best q-vector
#' for an item is identified by minimizing its RSS. Specifically, the RSS of 
#' item \eqn{j}{} for examinee \eqn{i}{} is defined as
#' \deqn{RSS_{j} =\sum_{m=1}^{2^K} \sum_{i \in C_{m}} (Y_{ij} - \eta_{jm})^2,}{}
#' where \eqn{C_m}{} for \eqn{m = 1, \ldots, 2^K} is the \eqn{m}{}th proficiency class, and \eqn{N}{}
#' is the number of examinees. Chiu (2013) proved that the expected value of
#' \eqn{RSS_j}{} corresponding to the correct q-vector is the minimum among the
#' \eqn{2^K - 1}{} candidates.
#'
#'
#' @param Y A \eqn{N\times J} matrix of binary responses (1=correct, 0=incorrect). Rows
#'     represent persons and columns represent items.
#' @param Q A \eqn{J\times K} provisional Q-matrix to be refined. Rows represent items and columns
#'     represent attributes.
#' @param gate A string, "\code{AND}" or "\code{OR}". "\code{AND}" is specified when a conjunctive 
#' relation between the attributes an examinee possesses and the attributes required by an item is assumed. 
#' "\code{OR}" is specified when a disjunctive relation between the attributes an examinee possesses  
#' and the attributes required by an item is assumed.
#' @param max.ite The number of iterations to run until the RSS's of all items
#'     are stationary.
#'
#' @return A list containing:
#' \item{initial.class}{Initial classifications of examinees}
#' \item{terminal.class}{Terminal classification of examinees}
#' \item{modified.Q}{The modified Q-matrix}
#' \item{modified.entries}{The modified q-entries}
#'
#' @seealso \code{\link{NPC}}
#' @references
#' Chiu, C. Y. (2013). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Applied Psychological Measurement, 37}(8), 598-618.
#' \doi{10.1177/0146621613488436}
#' @export
#'
#' @examples
#' \dontrun{
#' ## Generate data
#' library(GDINA)
#' N = 500
#' Q = sim30GDINA$simQ
#' J = nrow(Q)
#' K= ncol(Q)
#' gs = data.frame(guess = rep(0.2,J), slip = rep(0.2,J))
#' sim = simGDINA(N, Q, gs.parm = gs, model = "DINA")
#' Y = extract(sim,what = "dat")
#' 
#' ## Randomly generate a misspecified Q with 20% of misspecifications
#' mis.Q = matrix(0, J, K)
#' while (any(rowSums(mis.Q)==0)==T){
#'   mis.q = sample(J*K, J*K*0.2) ## percentage of misspecified q
#'   ind = arrayInd(mis.q, dim(Q))
#'   mis.Q = Q
#'   mis.Q[ind] = 1-mis.Q[ind]
#' }
#' 
#' ## Refine the misspecified Q-matrix
#' ref = QR(Y, mis.Q)
#' ref.Q = ref$modified.Q
#' 
#' ## Compute the entry-wise and item-wise recovery rates
#' rr = RR(ref.Q, Q)
#' rr$entry.wise
#' rr$item.wise
#' 
#' ## Compute the retention rate
#' retention.rate(ref.Q, mis.Q, Q)
#' 
#' ## Compute the correction rate
#' correction.rate(ref.Q, mis.Q, Q)
#' }

QR <- function(Y, Q, gate=c("AND", "OR"), max.ite=50)
{
  gate <- match.arg(gate)
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  initial.Q <- Q
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2 ^ natt
  index <- NULL
  RSS <- NULL

  check <- NULL
  check <- CheckInput(Y, Q)
  if (!is.null(check)) return(warning(check))

  pattern <- diag(natt)
  for (l in 2:natt) {
    pattern <- rbind(pattern, t(apply(utils::combn(natt, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  pattern.Q <- pattern[-1,]

  # We may need to run the searching a couple of times until all RSS of all items are stationary.

  for (m in 1:max.ite)
  {
    message(paste("Iteration:", m))
    message("\n")

    max.rss.item <- NULL

    for (k in 1:nitem)
    {
      classification <- NPC(Y, Q, distance = "hamming", gate = "AND")
      est.class <- classification$est.class
      est.ideal <- classification$est.ideal

      # Record the initial classification

      if (m==1 & k==1) initial.class <- est.class

      # Update Q by minimizing RSS

      diff <- Y - est.ideal
      rss <- apply(diff ^ 2, 2, sum)
      RSS <- rbind(RSS, c(rss, sum(rss)))

      # Start the algorithm with the item of largest RSS. Each item is visited only once.

      max.rss <- {if (k == 1) which(rss == max(rss)) else which(rss == max(rss[-max.rss.item]))}
      if (length(max.rss) > 1) max.rss <- sample(max.rss, 1)
      max.rss.item <- c(max.rss.item, max.rss)

      # Find the q-vector among the 2^K possible vectors that the corresponding ideal responses yield
      # minimal RSS with data of item of max RSS.

      update.rss <- NULL

      for (i in 1:(M - 1))
      {
        u <- apply(t(pattern[est.class, ]) ^ pattern.Q[i, ], 2, prod)
        temp.rss <- sum((Y[, max.rss] - u) ^ 2)
        update.rss <- c(update.rss, temp.rss)
      }

      min.update.rss <- which(update.rss == min(update.rss))
      if (length(min.update.rss) > 1) min.update.rss <- sample(min.update.rss, 1)
      update.q <- pattern.Q[min.update.rss, ]
      Q[max.rss, ] <- update.q
    }

    # Stopping criterion: if all of the RSSs between two iterations are identical (which means
    # it's not possible to improve), break the loop.

    if (sum(abs(RSS[((m - 1) * nitem + 1):(m * nitem), nitem + 1] - rep(RSS[((m - 1) * nitem + 1), nitem + 1], nitem))) == 0) break else m <- m + 1
  }

  # Record the terminal classification

  terminal.class <- est.class

  # Report the modified entries

  if (sum((Q - initial.Q) ^ 2) == 0) modified <- "NA" else
  {
    modified <- which((Q - initial.Q) != 0, arr.ind=TRUE)
    colnames(modified) <- c("Item", "Attribute")
  }

  output <- list(modified.Q=Q, modified.entries=modified, initial.class=initial.class, terminal.class=terminal.class)
  class(output) <- "Qrefine"
  return(output)
}

#' @title Print Summary of the Q-Matrix Refinement Result
#' @description Prints a summary of the Q-matrix refinement process, including 
#' which entries were modified and the final refined Q-matrix.
#' @param x An object of class \code{Qrefine}.
#' @param ... Additional arguments passed to print methods.
#' @method print Qrefine
#' @rdname print.Qrefine
#' @export
print.Qrefine <- function(x, ...) {
  
  # 1. Header
  cat("\nQ-Matrix Refinement\n")
  cat("================================\n")
  
  # 2. Check for modifications
  # We check if 'modified.entries' is the string "NA" 
  # (Note: In your QR code, you assign the string "NA", not the logical NA)
  no_changes <- (is.character(x$modified.entries) && x$modified.entries == "NA")
  
  if (no_changes) {
    cat("\nStatus: Converged. No modifications were suggested.\n")
  } else {
    n_changes <- nrow(x$modified.entries)
    cat(paste0("\nStatus: Converged. ", n_changes, " entries were modified.\n"))
    
    cat("\nModified Entries (Item, Attribute):\n")
    print(x$modified.entries)
  }
  
  # 3. Print the Final Matrix
  # We use 'print' specifically for matrices to keep the row/col names aligned
  cat("\nModified Q-Matrix:\n")
  print(x$modified.Q)
  
  # 4. Invisible Return
  # Standard R practice: print methods should return the object invisibly
  invisible(x)
}