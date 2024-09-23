#' @title Refine the Q-matrix by Minimizing the RSS
#'
#' @description We estimate memberships using the non-parametric classification
#'     method (weighted hamming), and comparisons of the residual sum of squares
#'     computed from the observed and the ideal item responses.
#'
#' @section The Q-Matrix Refinment (QR) Method:
#'
#' This function implements the Q-matrix refinement method developed by Chiu
#' (2013), which is also based on the aforementioned nonparametric classification
#' methods (Chiu & Douglas, 2013). This Q-matrix refinement method corrects
#' potential misspecified entries of the Q-matrix through comparisons of the
#' residual sum of squares computed from the observed and the ideal item responses.
#'
#' The algorithm operates by minimizing the RSS. Recall that \eqn{Y_{ij}} is the
#' observed response and \eqn{\eta_{ij}} is the ideal response.
#' Then the RSS of item \eqn{j} for examinee \eqn{i} is defined as
#' \deqn{RSS_{ij} = (Y_{ij} - \eta_{ij})^2}.
#' The RSS of item \eqn{j} across all examinees is therefor
#' \deqn{RSS_{j} = \sum_{i=1}^{N} (Y_{ij} - \eta_{ij})^2 = \sum_{m=1}^{2^k} \sum_{i \in C_{m}} (Y_{ij} - \eta_{jm})^2}
#' where \eqn{C_m} is the latent proficiency-class \eqn{m}, and \eqn{N}
#' is the number of examinees. Chiu(2013) proved that the expectation of
#' \eqn{RSS_j} is minimized for the correct q-vector among the
#' \eqn{2^K - 1} candidates. Please see the paper for the justification.
#'
#'
#' @param Y A matrix of binary responses (1=correct, 0=incorrect). Rows
#'     represent persons and columns represent items.
#' @param Q The Q-matrix of the test. Rows represent items and columns
#'     represent attributes.
#' @param gate A string, "AND" or "OR". "AND": the examinee needs to possess
#'     all related attributes to answer an item correctly.
#'     "OR": the examinee needs to possess only one of the related attributes
#'     to answer an item correctly.
#' @param max.ite The number of iterations to run until all RSS of all items
#'     are stationary.
#'
#' @return A list containing:
#' \item{initial.class}{Initial classification}
#' \item{terminal.class}{Terminal classification}
#' \item{modified.Q}{The modified Q-matrix}
#' \item{modified.entries}{The modified q-entries}
#'
#'
#' @references
#' Chiu, C. Y. (2013). Statistical Refinement of the Q-matrix in Cognitive Diagnosis. \emph{Applied Psychological Measurement, 37(8)}, 598-618.
#'
#' @export
#'
#'

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

  pattern <- AlphaPermute(natt)
  pattern.Q <- pattern[-1,]

  # We may need to run the searching a couple of times until all RSS of all items are stationary.

  for (m in 1:max.ite)
  {
    message(paste("Iteration:", m))
    message("\n")

    max.rss.item <- NULL

    for (k in 1:nitem)
    {
      classification <- NPCD::AlphaNP(Y, Q, gate, method="Hamming")
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

  output <- list(modified.Q=Q, modified.entries=modified, initial.class=initial.class, terminal.class=terminal.class, patterns=pattern, initial.Q=initial.Q, Y=Y, gate=gate)
  class(output) <- "Qrefine"
  return(output)
}
