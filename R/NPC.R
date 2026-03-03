#' @title Estimation of examinees' attribute profiles using the NPC method
#'
#' @description The function estimates examinees' attribute profiles
#'     using the nonparametric classification (NPC) method (Chiu & Douglas, 2013).
#'     An examinee's attribute profile is estimated by minimizing the distance between
#'     the observed and ideal item responses.
#'
#' @section Details:
#' The nonparametric classification (NPC) method (Chiu & Douglas, 2013) assigns examinees to the
#' proficiency classes they belong to by comparing their observed item response patterns with each of the ideal
#' item response patterns of the \eqn{2^K} proficiency classes. When there is no data perturbation, an 
#' examinee's ideal response pattern corresponding to the examinee's true attribute pattern and his/her 
#' observed item response patterns are identical, and thus the distance between them is 0. When data 
#' perturbations are small, this ideal response pattern remains the one most similar to the observed 
#' response pattern, which is exactly the setup of data conforming to the DINA or DINO model. Hence, based 
#' on this rationale, an examinee's attribute profile is obtained by minimizing the distance between the 
#' observed and the ideal item response patterns. The nonparametric nature of the NPC method furthermore 
#' makes it suitable for data obtained from small-scale settings.
#' 
#' @param Y A \eqn{N \times J} binary data matrix consisting of the responses from \eqn{N} examinees to
#' \eqn{J} items.
#' @param Q A \eqn{J \times K} binary Q-matrix where the entry \eqn{q_{jk}} describes whether the \eqn{k}th attribute
#' is required by the \eqn{j}th item.
#'
#' @param distance The type of distance used to compute the loss function. The possible options include  
#' (i) "\code{hamming}" representing the plain Hamming distance method, 
#' (ii) "\code{whamming}" representing the Hamming distance weighted by the inverse of item variance, and  
#' (iii) "\code{penalized}" representing the Hamming distance weighted by the inverse of item variance
#'       and specified penalizing weights for guess and slip.  
#' @param gate A character string specifying the type of gate. The possible options include "\code{AND}" and "\code{OR}" standing
#' for conjunctive and disjunctive gate, respectively.
#' @param wg Additional argument for the "penalized" method. It is a weight
#' assigned to guesses in the DINA or DINO models. A large
#' value of weight results in a stronger impact on the
#' distance (i.e., larger loss function values) caused by guessing.
#' @param ws Additional input for the "penalized" method.
#'     It is the weight assigned to slips in the DINA or DINO models.
#'     A large value of weight results in a stronger impact on the
#'     distance (i.e., larger loss function values) caused by slipping.
#'
#' @return The function returns a series of outputs, including:
#' \describe{
#'   \item{alpha.est}{A \eqn{N \times K}{}  matrix representing the estimated attribute profiles.
#'       1 = examinee masters the attribute, 0 = examinee does not master the attribute.}
#'   \item{est.ideal}{A \eqn{N \times J}{}  matrix indicating the estimated ideal response to all 
#'   items from all examinees. 1 = correct, 0 = incorrect.}
#'   \item{est.class}{A \eqn{N}{}-dimensional vector showing the class memberships for all examinees.}
#'   \item{n.tie}{The number of ties in the Hamming distance among the candidate
#'       attribute profiles for each person. When ties occur, one of
#'       the tied attribute profiles is randomly chosen.}
#'   \item{pattern}{All possible attribute profiles in the latent space.}
#'   \item{loss.matrix}{A \eqn{2^K \times N}{} matrix containing the values of the loss function 
#'   (the distances) between each examinee's observed response vector and the \eqn{2^K}{} ideal response vectors.}
#' }
#'
#' @references
#' Chiu, C. Y., & Douglas, J. A. (2013). A nonparametric approach to cognitive diagnosis by proximity to ideal response patterns.
#' \emph{Journal of Classification, 30}(2), 225-250.
#' \doi{10.1007/s00357-013-9132-9}
#'
#' @export
#' @seealso \code{\link{GNPC}}
#' @examples
#' \dontrun{
#' library(GDINA)
#' N <- 500
#' Q <- sim30GDINA$simQ
#' gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA")
#' Y <- extract(sim, what = "dat")
#' alpha <- extract(sim, what = "attribute")
#' 
#' # Estimate attribute profiles using NPC
#' result <- NPC(Y, Q, distance = "hamming", gate = "AND")
#' print(result)
#' result$alpha.est
#' 
#' # Check attributed agreement rate
#' PAR(alpha, result$alpha.est)
#' AAR(alpha, result$alpha.est)
#' }


NPC = function (Y, Q,  distance = c("hamming", "whamming", "penalized"), gate = c("AND", "OR"),wg = 1, ws = 1) 
{
  distance <- tolower(distance)
  gate <- toupper(gate)
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)
  if (!is.null(check)) 
    return(warning(check))
  gate <- match.arg(gate)
  distance <- match.arg(distance)
  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2^natt
  pattern <- diag(natt)
  for (l in 2:natt) {
    pattern <- rbind(pattern, t(apply(utils::combn(natt, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  
  Ideal <- matrix(NA, M, nitem)
  for (m in 1:M) {
    for (j in 1:nitem) {
      if (gate == "AND") {
        u <- prod(pattern[m, ]^Q[j, ])
      }
      else if (gate == "OR") {
        u <- 1 - prod((1 - pattern[m, ])^Q[j, ])
      }
      else {
        return(warning("Gate specification not valid."))
      }
      Ideal[m, j] <- u
    }
  }
  if (distance == "hamming") {
    if (ws != 1 | wg != 1) {
      stop("Guessing weight and slipping weight should be 1.")
    }
    else {
      weight <- rep(1, nitem)
    }
  }
  else if (distance == "whamming") {
    if (ws != 1 || wg != 1) {
      stop("Guessing weight and slipping weight should be 1.")
    }
    else {
      p.bar <- apply(Y, 2, mean)
      weight <- 1/(p.bar * (1 - p.bar))
      weight[weight > 1/(0.95 * 0.05)] <- 1/(0.95 * 0.05)
    }
  }
  else if (distance == "penalized") {
    p.bar <- apply(Y, 2, mean)
    weight <- 1/(p.bar * (1 - p.bar))
    weight[weight > 1/(0.95 * 0.05)] <- 1/(0.95 * 0.05)
    if (ws == wg) {
      warning("Penalizing weights for guess and slip are the same --> equivalent with the \"Hamming\" distance.")
    }
  }
  else {
    return(warning("distance specification not valid."))
  }
  loss.matrix <- matrix(NA, nrow = M, ncol = nperson)
  est.class <- NULL
  est.pattern <- NULL
  n.tie <- rep(0, nperson)
  for (i in 1:nperson) {
    Y.matrix <- matrix(rep(Y[i, ], M), M, nitem, byrow = TRUE)
    loss <- apply(matrix(rep(weight, M), M, nitem, byrow = TRUE) * 
                    (wg * abs(Y.matrix - Ideal) * Y.matrix + ws * abs(Y.matrix - 
                                                                        Ideal) * (1 - Y.matrix)), 1, sum)
    loss.matrix[, i] <- loss
    min.loss <- which(loss == min(loss))
    if (length(min.loss) != 1) {
      n.tie[i] <- length(min.loss)
      min.loss <- sample(min.loss, 1, prob = rep(1/length(min.loss), 
                                                 length(min.loss)))
    }
    est.class <- c(est.class, min.loss)
  }
  est.pattern <- pattern[est.class, ]
  est.ideal <- Ideal[est.class, ]
  output <- list(alpha.est = est.pattern, est.ideal = est.ideal, 
                 est.class = est.class, n.tie = n.tie, pattern = pattern, 
                 loss.matrix = loss.matrix, distance = distance,
                 gate = gate, Q = Q, Y = Y)
  class(output) <- "NPC"
  return(output)
}


#' @title Print the Summary of an NPC Object
#' @description Prints a summary of NPC classification results.
#' @param x An object of class \code{"NPC"}.
#' @param ... Additional arguments (not used).
#' @method print NPC
#' @rdname print.NPC
#' @export
print.NPC <- function(x, ...) {
  cat("NPC Classification Results\n")
  cat("=========================\n\n")
  
  # Basic dimensions
  N <- nrow(x$alpha.est)
  K <- ncol(x$alpha.est)
  J <- nrow(x$Q)
  M <- nrow(x$pattern)
  
  cat(sprintf("Number of examinees: %d\n", N))
  cat(sprintf("Number of items: %d\n", J))
  cat(sprintf("Number of attributes: %d\n", K))
  cat(sprintf("Number of possible attribute patterns: %d\n", M))
  
  # Gate (if stored) and distance
  if (!is.null(x$gate)) {
    cat(sprintf("Gate: %s\n", x$gate))
  }
  if (!is.null(x$distance)) {
    cat(sprintf("Distance: %s\n", x$distance))
  }
  cat("\n")
  
  # Class usage summary
  class.freq <- table(x$est.class)
  n.used <- length(class.freq)
  cat("Classification Summary:\n")
  cat(sprintf("  Number of non-empty classes: %d\n", n.used))
  cat(sprintf("  Largest class size: %d\n", max(class.freq)))
  cat(sprintf("  Smallest non-empty class size: %d\n", min(class.freq)))
  cat("\n")
  
  # Tie information
  n.tied <- sum(x$n.tie > 1)
  if (n.tied > 0) {
    cat("Tie Information:\n")
    cat(sprintf("  Examinees with tied minimum loss: %d\n", n.tied))
    cat(sprintf("  Max number of ties for one examinee: %d\n", max(x$n.tie)))
    cat("\n")
  } else {
    cat("Tie Information:\n")
    cat("  No ties in minimum loss classification.\n\n")
  }
  
  cat("Access results with: $alpha.est, $est.class, $est.ideal, $loss.matrix\n")
  invisible(x)
}

