#' @title Estimation of examinees' attribute profiles using the GNPC method
#'
#' @description
#' Function \code{GNPC} is used to estimate examinees' attribute profiles using
#'     the general nonparametric classification (GNPC) method
#'     (Chiu et al., 2018; Chiu & Köhn, 2019). It can be
#'     used with data conforming to any cognitive diagnosis models (CDMs).
#'
#' @section Details:
#' A weighted ideal response \eqn{\eta^{(w)}}, defined as the convex combination
#' of \eqn{\eta^{(c)}} and \eqn{\eta^{(d)}}, is used in the GNPC method to compute distances.
#' Suppose item \eqn{j} requires \eqn{K_{j}^* \leq {K}} attributes that, without loss of
#' generality, have been moved to the first \eqn{K_{j}^*} positions of the item
#' attribute vector \eqn{\boldsymbol{q_j}}. For each item \eqn{j} and latent class \eqn{\mathcal{C}_{l}},
#' the weighted ideal response \eqn{\eta_{lj}^{(w)}} is defined as the convex combination
#' \eqn{\eta_{lj}^{(w)} = w _{lj} \eta_{lj}^{(c)}+(1-w_{lj})\eta_{lj}^{(d)}}
#' where \eqn{0\leq w_{lj}\leq 1}. The distance between the observed responses
#' to item \eqn{j} and the weighted ideal responses \eqn{w_{lj}^{(w)}} of examinees
#' in \eqn{\mathcal{C}_{l}} is defined as the sum of squared deviations:
#' \eqn{d_{lj} = \sum_{i \in \mathcal {C}_{l}} (y_{ij} - \eta_{lj}^{(w)})^2}. 
#' \eqn{\hat{w}_{lj}} can then be obtained by minimizing \eqn{d_{lj}}, which can then be used to compute \eqn{\hat{\eta}_{lj}}. 
#'
#' After all the \eqn{\hat{\eta}_{lj}} are obtained, examinees' attribute profiles \eqn{\boldsymbol{\alpha}} 
#' can be estimated by minimizing the loss function \eqn{\hat{d}_{lj} = \sum_{i \in \mathcal{C}_{l}} (y_{ij} - \hat{\eta}_{lj}^{(w)})^2}
#' 
#' The algorithm iteratively updates the weighted ideal responses and reclassifies
#' examinees until convergence is achieved. The stopping criterion is based on the proportion
#' of examinees whose classification changes between consecutive iterations:
#' \eqn{\frac{\sum_{i=1}^{N} I\left[\alpha_i^{(t)} \neq \alpha_i^{(t-1)}\right]}{N} < \epsilon}{
#'      sum_{i=1}^N I[alpha_i^(t) != alpha_i^(t-1)]/N < epsilon
#' }
#' where \eqn{\epsilon}{} is the tolerance level (default = 0.001).
#' 
#' 
#' The default initial values of \eqn{\boldsymbol{\alpha}} are obtained by using the NPC method. Chiu et al. (2018)
#' suggested another viable alternative for obtaining initial estimates of the proficiency classes by
#' using an ideal response with fixed weights defined as
#' \eqn{\eta_{lj}^{(fw)}=\frac{\sum_{k=1}^{K}\alpha_{k}q_{jk}}{K}\eta_{lj}^{(c)}+(1-\frac{\sum_{k=1}^{K}\alpha_{k}q_{jk}}{K})\eta_{lj}^{(d)}}.
#'
#'
#' @param Y A \eqn{N \times J} binary data matrix consisting of the responses
#'     from \eqn{N} examinees to \eqn{J} items.
#' @param Q A \eqn{J \times K} binary Q-matrix where the entry \eqn{q_{jk}}
#'     describing whether the \eqn{k}th attribute is required by the \eqn{j}th item.
#' @param fixed.w \code{TRUE} or \code{FALSE}. When \code{TRUE} is specified, fixed weights 
#'     are used as the initial weights to compute fixed-weight ideal response and there is no 
#'     need to use NPC to produce the initial values for examinees' attribute profiles. 
#'     Hence, \code{initial.dis} and \code{initial.gate} are turned off. The default value is \code{FALSE}.
#' @param initial.dis The type of distance used in the \code{NPC} to carry
#'     out the initial attribute profiles for the GNPC method.
#'     Allowable options are \code{"hamming"} and \code{"whamming"} representing
#'     the Hamming and the weighted Hamming distances, respectively.
#' @param initial.gate The type of relation between examinees' attribute profiles
#'     and the items.
#'     Allowable relations are "\code{AND}" and "\code{OR}",
#'     representing the conjunctive and disjunctive relations, respectively.
#' @param max.iter Maximum number of iterations allowed. Default is 1000.
#' @param tol Convergence tolerance. The algorithm stops when the proportion of
#'     examinees whose classification changes is less than this value. Default is 0.001.
#' @param track.convergence Logical. If \code{TRUE}, convergence information is tracked
#'     and returned for diagnostic purposes. Default is \code{TRUE}.
#'
#' @return The function returns a list with the following components:
#' \describe{
#' \item{att.est}{A \eqn{N \times K} matrix of estimated attribute profiles for examinees}
#' \item{class}{A vector of length \eqn{N} containing the estimated class memberships}
#' \item{ideal.response}{A \eqn{2^K \times J} matrix of weighted ideal responses}
#' \item{weight}{A \eqn{2^K \times J} matrix of weights used to compute the weighted ideal responses}
#' \item{convergence}{(Only if \code{track.convergence = TRUE}) A list containing:
#'   \itemize{
#'     \item \code{iteration}: Vector of iteration numbers
#'     \item \code{prop.change}: Proportion of examinees whose classification changed at each iteration
#'     \item \code{total.distance}: Total squared distance between observed and weighted ideal responses
#'     \item \code{n.iter}: Total number of iterations until convergence
#'     \item \code{converged}: Logical indicating whether the algorithm converged within \code{max.iter}
#'   }
#' }
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage
#' library(GDINA)
#' set.seed(123)
#' N <- 500
#' Q <- sim30GDINA$simQ
#' gs <- data.frame(guess = rep(0.2, nrow(Q)), slip = rep(0.2, nrow(Q)))
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "DINA")
#' Y <- extract(sim, what = "dat")
#' alpha <- extract(sim, what = "attribute")
#' 
#' # Analyze data using GNPC
#' result <- GNPC(Y, Q, initial.dis = "hamming", initial.gate = "AND")
#' 
#' # View results
#' head(result$att.est)
#' table(result$class)
#' 
#' # Plot overall convergence 
#' plot(result)
#' 
#' # Plot individual examinee's convergence
#' plot(result, type = "individual", examinee.id = 1, true.alpha = alpha[1, ])
#' 
#' # Check attribute agreement rate
#' PAR(alpha, result$att.est)
#' AAR(alpha, result$att.est)
#' 
#' # Example 2: Without convergence tracking (Convergence tracking is only used for the GNPC plots.)
#' result2 <- GNPC(Y, Q, track.convergence = FALSE)
#' }
#'
#'
#' @references
# 'Chiu, C.-Y., Sun, Y., & Bian, Y. (2018). Cognitive diagnosis for small educational
#'   programs: The general nonparametric classification method. \emph{Psychometrika},
#'   \emph{83}(2), 355--375.
#'   \doi{10.1007/s11336-017-9595-4}
#'
# 'Chiu, C.-Y., & Köhn, H.-F. (2019). Consistency theory for the general nonparametric
#'   classification method. \emph{Psychometrika}, \emph{84}(3), 830--845.
#' \doi{10.1007/s11336-019-09660-x}
#'
#' @export
GNPC <- function(Y, Q, 
                 fixed.w = FALSE,
                 initial.dis = "hamming", 
                 initial.gate = "AND",
                 max.iter = 1000,
                 tol = 0.001,
                 track.convergence = TRUE) {
  
  initial.dis  <- tolower(initial.dis)
  initial.gate <- toupper(initial.gate)
  initial.dis  <- match.arg(initial.dis,  c("hamming", "whamming"))
  initial.gate <- match.arg(initial.gate, c("AND", "OR"))
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)
  if (!is.null(check)) stop(check)
  
  N <- dim(Y)[1]
  K <- dim(Q)[2]
  J <- dim(Q)[1]
  M <- 2^K
  
  # Generate all possible attribute patterns
  pattern <- diag(K)
  for (l in 2:K) {
    pattern <- rbind(pattern, t(apply(utils::combn(K, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  
  # Conjunctive Ideal Response
  Ideal <- pattern %*% t(Q)
  Ideal.conj <- 1 * (Ideal == (matrix(1, M) %*% t(rowSums(Q))))
  
  # Disjunctive Ideal Response
  Ideal.dis <- 1 * (Ideal >= 1)
  
  # Use fixed weight Ideal Response as the initial of the weight
  if (fixed.w == TRUE){
    weight <- Ideal / matrix(rep(colSums(t(Q)), M), M, J, TRUE)
    Ideal.weight <- Ideal.conj + (Ideal.dis - Ideal.conj) * weight
    
    w.class <- NULL
    for (i in 1:N) {
      ham <- rowSums(((matrix(Y[i, ], M, J, byrow = TRUE) - Ideal.weight))^2)
      min.ham <- which(ham == min(ham))
      if (length(min.ham) != 1) {
        min.ham <- sample(min.ham, 1, prob = rep(1 / length(min.ham), length(min.ham)))
      }
      w.class <- c(w.class, min.ham)
    }
    initial.class <- w.class
  } else {
    # If fixed weight is not used, obtain initial attribute profiles using NPC
    initial.gate <- match.arg(initial.gate)
    initial.dis <- match.arg(initial.dis)
    initial.class = NPC(Y, Q, distance = initial.dis, gate = initial.gate)$est.class
  }
  
  # Initialize convergence tracking
  if (track.convergence) {
    conv.prop.change <- numeric(0)
    conv.total.dist <- numeric(0)
    conv.class.size <- matrix(NA, nrow = 0, ncol = M)
    conv.class.avg.dist <- matrix(NA, nrow = 0, ncol = M)
    conv.ideal.response <- list()  # Store weighted ideal at each iteration
    conv.classifications <- matrix(NA, nrow = 0, ncol = N)  # Store classifications
  }
  
  # Iteration Starts
  d <- 1
  time <- 0
  converged <- FALSE
  
  while (d > tol && time < max.iter) {
    time <- time + 1
    
    # Compute the general weights using the closed form
    w <- NULL
    temp <- matrix(NA, M, 1)
    Ideal.comb <- Ideal.dis + Ideal.conj
    
    for (j in 1:J) {
      pQ <- pattern * matrix(Q[j, ] == 1, M, K, TRUE)
      upQ <- unique(pQ)
      ng <- nrow(upQ)
      
      for (g in 1:ng) {
        match <- apply(pQ, 1, identical, upQ[g, ])
        m <- which(match)
        c <- which(initial.class %in% m)
        c <- c[!is.na(c)]
        
        if (length(c) == 0) {
          temp[m, ] <- 0.01
        } else if (length(c) != 0) {
          temp[m, ] <- sum(Y[c, j]) / sum((1 - Ideal.conj[initial.class[c], j])^2)
        }
      }
      temp[which(Ideal.comb[, j] == 0)] <- 0.01
      temp[which(Ideal.comb[, j] == 2)] <- 0.99
      w <- cbind(w, temp)
    }
    
    Ideal.w <- (1 - w) * Ideal.conj + w * Ideal.dis
    
    w.class <- NULL
    for (i in 1:N) {
      ham <- rowSums(((matrix(Y[i, ], M, J, byrow = TRUE) - Ideal.w))^2)
      min.ham <- which(ham == min(ham))
      if (length(min.ham) != 1) {
        min.ham <- sample(min.ham, 1, prob = rep(1 / length(min.ham), length(min.ham)))
      }
      w.class <- c(w.class, min.ham)
    }
    
    # Calculate proportion of changed classifications
    d <- length(which(w.class - initial.class != 0)) / N
    
    # Track convergence metrics if requested
    if (track.convergence) {
      conv.prop.change <- c(conv.prop.change, d)
      
      # Calculate total squared distance
      total.dist <- 0
      for (i in 1:N) {
        total.dist <- total.dist + sum((Y[i, ] - Ideal.w[w.class[i], ])^2)
      }
      conv.total.dist <- c(conv.total.dist, total.dist)
      
      # Track class-level metrics
      class.size.vec <- numeric(M)
      class.avg.dist.vec <- numeric(M)
      
      for (m in 1:M) {
        members <- which(w.class == m)
        class.size.vec[m] <- length(members)
        
        if (length(members) > 0) {
          class.distances <- numeric(length(members))
          for (idx in seq_along(members)) {
            i <- members[idx]
            class.distances[idx] <- sum((Y[i, ] - Ideal.w[m, ])^2)
          }
          class.avg.dist.vec[m] <- mean(class.distances)
        } else {
          class.avg.dist.vec[m] <- NA
        }
      }
      
      conv.class.size <- rbind(conv.class.size, class.size.vec)
      conv.class.avg.dist <- rbind(conv.class.avg.dist, class.avg.dist.vec)
      
      # Store weighted ideal responses for this iteration
      conv.ideal.response[[time]] <- Ideal.w
      
      # Store classifications
      conv.classifications <- rbind(conv.classifications, w.class)
    }
    
    initial.class <- w.class
  }
  
  if (d <= tol) {
    converged <- TRUE
  }
  
  att.est <- pattern[w.class, ]
  
  # Build output list - always include Y for individual plotting
  output <- list(
    att.est = att.est, 
    class = w.class, 
    ideal.response = Ideal.w, 
    weight = w,
    Y = Y  # Store response data
  )
  
  # Add convergence information if tracked
  if (track.convergence) {
    class.labels <- apply(pattern, 1, function(x) paste0("(", paste(x, collapse = ","), ")"))
    colnames(conv.class.size) <- class.labels
    colnames(conv.class.avg.dist) <- class.labels
    
    output$convergence <- list(
      iteration = 1:time,
      prop.change = conv.prop.change,
      total.distance = conv.total.dist,
      class.size = conv.class.size,
      class.avg.distance = conv.class.avg.dist,
      attribute.patterns = pattern,
      ideal.response.history = conv.ideal.response,  # NEW: weighted ideal at each iteration
      classifications = conv.classifications,  # NEW: classifications at each iteration
      n.iter = time,
      converged = converged
    )
  }
  
  class(output) <- "GNPC"
  return(output)
}


#' @title Print the Summary of a GNPC Object
#'
#' @description
#' Prints a summary of the GNPC estimation results.
#'
#' @param x An object of class \code{GNPC}.
#' @param ... Additional arguments (not used).
#'
#'
#' @method print GNPC
#' @rdname print.GNPC
#' @export
print.GNPC <- function(x, ...) {
  cat("GNPC Estimation Results\n")
  cat("=======================\n\n")
  
  N <- nrow(x$att.est)
  K <- ncol(x$att.est)
  J <- ncol(x$ideal.response)
  M <- nrow(x$ideal.response)
  
  cat(sprintf("Number of examinees: %d\n", N))
  cat(sprintf("Number of items: %d\n", J))
  cat(sprintf("Number of attributes: %d\n", K))
  cat(sprintf("Number of possible attribute patterns: %d\n", M))
  cat("\n")
  
  # Print convergence information if available
  if (!is.null(x$convergence)) {
    cat("Convergence Information:\n")
    cat(sprintf("  Number of iterations: %d\n", x$convergence$n.iter))
    cat(sprintf("  Converged: %s\n", ifelse(x$convergence$converged, "Yes", "No")))
    cat(sprintf("  Final proportion of changes: %.4f\n", 
                x$convergence$prop.change[length(x$convergence$prop.change)]))
    cat("\n")
  }
  
  cat("Use plot() to visualize convergence.\n")
  cat("Access results with: $att.est, $class, $ideal.response, $weight\n")
  
  invisible(x)
}
