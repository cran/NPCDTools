#' @title Check the completeness status of a binary Q-matrix
#'
#' @description
#' \code{Q.completeness} is used to examine whether a given Q-matrix is 
#' complete when data conform to a specified CDM. A Q-matrix is said 
#' to be complete if it allows for the unique identification of all possible attribute profiles 
#' among examinees. So far, the function can only be used for a binary Q-matrix with binary responses.
#'
#' @details
#' The conditions for one Q-matrix completeness are model-dependent: a Q-matrix may be complete for one CDM
#' but incomplete for another. This function implements the theoretical work developed by
#' Chiu et al. (2009) and \enc{Köhn}{Koehn} and Chiu (2017).
#' 
#' \strong{For DINA and DINO models:}
#' 
#' A Q-matrix is complete if and only if it contains all \eqn{K} single-attribute items
#' (Chiu et al., 2009). 
#'
#' \strong{For More General CDMs:}
#' 
#' The function implements a sequential procedure based on Theorems 3-4 and
#' Propositions 1-2 in the work by \enc{Köhn}{Koehn} and Chiu (2017).
#' \enumerate{
#'   \item If Q contains all \eqn{K} single-attribute items, it is complete (Proposition 1).
#'   \item If Q has rank \eqn{< K}, it is incomplete (Theorem 3).
#'   \item For full-rank Q-matrices without all single-attribute items, the function
#'         examines non-nested attribute pairs using indicator vectors to determine
#'         if distinct expected response patterns \eqn{\mathbf{S}(\boldsymbol{\alpha})} can be guaranteed.
#' }
#'
#' The theoretical framework establishes the sufficient conditions for Q completeness, 
#' which means completeness implies distinct expected item response patterns for 
#' all \eqn{2^K} possible attribute profiles.
#'
#' @param raw.Q The Q-matrix that is to be checked, where \eqn{J} is the number
#' of items and \eqn{K} is the number of attributes. It must be a binary (0/1) matrix
#' or data frame that can be coerced to a matrix.
#' @param model Character string specifying the cognitive diagnosis model. Valid options
#'   are "\code{DINA}", "\code{DINO}", or "\code{General}" (for general CDMs including
#'   G-DINA, LCDM, etc.). Case-insensitive. If not specified or invalid, the default is
#'   "\code{General}" with a warning.
#'   
#' @return A list of class "\code{Qcompleteness}" containing:
#'   \item{is_complete}{Logical value indicating completeness: \code{TRUE} if complete,
#'     \code{FALSE} if incomplete, \code{NA} if uncertain.}
#'   \item{status}{Character string: "\code{complete}", "\code{incomplete}", or 
#'     "\code{uncertain}".}
#'   \item{message}{Character string with detailed explanation of the result.}
#'   \item{model}{The CDM used for assessment.}
#'   \item{K}{Number of attributes in the Q-matrix.}
#'   \item{J}{Number of items in the Q-matrix.}
#'   
#'   The function also prints the status message to the console as a side effect.
#'
#' @references
#' Chiu, C.-Y., Douglas, J. A., & Li, X. (2009). Cluster analysis for cognitive diagnosis:
#' Theory and applications. \emph{Psychometrika, 74}(4), 633-665.
#' \doi{10.1007/s11336-009-9125-0}
#' 
#' Köhn, H.-F., & Chiu, C.-Y. (2017). A procedure for assessing the completeness of the
#' Q-matrices of cognitively diagnostic tests. \emph{Psychometrika, 82}(1), 112-132.
#' \doi{10.1007/s11336-016-9536-7}
#' 
#' Köhn, H.-F., & Chiu, C.-Y. (2018). How to build a complete Q-matrix for a 
#' #' cognitively diagnostic test. \emph{Journal of Classification, 35}(2), 273-299.
#' \doi{10.1007/s00357-018-9255-0}
#'
#' @examples
#' \dontrun{
#' # Example 1: Complete Q-matrix for DINA model
#' # (contains all 3 single-attribute items)
#' Q1 <- matrix(c(1, 0, 0,
#'                0, 1, 0,
#'                0, 0, 1,
#'                1, 1, 0,
#'                1, 0, 1), ncol = 3, byrow = TRUE)
#' result1 <- Q.completeness(Q1, model = "DINA")
#' print(result1$is_complete)  # TRUE
#'
#' # Example 2: Incomplete Q-matrix for DINA model
#' # (missing single-attribute items)
#' Q2 <- matrix(c(1, 1, 0,
#'                1, 0, 1,
#'                0, 1, 1), ncol = 3, byrow = TRUE)
#' result2 <- Q.completeness(Q2, model = "DINA")
#' print(result2$is_complete)  # FALSE
#'
#' # Example 3: Check completeness for general CDM
#' Q3 <- matrix(c(1, 0, 0,
#'                0, 1, 0,
#'                0, 0, 1,
#'                1, 1, 0,
#'                1, 0, 1,
#'                0, 1, 1), ncol = 3, byrow = TRUE)
#' result3 <- Q.completeness(Q3, model = "General")
#' }
#' @export

Q.completeness <- function(raw.Q, model = NULL) {
  
  # ==================== INPUT VALIDATION ====================
  
  # Check if raw.Q is provided
  if (missing(raw.Q)) {
    stop("Argument 'raw.Q' is required.", call. = FALSE)
  }
  
  # Check if raw.Q is a matrix or can be converted to one
  if (!is.matrix(raw.Q) && !is.data.frame(raw.Q)) {
    stop("Argument 'raw.Q' must be a matrix or data.frame.", call. = FALSE)
  }
  
  # Convert to matrix
  raw.Q <- as.matrix(raw.Q)
  
  # Check for empty matrix
  if (nrow(raw.Q) == 0L || ncol(raw.Q) == 0L) {
    stop("Q-matrix must have at least one row and one column.", call. = FALSE)
  }
  
  # Check for binary values (0 or 1 only)
  if (!all(raw.Q %in% c(0, 1))) {
    stop("Q-matrix must contain only binary values (0 or 1).", call. = FALSE)
  }
  
  # Check for all-zero columns (attributes never required)
  zero_cols <- which(colSums(raw.Q) == 0)
  if (length(zero_cols) > 0) {
    stop("Q-matrix contains attribute(s) that are never required by any item (all-zero column(s): ",
         paste(zero_cols, collapse = ", "), "). Please remove these attributes.", call. = FALSE)
  }
  
  # --------- model argument handling ----------
  if (missing(model) || is.null(model) || length(model) == 0L) {
    warning(
      "\nArgument 'model' not supplied; defaulting to \"GENERAL\". Valid options are \"DINA\", \"DINO\", or \"GENERAL\".",
      call. = FALSE
    )
    model <- "GENERAL"
  }
  
  model <- toupper(as.character(model))
  if (length(model) > 1L) {
    warning("\nArgument 'model' has length > 1; using the first value only.", call. = FALSE)
    model <- model[1L]
  }
  
  if (!model %in% c("DINA", "DINO", "GENERAL")) {
    warning(
      "\nUnrecognized 'model'; defaulting to \"GENERAL\". Valid options are \"DINA\", \"DINO\", or \"GENERAL\".",
      call. = FALSE
    )
    model <- "GENERAL"
  }
  
  # Store dimensions
  J <- nrow(raw.Q)
  K <- ncol(raw.Q)
  
  # ==================== HELPER FUNCTION ====================
  create_result <- function(is_complete, status, msg, model, K, J) {
    result <- list(
      is_complete = is_complete,
      status = status,
      message = msg,
      model = model,
      K = K,
      J = J
    )
    class(result) <- "Qcompleteness"
    message(msg)  # Print message for user
    return(invisible(result))
  }
  
  # ======================== DINA / DINO branch ========================
  if (model %in% c("DINA", "DINO")) {
    if (is.null(K) || K == 0L) {
      return(create_result(FALSE, "incomplete", "Q is incomplete.", model, K, J))
    }
    
    units <- diag(1, K)
    have_each <- sapply(1:K, function(k) any(apply(raw.Q, 1L, function(r) all(r == units[k, ]))))
    
    if (all(have_each)) {
      return(create_result(TRUE, "complete", "Q is complete.", model, K, J))
    } else {
      return(create_result(FALSE, "incomplete", "Q is incomplete.", model, K, J))
    }
  }
  
  # ======================== General CDMs branch ========================
  M <- 2^K
  
  pattern <- diag(K)
  for (l in 2:K) {
    pattern <- rbind(pattern, t(apply(utils::combn(K, l), 2, function(x) {
      apply(pattern[x, , drop = FALSE], 2, sum)
    })))
  }
  pattern <- rbind(0, pattern)
  
  non.nested <- NULL
  for (m in 2:(M-2)) {
    ind <- NULL
    for (mm in (m+1):(M-1))
      if (length(which(which((pattern[m, ]==1)==TRUE) %in% which((pattern[mm, ]==1)==TRUE))==TRUE) != sum(pattern[m, ]))
        ind <- rbind(ind, c(m, mm))
    non.nested <- rbind(non.nested, ind)
  }
  
  Q <- unique(raw.Q)
  row <- NULL
  for (i in 1:nrow(Q)) {
    for (m in 1:M) {
      if (sum(abs(Q[i, ] - pattern[m, ])) == 0) {
        row <- c(row, m)
        break
      }
    }
  }
  Q <- pattern[sort(row), ]
  
  if (Matrix::rankMatrix(Q) != K) {
    return(create_result(FALSE, "incomplete", "Q is incomplete.", model, K, J))
  } else {
    V <- NULL
    for (m in 1:M) {
      v <- NULL
      for (mm in 1:nrow(Q)) {
        v.temp <- pattern[m, ] * Q[mm, ]
        vv <- matrix(rep(v.temp, M), M, K, byrow = TRUE)
        v <- c(v, which(apply(abs(vv - pattern), 1, sum) == 0))
      }
      V <- rbind(V, v)
    }
    
    sub.ind <- t(utils::combn(nrow(Q), K))
    f <- 0
    for (j in 1:nrow(sub.ind)) {
      if (Matrix::rankMatrix(Q[sub.ind[j, ], ]) != K) j <- j + 1 else {
        sub.V <- V[, sub.ind[j, ]]
        for (i in 1:nrow(non.nested)) {
          non <- NULL
          for (m in 1:ncol(sub.V)) {
            if (sub.V[non.nested[i, 1], m] == sub.V[non.nested[i, 2], m]) {
              non <- c(non, m)
            } else if (min(sub.V[non.nested[i, ], m]) != 1 && max(sub.V[non.nested[i, ], m]) != M) {
              sub.non.nested <- non.nested[which(non.nested[, 1] == min(sub.V[non.nested[i, ], m])), ]
              if (is.matrix(sub.non.nested) == FALSE) {
                match.non <- match(sub.non.nested, sort(sub.V[non.nested[i, ], m]), nomatch = 0)
                non.ind <- which(match.non[2] == 2)
              } else {
                match.non <- matrix(match(sub.non.nested, sort(sub.V[non.nested[i, ], m]), nomatch = 0),
                                    nrow(sub.non.nested), ncol(sub.non.nested))
                non.ind <- which(match.non[, 2] == 2)
              }
              if (length(non.ind) != 0) non <- c(non, m)
            }
          }
          if (length(non) < K) i <- i + 1 else { f <- f + 1; break }
        }
        if (f == 0) { 
          return(create_result(TRUE, "complete", "Q is complete.", model, K, J))
          break 
        } else j <- j + 1
      }
    }
    if (f == 0 && j == (nrow(sub.ind) + 1)) {
      return(create_result(FALSE, "incomplete", "Q is incomplete.", model, K, J))
    }
    if (f != 0 && j == (nrow(sub.ind) + 1)) {
      return(create_result(NA, "uncertain", "Q may be incomplete.", model, K, J))
    }
  }
}

#' @title Print Summary of a Qcompleteness Object
#'
#' @description
#' Print method for objects of class \code{"Qcompleteness"}.
#'
#' @param x An object of class \code{"Qcompleteness"}
#' @param ... Additional arguments (not used)
#'
#' @method print Qcompleteness
#' @rdname print.Q.completeness
#' @export
print.Qcompleteness <- function(x, ...) {
  cat("\nQ-Matrix Completeness Check\n")
  cat("===========================\n\n")
  cat("Model:      ", x$model, "\n", sep = "")
  cat("Dimensions: J = ", x$J, " items, K = ", x$K, " attributes\n", sep = "")
  cat("Result: ", x$message, "\n\n", sep = "")
  invisible(x)
}