#' @title Estimation of  examinees' attribute profiles using the NPC method
#'
#' @description The function is used to estimate examinees' attribute profiles
#'     using the nonparametric classification (NPC) method (Chiu, & Douglas, 2013).
#'     It uses a distance-based algorithm on the observed item responses for
#'     classifying examiness. This function estimates attribute profiles using
#'     nonparametric approaches for both the "AND gate" (conjunctive) and the
#'     "OR gate" (disjunctive) cognitive diagnostic models. These algorithms select
#'     the attribute profile with the smallest loss function value (plain,
#'     weighted, or penalized Hamming distance, see below for details) as the
#'     estimate. If more than one attribute profiles have the smallest loss
#'     function value, one of them is randomly chosen.
#'
#' @section NPC algorithm with three distacne methods:
#' Proficiency class membership is determined by comparing an examinee's
#' observed item response vector \eqn{\boldsymbol{Y}} with each of the ideal
#' item response vectors of the realizable \eqn{2^K=M} proficiency classes.
#' The ideal item responses are a function of the Q-matrix and the attribute
#' vectors characteristic of the different proficiency classes. Hence, an
#' examinee’s proficiency class is identified by the attribute vector
#' \eqn{\boldsymbol{\alpha}_{m}} underlying that ideal item response vector
#' which is closest—or most similar—to an examinee’s observed item response
#' vector. The ideal response to item j is the score that would be obtained
#' by an examinee if no perturbation occurred.
#'
#' Let \eqn{\boldsymbol{\eta}_{i}} denote the J-dimensional ideal item response
#' vector of examinee i, and the \eqn{\hat{\boldsymbol{\alpha}}} of an
#' examinee’s attribute vector is defined as the attribute vector
#' underlying the ideal item response vector that among all ideal item response
#' vectors minimizes the distance to an examinee’s observed item response vector:
#' \eqn{\hat{\boldsymbol{\alpha}} = \arg\min_{m \in \{1,2,\ldots,M\}} d(\boldsymbol{y_i}, \boldsymbol{\eta}_{m})}
#'
#' A distance measure often used for clustering binary data is the Hamming
#' distance that simply counts the number of disagreements between two vectors:
#' \eqn{d_H(\boldsymbol{y,\eta}) =  \sum_{j=1}^{J} | y_j - \eta_j |}
#'
#' If the different levels of variability in the item responses are to be
#' incorporated, then the Hamming distances can be weighted, for example, by the
#' inverse of the item sample variance, which allows for larger impact on the
#' distance functions of items with smaller variance:
#' \eqn{d_{wH} (\boldsymbol{y,\eta}) = \sum_{j=1}^{J} \frac{1}{\overline{p_j}(1-\overline{p_j})} |y_j-\eta_j|}
#'
#' Weighting weighting differently for departures from the ideal response model
#' that would result from slips versus guesses is also considered:
#' \eqn{d_{gs}(\boldsymbol{y,\eta})=\sum_{j=1}^{J} w_gI[y_j=1]|y_j-\eta_j|+\sum_{j=1}^{J}w_sI[y_j=0]|y_j-\eta_j|}
#'
#' @param Y A matrix of binary responses. Rows represent persons and columns
#'     represent items. 1=correct, 0=incorrect.
#' @param Q The Q-matrix of the test. Rows represent items and columns represent
#'     attributes. 1=attribute required by the item, 0=attribute not required
#'     by the item.
#' @param gate A character string specifying the type of gate. It can be one of
#'     the following:
#' \describe{
#'   \item{"AND"}{The examinee needs to possess all required attributes of an
#'       item in order to answer it correctly.}
#'   \item{"OR"}{The examinee needs to possess only one of the required
#'       attributes of an item in order to answer it correctly.}
#' }
#' @param method The method of nonparametric estimation.
#' \describe{
#'   \item{"Hamming"}{The plain Hamming distance method}
#'   \item{"Weighted"}{The Hamming distance weighted by inversed item variance}
#'   \item{"Penalized"}{The Hamming distance weighted by inversed item variance
#'       and specified penalizing weights for guess and slip.}
#' }
#' @param wg Additional argument for the "penalized" method.
#'     It is the weight assigned to guessing in the DINA or DINO models. A large
#'     value of weight results in a stronger impact on
#'     Hamming distance (larger loss function values) caused by guessing.
#' @param ws Additional input for the "penalized" method.
#'     It is the weight assigned to slipping in the DINA or DINO models.
#'     A large value of weight results in la stronger impact on Hamming
#'     distance (larger loss function values) caused by slipping.
#'
#' @return The function returns a series of outputs, including:
#' \describe{
#'   \item{alpha.est}{Estimated attribute profiles. Rows represent persons
#'       and columns represent attributes.
#'       1=examinee masters the attribute, 0=examinee does not master the attribute.}
#'   \item{est.ideal}{Estimated ideal response to all items by all examinees.
#'       Rows represent persons and columns represent items.
#'       1=correct, 0=incorrect.}
#'   \item{est.class}{The class number (row index in pattern) for each person's
#'       attribute profile.
#'       It can also be used for locating the loss function value in loss.matrix
#'       for the estimated attribute profile for each person.}
#'   \item{n.tie}{Number of ties in the Hamming distance among the candidate
#'       attribute profiles for each person. When we encounter ties, one of
#'       the tied attribute profiles is randomly chosen.}
#'   \item{pattern}{All possible attribute profiles in the search space.}
#'   \item{loss.matrix}{The matrix of the values for the loss function (the
#'       plain, weighted, or penalized Hamming distance).
#'       Rows represent candidate attribute profiles in the same order with
#'       the pattern matrix; columns represent different examinees.}
#' }
#'
#' @references
#' Chiu, C. (2011). Flexible approaches to cognitive diagnosis: nonparametric methods and small sample techniques.
#' Invited session of cognitive diagnosis and item response theory at 2011 Joint Statistical Meeting.
#'
#' Chiu, C. Y., & Douglas, J. A. (2013). A nonparametric approach to cognitive diagnosis by proximity to ideal response patterns.
#' Journal of Classification 30(2), 225-250.
#'
#'
#' @export
#' @examples
#' # Generate item and examinee profiles
#'
#' natt <- 3
#' nitem <- 4
#' nperson <- 5
#' Q <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
#' alpha <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
#'
#' # Generate DINA model-based response data
#' slip <- c(0.1, 0.15, 0.2, 0.25)
#' guess <- c(0.1, 0.15, 0.2, 0.25)
#' my.par <- list(slip=slip, guess=guess)
#' data <- matrix(NA, nperson, nitem)
#' eta <- matrix(NA, nperson, nitem)
#' for (i in 1:nperson) {
#'   for (j in 1:nitem) {
#'     eta[i, j] <- prod(alpha[i,] ^ Q[j, ])
#'     P <- (1 - slip[j]) ^ eta[i, j] * guess[j] ^ (1 - eta[i, j])
#'     u <- runif(1)
#'     data[i, j] <- as.numeric(u < P)
#'     }
#' }
#'
#' # Using the function to estimate examinee attribute profile
#' alpha.est.NP.H <- NPC(data, Q, gate="AND", method="Hamming")
#' alpha.est.NP.W <- NPC(data, Q, gate="AND", method="Weighted")
#' alpha.est.NP.P <- NPC(data, Q, gate="AND", method="Penalized", wg=2, ws=1)
#'
#' nperson <- 1   # Choose an examinee to investigate
#' print(alpha.est.NP.H) # Print the estimated examinee attribute profiles



NPC = function(Y, Q, gate=c("AND", "OR"), method=c("Hamming", "Weighted", "Penalized"), wg=1, ws=1){

  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats

  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  check <- NULL
  check <- CheckInput(Y, Q)
  if (!is.null(check)) return(warning(check))

  gate <- match.arg(gate)
  method <- match.arg(method)

  #####
  # 2 #
  ##### Estimation

  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2 ^ natt

  # all possible alpha vectors

  pattern <- AlphaPermute(natt)

  # ideal responses for each possible alpha vector

  Ideal <- matrix(NA, M, nitem)
  for (m in 1:M) {
    for (j in 1:nitem){
      if (gate == "AND"){
        u <- prod(pattern[m, ] ^ Q[j, ])
      }
      else if (gate == "OR"){
        u <- 1 - prod((1 - pattern[m, ]) ^ Q[j, ])
      }
      else {
        return(warning("Gate specification not valid."))
      }
      Ideal[m, j] <- u
    }
  }

  if (method == "Hamming"){
    if (ws !=1 | wg != 1){
      stop("Guessing weight and slipping weight should be 1.")
    }
    else{weight <- rep(1, nitem)}
  }

  else if (method == "Weighted"){
    if (ws!=1 ||wg!=1){
      stop("Guessing weight and slipping weight should be 1.")
    }
    else{
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))
    weight[weight > 1 / (0.95 * 0.05)] <- 1 / (0.95 * 0.05)
    }
  }

  else if (method == "Penalized"){
    p.bar <- apply(Y, 2, mean)
    weight <- 1 / (p.bar * (1 - p.bar))
    weight[weight > 1 / (0.95 * 0.05)] <- 1 / (0.95 * 0.05)
    if (ws == wg){
      warning("Penalizing weights for guess and slip are the same --> equivalent with the \"Hamming\" method.")
    }
  }
  else {
    return(warning("Method specification not valid."))
  }

  loss.matrix <- matrix(NA, nrow=M, ncol=nperson)
  est.class <- NULL
  est.pattern <- NULL
  n.tie <- rep(0, nperson)

  for (i in 1:nperson){
    Y.matrix <- matrix(rep(Y[i,], M), M, nitem, byrow=TRUE)
    loss <- apply(matrix(rep(weight, M), M, nitem, byrow=TRUE) * (wg * abs(Y.matrix - Ideal) * Y.matrix + ws * abs(Y.matrix - Ideal) * (1 - Y.matrix)), 1, sum)
    loss.matrix[,i] <- loss

    min.loss <- which(loss == min(loss))

    if (length(min.loss) != 1) {
      n.tie[i] <- length(min.loss)
      min.loss <- sample(min.loss, 1, prob=rep(1 / length(min.loss), length(min.loss)))
    }

    est.class <- c(est.class, min.loss)
  }

  est.pattern <- pattern[est.class,]
  est.ideal <- Ideal[est.class,]
  output <- list(alpha.est=est.pattern, est.ideal=est.ideal, est.class=est.class, n.tie=n.tie, pattern=pattern, loss.matrix=loss.matrix, method=method, Q=Q, Y=Y)
  class(output) <- "NPC"
  return(output)
}



