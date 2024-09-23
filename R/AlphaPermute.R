AlphaPermute <- function(dim) {

  alpha <- matrix(c(0, 1), 2, 1)

  for (i in 1:(dim - 1)) {
    alpha <- rbind(alpha, alpha)
    alpha <- cbind(alpha, c(rep(0, 2 ^ i), rep(1, 2 ^ i)))
  }

  return(alpha)
}
