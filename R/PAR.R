#' @title Pattern-wise agreement rate
#'
#' @description The function is used to compute the pattern-wise agreement rate between two sets of
#' attribute profiles. They need to have the same dimensions.
#'
#' @param x One set of attribute profiles
#' @param y The other set of attribute profiles
#'
#' @return The function returns the pattern-wise agreement rate between two sets of attribute profiles.
#'
#' @export
#'
PAR = function(x, y) {
  if (!all(dim(x) == dim(y))) {
    warning("The dimensions of the attribute profiles are not the same.")
    return(NULL)
  } else {
    out = mean(1 * (rowSums(abs(x - y)) == 0))
    return(out)
  }
}
