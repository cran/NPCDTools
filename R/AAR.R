#' @title Attribute-wise agreement rate
#'
#' @description The function is used to compute the attribute-wise agreement rate between
#' two sets of attribute profiles. They need to have the same dimensions.
#'
#' @param x One set of attribute profiles
#' @param y The other set of attribute profiles
#'
#' @return The function returns the attribute-wise agreement rate between two sets of attribute profiles.
#'
#' @export
#'
AAR = function(x, y) {
  if (!all(dim(x) == dim(y))) {
    warning("The dimensions of the attribute profiles are not the same.")
    return(NULL)
  } else {
    aar = 1 - mean(abs(x - y))
    return(aar)
  }
}
