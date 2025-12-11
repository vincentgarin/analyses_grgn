#' nm_0_lead
#'
#' @description function to create vector with leading 0
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

nm_0_lead <- function(n){

  x <- 1:n
  str_pad(x, width = max(nchar(x)), pad = "0")

}
