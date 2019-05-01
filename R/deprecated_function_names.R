#' reverse-compatibility wrapper around de_selected_pairs()
#' @seealso \code{link{de_selected_pairs}}
DE_genes_pairs <- function(...) {
  de_selected_pairs(...)
}

#' reverse-compatibility wrapper around de_all_pairs()
#' @seealso \code{link{de_all_pairs}}
DE_genes_pw <- function(...) {
  de_all_pairs(...)
}

#' reverse-compatibility wrapper around compute_pair_deScore()
#' @seealso \code{link{compute_pair_deScore}}
de_pair <- function(...) {
  compute_pair_deScore(...)
}

#' reverse-compatibility wrapper around find_vg()
#' @seealso \code{link{find_vg}}
findVG <- function(...) {
  find_vg(...)
}
