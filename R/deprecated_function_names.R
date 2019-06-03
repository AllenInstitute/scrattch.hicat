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

#' reverse-compatibility wrapper around score_pair_stats()
#' @seealso \code{link{score_pair_stats}}
de_pair <- function(...) {
  de_stats_pair(...)
}

#' reverse-compatibility wrapper around de_selected_pairs_stats()
#' @seealso \code{link{de_selected_pairs_stats}}
de_score_pairs <- function(...) {
  de_stats_selected_pairs(...)
}

#' reverse-compatibility wrapper around de_all_pairs_stats()
#' @seealso \code{link{de_all_pairs_stats}}
de_score <- function(...) {
  de_stats_all_pairs(...)
}

#' reverse-compatibility wrapper around find_vg()
#' @seealso \code{link{find_vg}}
findVG <- function(...) {
  find_vg(...)
}
