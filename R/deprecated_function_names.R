#' reverse-compatibility wrapper around de_selected_pairs()
#' @param ... parameters passed to \code{link{de_selected_pairs}}
#' 
#' @seealso \code{link{de_selected_pairs}}
#' @export
DE_genes_pairs <- function(...) {
  de_selected_pairs(...)
}

#' reverse-compatibility wrapper around de_all_pairs()
#' @param ... parameters passed to \code{link{de_all_pairs}}
#' 
#' @seealso \code{link{de_all_pairs}}
#' @export
DE_genes_pw <- function(...) {
  de_all_pairs(...)
}

#' reverse-compatibility wrapper around score_pair_stats()
#' @param ... parameters passed to \code{link{score_pair_stats}}
#' 
#' @seealso \code{link{score_pair_stats}}
#' @export
de_pair <- function(...) {
  de_stats_pair(...)
}

#' reverse-compatibility wrapper around de_selected_pairs_stats()
#' @param ... parameters passed to \code{link{de_selected_pairs_stats}}
#' 
#' @seealso \code{link{de_selected_pairs_stats}}
#' @export
de_score_pairs <- function(...) {
  de_stats_selected_pairs(...)
}

#' reverse-compatibility wrapper around de_all_pairs_stats()
#' @param ... parameters passed to \code{link{de_all_pairs_stats}}
#' 
#' @seealso \code{link{de_all_pairs_stats}}
#' @export
de_score <- function(...) {
  de_stats_all_pairs(...)
}

#' reverse-compatibility wrapper around find_vg()
#' @param ... parameters passed to \code{link{find_vg}}
#' 
#' @seealso \code{link{find_vg}}
#' @export
findVG <- function(...) {
  find_vg(...)
}
