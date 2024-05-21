#' Get pairwise distance
#'
#' Funtion to get pairwise distance of two codesets
#'
#' @param codesA data.frame
#' @param codesB data.frame
#' @return Pairwise distance between code A and B
#' @export
#' @examples
#' # Load necessary libraries and data
#' library(loxcodeR)
#' lox <- readRDS("~/Desktop/loxcodeR/LoxcodeR_app/Week2.rds")
#'
#' # Assuming you have two codesets
#' #codeset_A <- data[data$codeset == "A", ]
#' #codeset_B <- data[data$codeset == "B", ]
#'
#' # Get pairwise distance
#' #pairwise_distance <- get_pairwise_distance(codeset_A, codeset_B)
get_pair_dist <- function(codesA, codesB) {
    # if (!all(c(
    #     all(codesA$is_valid),
    #     all(codesB$is_valid),
    #     all(codesA$size %in% c(13, 9), all(codesB$size %in% c(13, 9)))
    # ))) {
    #     stop("Error: either some of your codes aren't valid, or not size 9, 13")
    # }
    vec_A <- get_cass_vec(codesA$code)
    vec_B <- get_cass_vec(codesB$code)
    return(retrieve_dist_pair(vec_A, vec_B))
}
