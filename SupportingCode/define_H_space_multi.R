#' Define the hypothesis space
#' @return A list containing the H_space and H_annot (tells whether each row
#' comes from the alternative (1) or null (0)).
#' @export
#'
#' @examples
#' define_H_space(3)
define_H_space <- function(K) {
    H_space <- expand.grid( rep(list(c(-1,0,1)), K) )
    s <- rep(0, nrow(H_space))
    for (col_it in 1:ncol(H_space)) {
        s <- s + abs(H_space[, col_it])
    }
    H_space <- H_space %>% mutate(s = s) %>%
        arrange(s) %>%
        select(-s) %>%
        as.matrix(.)

    H_annot <- rep(0, nrow(H_space))
    H_annot[which(apply(abs(H_space), 1, sum) == K)] <- 1

    return( list(H_space=H_space, H_annot=H_annot) )
}

