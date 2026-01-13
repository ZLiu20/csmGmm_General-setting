#' Calculate the conditional pmf of z given h
#'
#' @param summary_tab A J*K table of the J summary statistics from K studies.
#' @param H_space An num_bins*(2*K) table, two columns from each study, first column
#' is bin midpoint, second is marginal probability of z Pr(Z=z) at that point
#' @param pi0_estim Estimated marginal probability of null, Pr(H=0), for each study.
#' @param checkpoint TRUE to print diagnostic information, FALSE to silence.
#'
#' @return A list with binned_tab, null_pmf_tab, neg_pmf_tab, pos_pmf_tab
#'
#' @export
#'
#' @examples
#' summary_tab <- cbind(rnorm(1000), rnorm(1000))
#' locfdr_output <- locfdr_estim(summary_tab)
#' calc_cond_pmfs(summary_tab, locfdr_output$marg_pmf_tab, locfdr_output$pi0_estim)
#'
calc_cond_pmfs <- function(summary_tab, marg_pmf_tab, pi0_estim, checkpoint=TRUE) {
    if (checkpoint) {
        cat('Calculating null and alternative pmf from marginal... \n')
    }

    # We need to bin the summary tab
    binned_tab <- matrix(data=NA, nrow=nrow(summary_tab), ncol=ncol(summary_tab))
    null_pmf_tab <- matrix(data=NA, nrow=nrow(marg_pmf_tab), ncol=ncol(marg_pmf_tab)/2)
    neg_pmf_tab <- matrix(data=NA, nrow=nrow(marg_pmf_tab), ncol=ncol(marg_pmf_tab)/2)
    pos_pmf_tab <- matrix(data=NA, nrow=nrow(marg_pmf_tab), ncol=ncol(marg_pmf_tab)/2)
    num_bins <- nrow(marg_pmf_tab)
    for (study_it in 1:ncol(summary_tab)) {
        midpts <- marg_pmf_tab[, (2*study_it-1)]
        bin_start <- min(midpts)
        bin_end <- max(midpts)
        # some of these variable names are a bit misleading - for locfdr the default is 120 bins,
        # which returns 119 "midpoints" and between these 119 there are 118 bins, plus two more outside
        # the extreme values. this code is correct in dividing by 118, the wording is just a little confusing. 
        bin_interval <- (bin_end - bin_start) / (num_bins - 1)
        bin_start <- bin_start - bin_interval / 2

        # Get the bins
        # Note that the way I do it means there are actually only 119 bins, I kind of discard the most negative one.
        # But that also makes sense if you think about the breaks as midpoints, there are only 119 midpoints, and you
        # can only assign something to each midpoint. 
        binned_tab[, study_it] <- ceiling((summary_tab[, study_it] - bin_start) / bin_interval)

        # Sometimes one of the endpoints will fall out of the intervals
        too_small <- which(binned_tab[, study_it] < 1)
        too_large <- which(binned_tab[, study_it] > num_bins)
        if (length(too_small) > 0) {binned_tab[too_small, study_it] <- 1}
        if (length(too_large) > 0) {binned_tab[too_large, study_it] <- num_bins}

        # Get the null and alternative pmf from marginal and pi0
        temp_marginal <- marg_pmf_tab[, 2*study_it] / sum(marg_pmf_tab[, 2*study_it])
        null_pmf_tab[, study_it] <- dnorm(midpts)

        # Standardize so it sums to 1
        null_pmf_tab[, study_it] <- null_pmf_tab[, study_it] / sum(null_pmf_tab[, study_it])

        # Then alternative and standardize
        neg_pmf_tab[, study_it] <- (temp_marginal - pi0_estim[study_it]*null_pmf_tab[, study_it]) / (1-pi0_estim[study_it])
        neg_pmf_tab[which(neg_pmf_tab[, study_it] < 0), study_it] <- 0
        neg_pmf_tab[which(midpts > 0), study_it] <- 0
        neg_pmf_tab[, study_it] <- neg_pmf_tab[, study_it] / sum(neg_pmf_tab[, study_it])

        pos_pmf_tab[, study_it] <- (temp_marginal - pi0_estim[study_it]*null_pmf_tab[, study_it]) / (1-pi0_estim[study_it])
        pos_pmf_tab[which(pos_pmf_tab[, study_it] < 0), study_it] <- 0
        pos_pmf_tab[which(midpts < 0), study_it] <- 0
        pos_pmf_tab[, study_it] <- pos_pmf_tab[, study_it] / sum(pos_pmf_tab[, study_it])
    }

    return(list(binned_tab=binned_tab, null_pmf_tab=null_pmf_tab, neg_pmf_tab=neg_pmf_tab, pos_pmf_tab=pos_pmf_tab))
}
