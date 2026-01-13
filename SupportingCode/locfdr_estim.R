#' Estimate marginal distribution and proportion of nulls
#'
#' @param summary_tab An M*N table of summary statistics where M is the number of SNPs
#' and N is the number of studies.
#'
#' @return A list with pi0_vec (N*1 vector estimating proportion of nulls) and marg_pmf_tab
#' (120*N table of marginal pmfs)
#'
#' @export
#'
#' @examples locfdr_out <- locfdr_estim(summary_tab=cbind(rnorm(100), rnorm(100)))
#'

locfdr_estim <- function(summary_tab, kernel=TRUE, joint=TRUE, ind=TRUE, dfFit=7, nulltype=1) {
    marg_pmf_tab <- c()
    pi0_estim <- c()
    
    # find the gridpoints to evaluate the density
    gridEnds <- matrix(data=NA, nrow=2, ncol=ncol(summary_tab)) 
    for (col_it in 1:ncol(summary_tab)) {
        # there will be 119 midpoints, where the extreme min/max values are exactly half a bin width
        # outside the last midpoints. 120 bins total, 119 bin lengths between min and max observed value.
        binWidth <- (max(summary_tab[, col_it]) - min(summary_tab[, col_it])) / 119
        gridEnds[1, col_it] <- min(summary_tab[, col_it]) + binWidth / 2
        gridEnds[2, col_it] <- max(summary_tab[, col_it]) - binWidth / 2
    }
    
    # kernel joint doesn't need a loop
    if (kernel & joint) {
        # using the KS package to do 2D estimation
        # Hpi is unstructured H matrix, Hpi.diag is diagonal H matrix (for uncorrelated data)
        # remember the rule of thumb is 1.06 * sd * n^(-1/5)
        if (ind) {
            HpiEst <- Hpi.diag(x=summary_tab)
        } else {
            HpiEst <- Hpi(x=summary_tab)
        }
        # estimation
        fhatJoint <- kde(x=summary_tab, H=HpiEst, gridsize=119, xmin=gridEnds[1, ], xmax=gridEnds[2, ], binned=FALSE)
        marg_pmf_tab <- fhatJoint$estimate
        evalsPoints <- cbind(fhatJoint$eval.points[[1]], fhatJoint$eval.points[[2]])
    } else {
        # kernel separately or locfdr separately
        marg_pmf_tab <- c()
        evalPoints <- NA
        for (study_it in 1:ncol(summary_tab)) {

            # locfdr for pi0
            tempLocfdr <- suppressWarnings( tryCatch(locfdr::locfdr(zz=summary_tab[, study_it], bre=120, df=dfFit,
                                                                        pct=0, pct0=1/4, nulltype=nulltype, type=0,
                                                                        plot=0, sw=2),  error=function(e) e) )
            # if error, then skip
            if (class(tempLocfdr)[1] %in% c('simpleError')) {
                pi0_estim <- c(pi0_estim, 0.99)
                next
            } else {
                pi0_estim <- c(pi0_estim, tempLocfdr$pds[1])
                if (!kernel) {
                    marg_pmf_tab <- cbind(marg_pmf_tab, tempLocfdr$x, tempLocfdr$f)
                }
            }

            # kernel
            if (kernel) {
                hEst <- ks::hpi(x=summary_tab[, study_it])
                tempKde <- ks::kde(x=summary_tab[, study_it], H=hEst, gridsize=119, xmin=gridEnds[1, study_it], 
                               xmax=gridEnds[2, study_it], binned=FALSE)
                marg_pmf_tab <- cbind(marg_pmf_tab, tempKde$eval.points, tempKde$estimate)
            } 
        } # done looping through columns of summary_tab
    }

    # teturn
    return(list(marg_pmf_tab = marg_pmf_tab, evalPoints = evalPoints, pi0_estim = pi0_estim)) 
}
