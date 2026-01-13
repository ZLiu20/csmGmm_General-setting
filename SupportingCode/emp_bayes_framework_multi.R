#' Main function to run empirical bayes FDR inference procedure.
#'
#' @param summary_tab An M*N table of summary statistics with M SNPs from N studies.
#' @param q_threshold The q-value we want to control.
#' @param pi0_epsilon EM converges when L1 norm of pi0 difference is less than this value.
#' @param Hdist_epsilon 'Inner EM' stops when ll change less than this value. Then recalculate
#' pi0 and f_hat.
#' @param q_threshold The q-value we want to control.
#' @param checkpoint TRUE to print diagnostic information, FALSE to silence.
#'
#' @return A list with acc_tab (all accepted z vectors), rej_tab (all rejected z vectors),
#' num_rej, pi0_estim, and H_dist.
#'
#' @export
#'
#' @examples
#' #' H_space <- define_H_mediation()$H_space
#' H_annot <- define_H_mediation()$H_annot
#' null_rows <- which(H_annot == 0)
#' summary_tab <- cbind(rnorm(100), rnorm(100)
#' locfdr_output <- locfdr_estim(summary_tab=summary_tab)
#' pi0_estim <- locfdr_output$pi0_vec
#' marg_pmf_tab <- locfdr_output$marg_pmf_tab
#' H_dist <- calc_H_dist_indep(pi0_estim=pi0_estim, H_space=H_space)
#' cond_pmfs <- calc_cond_pmfs(summary_tab=summary_tab, marg_pmf_tab=marg_pmf_tab, pi0_estim=pi0_estim)
#' null_pmf_tab <- cond_pmfs$null_pmf_tab
#' neg_pmf_tab <- cond_pmfs$neg_pmf_tab
#' pos_pmf_tab <- cond_pmfs$pos_pmf_tab
#' binned_tab <- cond_pmfs$binned_tab
#' calc_bayes_FDR_locfdr(binned_tab=binned_tab,
#' H_space=H_space, null_rows=null_rows, null_pmf_tab=null_pmf_tab,
#' neg_pmf_tab=neg_pmf_tab, pos_pmf_tab=pos_pmf_tab, checkpoint=TRUE)
#'
emp_bayes_framework <- function(summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE, ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
    # Step 1 
    H_space_output <- define_H_space(K = ncol(summary_tab))
    H_space <- as.matrix(H_space_output$H_space)
    H_annot <- H_space_output$H_annot
    if (sameDirAlt) {
      null_rows <- which(apply(H_space, 1, sum) != ncol(summary_tab) & apply(H_space, 1, sum) != -ncol(summary_tab))
    } else { 
      null_rows <- which(H_annot == 0)
    }

    ###################################
    # Step 2
    # Find the distribution f and pi(H)
    locfdr_output <- locfdr_estim(summary_tab=summary_tab, kernel=kernel, joint=joint, ind=ind, dfFit = dfFit, nulltype=nulltype)
    
    # can be negative - see first real data analysis
    pi0_estim <- locfdr_output$pi0_estim
    if (length(which(pi0_estim >= 1 | pi0_estim <= 0)) > 0) {
     pi0_estim[which(pi0_estim >= 1 | pi0_estim <= 0)] <- 0.99
    }
    marg_pmf_tab <- locfdr_output$marg_pmf_tab
    # was an error with locfdr
    if (is.null(marg_pmf_tab)) {return(-1)} 
    if (ncol(marg_pmf_tab) < 2*ncol(summary_tab)) {
        return(-1)
    }

    # calculate conditional distributions
    H_dist <- calc_H_dist_indep(pi0_estim=pi0_estim, H_space=H_space)
    cond_pmfs <- calc_cond_pmfs(summary_tab=summary_tab, marg_pmf_tab=marg_pmf_tab, pi0_estim=pi0_estim)
    null_pmf_tab <- cond_pmfs$null_pmf_tab
    neg_pmf_tab <- cond_pmfs$neg_pmf_tab
    pos_pmf_tab <- cond_pmfs$pos_pmf_tab

    # check
    if (is.na(sum(neg_pmf_tab)) | is.na(sum(neg_pmf_tab))) {return(-1)}

    # run EM
    binned_tab <- cond_pmfs$binned_tab
    EMoutput <- run_EM_pmf(binned_dat=binned_tab, H_dist=H_dist, H_space=H_space,
                           null_pmf_tab=null_pmf_tab, neg_pmf_tab=neg_pmf_tab,
                           pos_pmf_tab=pos_pmf_tab, pi0_estim=pi0_estim, epsilon=Hdist_epsilon)

    ###########################
    # Step 3: Calculate the local Bayes FDR
    loc_fdr_num <- apply(EMoutput$jointMat[, null_rows], 1, sum)
    lfdrVec <- loc_fdr_num / EMoutput$probZ

    return(list(marg_pmf_tab = marg_pmf_tab, null_pmf_tab = null_pmf_tab, pi0_estim = pi0_estim, Hdist_final = EMoutput$H_dist,
           neg_pmf_tab = neg_pmf_tab, pos_pmf_tab = pos_pmf_tab, lfdrVec = lfdrVec))
}
