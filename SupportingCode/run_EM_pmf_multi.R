#' Run the EM algorithm to estimate the distribution of H. Use this one for iterations
#' after the first, when we have approximated the densities f_hat as pmfs.
#' Takes in B*M tables where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins).
#'
#' @param binned_dat An M*N table of binned summary statistics from M SNPs and N studies
#' @param H_dist The probability of each element in H.
#' @param H_space Enumerates all possible values of h.
#' @param null_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) under the null.
#' @param neg_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) when the hypothesis comes from the negative alternative.
#' @param pos_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) when the hypothesis comes from the positive alternative.
#' @param pi0_estim Estimated marginal probabilities of null hypothesis for each study
#' @param epsilon EM algorithm stops when difference in loglik falls below epsilon
#' @param checkpoint Print diagnostic information.
#' @return A list with the objects H_dist (distribution of H), init_loglik (initial loglikelihood).
#'  end_loglik (ending log likelihood), and iterations (number of interations to converge)
#' @export
#'
#' @examples
#' binned_dat <- cbind(40:50, 80:70)
#' H_space <- define_H_mediation$H_space
#' H_dist <- H_dist_indep(h_row=c(1,0), pi0_estim=pi0_estim)
#' null_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120))
#' neg_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120), mean=-1)
#' pos_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120), mean=1)
#' run_EM(binned_dat, H_dist, H_space, pi0_estim=c(0.8,0.8),
#' null_pmf_tab, neg_pmf_tab, pos_pmf_tab)
#'
run_EM_pmf <- function(binned_dat, H_dist, H_space, pi0_estim, epsilon=10^(-2), null_pmf_tab, neg_pmf_tab, pos_pmf_tab, checkpoint=TRUE) {

    # Loop until convergence
    diff <- 1
    num_h <- nrow(H_space)
    counter <- 0
    while(diff > epsilon) {

        # this is the E step
        # get the joint Pr(H=h & Z=z)
        jointMat <- sapply(X=as.data.frame(t(H_space)), FUN=calc_fz_given_h_pmf, z_binned = binned_dat,
                           null_pmf_tab = null_pmf_tab, neg_pmf_tab = neg_pmf_tab, pos_pmf_tab = pos_pmf_tab) %>%
            sweep(., MARGIN=2, STATS=H_dist, FUN="*")
        probZ <- apply(jointMat, 1, sum)
        AikMat <- jointMat %>% sweep(., MARGIN=1, STATS=probZ, FUN="/")
        # this is the M step
        new_H_dist <- apply(AikMat, 2, mean)

        # composite total loglikelihood - now use this as convergence criterion
        totLoglik <- sum(log(probZ))
        if (counter == 0) {oldLoglik <- totLoglik - 10}
        diff <- totLoglik - oldLoglik

        # L2 difference in H_dist, was the old convergence criterion, now we don't use it anymore
        diffSecondary <- sum((new_H_dist - H_dist)^2)
        # new pi, not necessary but just to show it on the screen
        newPi <- calc_pi0_from_H_dist(H_dist = new_H_dist, H_space = H_space)

        # Checkpoint
        counter <- counter + 1
        if (checkpoint) {
            cat('Run:', counter, 'New H dist:', new_H_dist, 'Old H dist:', H_dist, '\n')
            cat('Diff Primary:', diff, '\n')
            cat('Diff Secondary:', diffSecondary, '\n')
            cat('New pi:', newPi, '\n')
            cat('New loglik:', totLoglik, 'Old loglik:', oldLoglik, '\n')
        }

        # Update
        H_dist <- new_H_dist
        oldLoglik <- totLoglik
    }

    return( list(H_dist=H_dist, newPi = newPi, counter=counter, jointMat = jointMat, probZ = probZ) )
}



#' Calculate the density of z_j given h_j, i.e. the
#' density of the test statistics observed for a single SNP given the
#' hypothesis vector, assuming the elements of z_j are independent given h_j.
#' @param h Observed value of hypothesis states
#' @param z_binned Vector of binned test statistics for one SNP
#' @param null_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) under the null.
#' @param neg_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) when the hypothesis comes from the negative alternative.
#' @param pos_pmf_tab B*M table where B is the number of bins and M is the number of studies.
#' The Bth row corresponds to the PMF for the Bth bin (assume each study given same
#' number of bins) when the hypothesis comes from the positive alternative.
#'
#' @export
#' @examples
#' null_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120))
#' neg_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120), mean=-1)
#' pos_pmf_tab=dnorm(x=seq(from=-5,to=5,length.out=120), mean=1)
#' calc_fz_given_h_pmf(h=c(-1,1), z_binned=c(0.5, 1.5), null_pmf_tab=null_pmf_tab,
#' neg_pmf_tab=neg_pmf_tab, pos_pmf_tab=pos_pmf_tab)
#'
calc_fz_given_h_pmf <- function(h, z_binned, null_pmf_tab, neg_pmf_tab, pos_pmf_tab) {

    tempProd <- rep(1, nrow(z_binned))
    for (h_it in 1:length(h)) {
        if (h[h_it] != 0) {
            if (h[h_it] == 1) {
                tempProd <- tempProd * pos_pmf_tab[z_binned[, h_it], h_it]
            } else {
                tempProd <- tempProd * neg_pmf_tab[z_binned[, h_it], h_it]
            }
        } else {
            tempProd <- tempProd * null_pmf_tab[z_binned[, h_it], h_it]
        }
    }

    return(tempProd)
}


#' Calculates the probability of h (an element of H) if all the studies are independent,
#' e.g. Pr(H_j=1,H_k=1) = Pr(H_j=1)*Pr(H_k=1). Used to make an initial guess
#' at Pr(H), although obviously not true because, for example, a top SNP
#' will probably replicate in multiple studies.
#'
#' @param pi0_estim Estimated marginal probabilities that each dimension of the hypothesis space is null.
#' @param H_space Each row of this matrix holds a possible value of h. Each column represents
#' a different dimension of the hypothesis space.
#' @export
#'
#' @examples
#' H_space <- define_H_mediation()$H_space
#' H_dist_indep(pi0_estim=c(0.9,0.9), H_space=H_space)
#'
calc_H_dist_indep <- function(pi0_estim, H_space) {

    multiply_pi0 <- function(h, pi0_estim) {
        prod(pi0_estim^(1-h)) * prod(((1-pi0_estim)/2)^h)
    }

    H_dist_indep <- apply(abs(H_space), 1, multiply_pi0, pi0_estim=pi0_estim)

    return(H_dist_indep)
}

#' Calculates Pr(H_j=0) the marginal probability of null by summing all the
#' probabilities for elements of H_space where H_j = 0.
#'
#' @param H_space Matrix enumerates all possible values of h
#' @param H_dist Estimated distribution of H_space, each element is the probability
#' of the corresponding row in H_space.
#' @export
#'
#' @examples
#' H_space <- define_H_mediation()$H_space
#' pi0_estim <- c(0.99, 0.99)
#' H_dist <- apply(H_space, 1, H_dist_indep, pi0_estim=pi0_estim)
#' pi0_from_Hdist(H_space=H_data$H_space, H_dist=H_dist)
calc_pi0_from_H_dist <- function(H_dist, H_space) {

    new_pi0 <- rep(NA, ncol(H_space))
    for(study_it in 1:ncol(H_space)) {
        new_pi0[study_it] <- sum(H_dist[which(H_space[, study_it] == 0)])
    }

    new_pi0
}

