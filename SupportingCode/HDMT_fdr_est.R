# Largely copied from HDMT package, some speed-ups implemeneted
fdr_est <-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues){

  ## alpha10,alpha01,alpha00 are estimated three types of null proportions
  ## alpha1 is the marginal null proportion for first p-value
  ## alpha2 is the marginal null proportion for second p-value
  ## input pvalues are two columns of p-values
  ## alpha is the level of FWER to be control at
  ## exact=0 corresponding to the approximation used in section 2.2-2.3 in the paper, the default value for exact is 0
  ## exact=1 corresponding to the exact used in section 2.4 in the paper
  ## check input

  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")

  # new faster calculation - assume exact = 0 always
  pmax <- apply(input_pvalues,1,max)
  orderedP <- data.frame(input_pvalues) %>% mutate(pmax = pmax) %>%
    mutate(origIdx = 1:nrow(.)) %>%
    arrange(pmax) %>% 
    mutate(avgRej = 1:nrow(.) / nrow(.)) %>%
    mutate(fdr = (pmax * (alpha01 + alpha10) + pmax^2 * alpha00) / avgRej) %>%
    #mutate(fdrAvg = cummean(fdr)) %>%
    arrange(desc(pmax))
    
  # force decreasing fdr
  fixedFdr <- orderedP$fdr 
  for (row_it in 2:length(fixedFdr)) {
    fixedFdr[row_it] <- min(fixedFdr[row_it - 1], fixedFdr[row_it])
  } 

  # rearrange
  fixedP <- orderedP %>% mutate(fixedFdr = fixedFdr) %>%  
    arrange(origIdx)

  return(fixedP)
}


