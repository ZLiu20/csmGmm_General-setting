# get rid of the n-1 in base R
myvar_fun <- function(x) {
    sum(x^2) / length(x)
}

# function to allocate the signals
allocate_sigs <- function(hMat, nSNPs, nDims) {

  # total number of signal locations
  hSigsDF <- hMat %>% filter(s != 0)
  totSigs <- sum(hSigsDF$number)

  # sample the locations
  tempPos <- sample(x=1:nSNPs, size=totSigs, replace=F)

  # will put the signal locations here
  sigLocsMat <- matrix(data=NA, nrow=length(tempPos), ncol=nDims)

  # loop through all types of signals
  counter <- 1
  for (row_it in 1:nrow(hSigsDF)) {
    tempNum <- hSigsDF$number[row_it]
    if (tempNum == 0) {next}
    for (col_it in 1:nDims) {
      sigLocsMat[counter:(counter + tempNum - 1), col_it] <- hSigsDF[row_it, col_it] * tempPos[counter:(counter + tempNum - 1)]
    }
    # increment counter
    counter <- counter + tempNum
  }

  return(sigLocsMat)
}


# place signals
set_beta <- function(sigLocsMat, set_it, setSize, betaMin, betaMax) {

  # return this 
  betaMat <- matrix(data=0, nrow=setSize, ncol=ncol(sigLocsMat))

  # indices of this set 
  tempIdx <- ((set_it - 1) * setSize + 1):(set_it * setSize)

  # loop through dimensions to add signals 
  for (col_it in 1:ncol(sigLocsMat)) {
    whichPos <- which(tempIdx %in% sigLocsMat[, col_it])
    whichNeg <- which((-1 * tempIdx) %in% sigLocsMat[, col_it])
    if (length(whichPos) > 0) {
      betaMat[whichPos, col_it] <- rep(1, length(whichPos))
    }
    if (length(whichNeg) > 0) {
      betaMat[whichNeg, col_it] <- rep(-1, length(whichNeg))
    }
    
    if (length(whichPos) + length(whichNeg) == 0) {next}

    # beta values
    betaMat[, col_it] <- betaMat[, col_it] * runif(n = setSize, min=betaMin[col_it], max=betaMax[col_it])
  }

  # return
  return(betaMat)
}

