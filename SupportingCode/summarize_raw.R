# Function to summarize raw results
summarize_raw <- function(fullDat, full=FALSE, cor=FALSE, maxP=FALSE, FDP2=FALSE) {
  # summarize raw output results
  outDF <- c()
  effSizes <- sort(unique(fullDat$minEff1))
  for (eff_it in 1:length(effSizes)) {
    # loop through each effect size
    tempEff <- effSizes[eff_it]
    tempDat <- fullDat %>% filter(minEff1 == tempEff)  %>%
      as.data.frame(.) %>%
      mutate(fdpDACT = ifelse(nRejDACT == 0, 0, fdpDACT)) %>%
      mutate(fdpHDMT = ifelse(nRejHDMT == 0, 0, fdpHDMT)) %>%
      mutate(fdpKernel = ifelse(nRejKernel == 0, 0, fdpKernel)) %>%
      mutate(fdp7df = ifelse(nRej7df == 0, 0, fdp7df)) %>%
      mutate(fdp50df = ifelse(nRej50df == 0, 0, fdp50df)) %>%
      mutate(fdpNew = ifelse(nRejNew == 0, 0, fdpNew))

    if (cor) {
      tempDat <- tempDat %>% mutate(fdpCor = ifelse(nRejCor == 0, 0, fdpCor)) %>%
        mutate(fdpNew = fdpCor) %>%
        mutate(powerNew = powerCor) %>%
        mutate(nRejNew = nRejCor)
    }

    if (full) {
      tempDat <- tempDat %>% mutate(fdpFull = ifelse(nRejFull == 0, 0, fdpFull)) %>%
        mutate(fdpDACT = fdpFull) %>%
        mutate(powerDACT = powerFull) %>%
        mutate(nRejDACT = nRejFull)
    }
    
    # summarize
    summaryOut <- data.frame(minEff1 = tempDat$minEff1[1],
                           Method = c("DACT", "HDMT", "Kernel", "df7", "df50", "New"))
    summaryOut$nRej <- c(mean(tempDat$nRejDACT, na.rm=T), mean(tempDat$nRejHDMT, na.rm=T),
                     mean(tempDat$nRejKernel, na.rm=T), mean(tempDat$nRej7df, na.rm=T),
                     mean(tempDat$nRej50df, na.rm=T), mean(tempDat$nRejNew, na.rm=T))
    summaryOut$Power <- c(mean(tempDat$powerDACT, na.rm=T),
                        mean(tempDat$powerHDMT, na.rm=T), mean(tempDat$powerKernel, na.rm=T),
                      mean(tempDat$power7df, na.rm=T), mean(tempDat$power50df, na.rm=T),
                      mean(tempDat$powerNew, na.rm=T))
    summaryOut$FDP <- c(mean(tempDat$fdpDACT, na.rm=T),
                      mean(tempDat$fdpHDMT, na.rm=T), mean(tempDat$fdpKernel, na.rm=T),
                    mean(tempDat$fdp7df, na.rm=T), mean(tempDat$fdp50df, na.rm=T),
                    mean(tempDat$fdpNew, na.rm=T))
    summaryOut$Incongruous <- c(NA, NA, mean(tempDat$inconKernel, na.rm=T), mean(tempDat$incon7df, na.rm=T),
                            mean(tempDat$incon50df, na.rm=T), mean(tempDat$inconNew, na.rm=T))
    summaryOut$numNA <- c(length(which(is.na(tempDat$powerDACT))),
                      length(which(is.na(tempDat$powerHDMT))), length(which(is.na(tempDat$powerKernel))),
                      length(which(is.na(tempDat$power7df))), length(which(is.na(tempDat$power50df))), length(which(is.na(tempDat$powerNew))))
    
    if (maxP) {
      summaryOut <- data.frame(minEff1 = tempDat$minEff1[1],
                           Method = c("MaxP", "DACT", "HDMT", "Kernel", "df7", "df50", "New"))
      summaryOut$nRej <- c(mean(tempDat$nRejMaxp, na.rm=T), mean(tempDat$nRejDACT, na.rm=T), mean(tempDat$nRejHDMT, na.rm=T),
                     mean(tempDat$nRejKernel, na.rm=T), mean(tempDat$nRej7df, na.rm=T),
                     mean(tempDat$nRej50df, na.rm=T), mean(tempDat$nRejNew, na.rm=T))
      summaryOut$Power <- c(mean(tempDat$powerMaxp, na.rm=T), mean(tempDat$powerDACT, na.rm=T),
                        mean(tempDat$powerHDMT, na.rm=T), mean(tempDat$powerKernel, na.rm=T),
                      mean(tempDat$power7df, na.rm=T), mean(tempDat$power50df, na.rm=T),
                      mean(tempDat$powerNew, na.rm=T))
      summaryOut$FDP <- c(mean(tempDat$fdpMaxp, na.rm=T), mean(tempDat$fdpDACT, na.rm=T),
                      mean(tempDat$fdpHDMT, na.rm=T), mean(tempDat$fdpKernel, na.rm=T),
                    mean(tempDat$fdp7df, na.rm=T), mean(tempDat$fdp50df, na.rm=T),
                    mean(tempDat$fdpNew, na.rm=T))
    }

    if (FDP2) {
      summaryOut$sig2Pow <- c(mean(tempDat$sig2PowDACT, na.rm=T),
                        mean(tempDat$sig2PowHDMT, na.rm=T), mean(tempDat$sig2PowKernel, na.rm=T),
                      mean(tempDat$sig2Pow7df, na.rm=T), mean(tempDat$sig2Pow50df, na.rm=T),
                      mean(tempDat$sig2PowNew, na.rm=T))
      summaryOut$sig2FDP <- c(mean(tempDat$sig2fdpDACT, na.rm=T),
                      mean(tempDat$sig2fdpHDMT, na.rm=T), mean(tempDat$sig2fdpKernel, na.rm=T),
                    mean(tempDat$sig2fdp7df, na.rm=T), mean(tempDat$sig2fdp50df, na.rm=T),
                    mean(tempDat$sig2fdpNew, na.rm=T))
    }
    outDF <- rbind(outDF, summaryOut)
  }

  return(outDF)
}

