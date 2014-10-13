##Function to find an optimised binning to discretise input data
##It uses Bayesian model selection to identify symmetric percentiles to divide
##the data into three bins, corresponding (for gene expression) to over/normally/under expressed
##
##This approach is only exploring a limited subset of possible schemes, but experience suggests that
##it works well.
##
##This function expects as input a (nItems*nFeatures) array
##
FindOptimalBinning <- function(inputData, sampleNames, transposeData=FALSE, verbose=FALSE) {
  ##----------------------------------------------------
  ##-- FIND USEFUL VALUES ------------------------------
  ##----------------------------------------------------
  logEvidence.best <- -Inf
  percentiles.best <- c(0.33, 0.33, 0.33)
  ##-----------------------------------------------------
  ##-- 1ST PASS TO SEARCH FOR THE OPTIMISED PERCENTILES -
  ##-----------------------------------------------------
  for (i in 2:7){
    v                    <- 5*i 
    percentiles.current  <- 0.01*c(v, 100-2*v, v)
    discreteData         <- DiscretiseData(inputData, percentiles=percentiles.current, verbose=verbose)
    if (transposeData){
      discreteData <- t(discreteData)
      sampleNames  <- vector("character", nrow(discreteData))
    }
    dendrogram           <- bhc(discreteData, sampleNames, 3, verbose=verbose)
    logEvidence.bhc      <- attr(dendrogram, "logEvidence")
    logEvidence.discrete <- attr(discreteData, "logEvidence")
    logEvidence.current  <- logEvidence.bhc + logEvidence.discrete
    ##IF THIS IS A BETTER SOLUTION, KEEP IT
    if (logEvidence.current>=logEvidence.best){
      if (is.finite(logEvidence.current)){
        logEvidence.best <- logEvidence.current
        percentiles.best <- percentiles.current
      }
    }
  }
  ##-----------------------------------------------------
  ##-- 1ST PASS TO SEARCH FOR THE OPTIMISED PERCENTILES -
  ##-----------------------------------------------------
  rangeArray <- 100*percentiles.best[1] - 5 + 1:9
  for (i in rangeArray){
    percentiles.current  <- 0.01*c(i, 100-2*i, i)
    discreteData         <- DiscretiseData(inputData, percentiles=percentiles.current, verbose=verbose)
    if (transposeData){
      discreteData <- t(discreteData)
      sampleNames  <- vector("character", nrow(discreteData))
    }
    dendrogram           <- bhc(discreteData, sampleNames, 3, verbose=verbose)
    logEvidence.bhc      <- attr(dendrogram, "logEvidence")
    logEvidence.discrete <- attr(discreteData, "logEvidence")
    logEvidence.current  <- logEvidence.bhc + logEvidence.discrete 
    ##IF THIS IS A BETTER SOLUTION, KEEP IT
    if (logEvidence.current>logEvidence.best){
      if (is.finite(logEvidence.current)){
        logEvidence.best <- logEvidence.current
        percentiles.best <- percentiles.current
      }
    }
  }
  ##----------------------------------------------------
  ##-- PRINT CONFIRMING INFO TO SCREEN -----------------
  ##----------------------------------------------------
  cat('', fill=TRUE)
  cat('OPTIMISED DISCRETISATION', fill=TRUE)
  cat('------------------------', fill=TRUE)
  cat('Percentiles:', percentiles.best, fill=TRUE)
  cat('LogEvidence:', logEvidence.best, fill=TRUE)
  ##----------------------------------------------------
  ##-- RETURN THE OPTIMISED PERCENTILES ----------------
  ##----------------------------------------------------
  invisible(percentiles.best)
}
##*********************************************************
##*********************************************************
##----------------------------------------------------
##--  ------------------------------------------------
##----------------------------------------------------
