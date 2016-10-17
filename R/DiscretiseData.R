##Function to discretise a set of input gene expression data, for subsequent use
##with the BHC clustering algorithm.
##The data are discretised on a gene-by-gene basis, assigning quantiles of the data into each bin
##(ie. we bin on the basis of number of data, not bin width)
##
##The inputs should be as follows:
##
## inputData   - a (nGenes * nExperiments) array or data frame containing log gene expression values
## percentiles - the percentiles into which to divide the data (one per bin) 
##
DiscretiseData <- function(inputData, percentiles=c(0.15, 0.7, 0.15), verbose=TRUE) {
  ##----------------------------------------------------
  ##-- ERROR-CHECKING ----------------------------------
  ##----------------------------------------------------
  stopifnot(sum(percentiles)==1)  
  ##----------------------------------------------------
  ##-- FIND USEFUL VALUES ------------------------------
  ##----------------------------------------------------
  nGenes                <- nrow(inputData)
  nExperiments          <- ncol(inputData)
  discreteData          <- as.matrix(inputData)
  logEvidence           <- 0
  nBins                 <- length(percentiles)
  cumulativePercentiles <- cumsum(percentiles)
  cumulativePercentiles <- c(0, cumulativePercentiles)##include lower edge
  ##convert these percentiles into indices of the sorted data
  binCounts             <- round(nExperiments * cumulativePercentiles)
  binCounts             <- binCounts[2:(nBins+1)] - binCounts[1:nBins]
  edgeIndex             <- 1 + cumsum(binCounts)
  edgeIndex             <- c(1, edgeIndex) 
  edgeIndex[nBins+1]    <- nExperiments##this must be the final element
  ##----------------------------------------------------
  ##-- PRINT CONFIRMING INFO TO SCREEN -----------------
  ##----------------------------------------------------
  if (verbose){
    cat('', fill=TRUE)
    cat('DATA DISCRETISATION', fill=TRUE)
    cat('-------------------', fill=TRUE)
    cat('Percentiles:', percentiles, '\n')
    cat('We have the following parameters for the data array:', fill=TRUE)
    cat(paste('nGenes:      ', nGenes), fill=TRUE)
    cat(paste('nExperiments:', nExperiments), fill=TRUE)
    cat('***Please check that these are the right way round! (it affects the discretisation)***', fill=TRUE)
  }
  ##----------------------------------------------------
  ##-- DISCRETISE EACH GENE IN TURN --------------------
  ##----------------------------------------------------
  for (i in 1:nGenes){
    ##EXTRACT THE DATA FOR THIS GENE; FIND ITS SORTED ORDER
    currentData <- discreteData[i,]
    sortedData  <- sort(currentData)
    newData     <- vector('numeric', nExperiments)
    dataCounter <- vector('numeric', nBins)
    ##FIND THE REQUIRED BIN EDGES, WIDTHS FOR THESE DATA (based on the percentiles)
    binEdges    <- sortedData[edgeIndex]
    binWidths   <- binEdges[2:(nBins+1)] - binEdges[1:nBins]
    ##SET A MINIMUM BIN WIDTH, TO AVOID INFINITIES
    ##this is a bit hacky, but should be okay
    minWidth <- max(currentData) - min(currentData)
    minWidth <- minWidth / length(currentData)
    for (j in 1:nBins)
      binWidths[j] <- max(binWidths[j], minWidth)
    ##HENCE DISCRETISE THE CURRENT DATA
    for (j in 1:nBins){
      index          <- which(currentData >= binEdges[j])
      newData[index] <- j-1
    }
    ##STORE THE DISCRETE DATA
    discreteData[i,] <- newData
    ##COUNT THE NUMBER OF DATA IN EACH BIN
    for (j in 1:nBins) dataCounter[j] = sum(newData==(j-1))
    ##ADD THE LOG-EVIDENCE CONTRIBUTION FOR THIS GENE
    newLogEv    <- - sum(binCounts * log(binWidths))
    logEvidence <- logEvidence + newLogEv
    ##IF THERE'S AN INFINITY, PRINT OUT DIAGNOSTIC INFO
    if (is.finite(newLogEv)==0){
      cat("----------------", fill=TRUE)
      cat("(infinite logEvidence found)", fill=TRUE)
      cat("gene number: ", i, fill=TRUE)
      cat("newLogEv:    ", newLogEv, fill=TRUE)
      cat("binWidths:   ", as.numeric(binWidths), fill=TRUE)
      cat("binEdges:    ", as.numeric(binEdges), fill=TRUE)
      cat("currentData: ", as.numeric(currentData), fill=TRUE)
      cat("nZeros:      ", sum(currentData==0), fill=TRUE)
    }
  }
  ##----------------------------------------------------
  ##-- PRINT THE LogEvidence TO SCREEN -----------------
  ##----------------------------------------------------
  if (verbose){
    cat('', fill=TRUE)
    cat(paste('Discretisation logEvidence:', logEvidence), fill=TRUE)
    cat('(Need to add this to the model logEvidence)', fill=TRUE)
    cat('-------------------', fill=TRUE)
  }
  ##----------------------------------------------------
  ##-- RETURN THE DISCRETISED DATA ---------------------
  ##----------------------------------------------------
  attr(discreteData, "logEvidence") <- logEvidence
  discreteData
}
##*********************************************************
##*********************************************************
##----------------------------------------------------
##--  ------------------------------------------------
##----------------------------------------------------
