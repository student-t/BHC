##Function to write out the items labels for each cluster.
WriteOutClusterLabels <- function(dendro, outputFile="", verbose=FALSE){
  ##----------------------------------------------------------------------
  ## DEFINE SOME FUNCTIONS TO USE RECURSIVELY ON THE DENDROGRAM NODES ----
  ##----------------------------------------------------------------------
  ##for ease, we use discrete height labels here
  ##this hardwires the logEvidence threshold at zero, which is the point
  ##where merge and not-merge hypotheses are equally likely  
  WhereToCut <- function(n){
    attr(n,"height") <- 1##default
    if (!is.leaf(n)){
      attr(n,"height") <- 2   
      if (attr(n, "logEvidence")<0)
        attr(n,"height") <- 3
    }
    n
  }
  ##----------------------------------------------------------------------
  ## PROCESS THE DENDROGRAM NODES RECURSIVELY ----------------------------
  ##----------------------------------------------------------------------
  dendro <- dendrapply(dendro, WhereToCut);
  ##----------------------------------------------------------------------
  ## CUT THE DENDROGRAM AND PRINT THE LABELS IN EACH CLUSTER -------------
  ##----------------------------------------------------------------------
  cutDendro     <- cut(dendro, 2)
  nClusters     <- length(cutDendro$lower)
  nTotalLabels  <- length(labels(dendro))
  outputStrings <- rep("", nTotalLabels+nClusters)
  counter       <- 1
  
  for (i in 1:nClusters) {
    ##extract the current dendrogram
    currentCluster <- cutDendro$lower[[i]]
    currentLabels  <- labels(currentCluster) 
    nLabels        <- length(currentLabels) 
    ##for each cluster, construct and store the labels
    outputStrings[counter] <- paste("---CLUSTER", i, "---")
    counter                <- counter + 1
    for (j in 1:nLabels){
      outputStrings[counter] <- currentLabels[j]
      counter                <- counter + 1
    }
  }
  ##----------------------------------------------------------------------
  ## IF REQUIRED, WRITE OUT THE CLUSTER LABELS TO A FILE -----------------
  ##----------------------------------------------------------------------
  if (outputFile!="") write.table(outputStrings, file=outputFile, quote=FALSE, row.names=FALSE)
  ##----------------------------------------------------------------------
  ## IF REQUIRED, PRINT THE CLUSTER LABELS OUT TO SCREEN -----------------
  ##----------------------------------------------------------------------
  if (verbose) for (i in 1:length(outputStrings)) print(outputStrings[i], quote=FALSE)
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------

