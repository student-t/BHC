##Function to handle construction of a "dendrogram" object.
##It should just be passed the output from the C++ BHC code.
##Once the output information is placed in a dendrogram object, various pre-existing object methods in R can be used to generate
##output plots and manipulate the dendrogram.
ConstructDendrogramObject <- function(out, nDataItems, nFeatures, itemLabels){
    ##there's code that does this sort of thing in the as.dendrogram.r source code (on the web)
  leftChild   <- out$node1;
  rightChild  <- out$node2;
  ##find the labels for the nodes (logEv associated with that merge)
  ##also find the labels for the items being clustered
  labelArray <- as.character(format(out$mergeWeight, digits=3, trim=TRUE, scientific=FALSE))
  ##extract the merge weights (log Evidence ratios)
  mergeWeight <- out$mergeWeight
  ##find some sensible heights for the dendrograms
  mergeHeight <- max(mergeWeight) - mergeWeight
  mergeHeight <- sqrt(mergeHeight)
  ##initialise a counter for the number of mergers
  nTotal                 <- 2*nDataItems - 1
  nMembers               <- rep(0, times=nTotal)
  nMerges                <- rep(0, times=nTotal)
  nMembers[1:nDataItems] <- 1L
  ##initialise the tree list (stores all the intermediate mergers); define the leaf nodes
  tree   <- rep(list(0), nTotal)
  for (i in 1:nDataItems) {#leaf nodes
    leaf                      <- i
    attr(leaf, 'logEvidence') <- Inf
    attr(leaf, "members")     <- 1L
    attr(leaf, "height")      <- 0.
    attr(leaf, "leaf")        <- TRUE
    attr(leaf, "midpoint")    <- 0.
    attr(leaf, "label")       <- itemLabels[i]
    attr(leaf, "edgePar")     <- list(lab.col = "blue", lab.cex=.7, col="blue", pch=16)
    tree[[i]]                 <- leaf  
  }
  ##loop over the mergers
  for (i in 1:(nDataItems-1)) {
    newIndex   <- i + nDataItems
    leftMerge  <- as.integer(leftChild[i])
    rightMerge <- as.integer(rightChild[i])
    ##form the new merged list
    leftList  <- tree[[leftMerge]]
    rightList <- tree[[rightMerge]]
    newMerge  <- list(leftList, rightList)
    ##enforce a minimum height, so the dendrograms look sensible
    if (mergeHeight[i] < attr(leftList, "height")) mergeHeight[i] <- 1.02 * attr(leftList, "height")
    if (mergeHeight[i] < attr(rightList,"height")) mergeHeight[i] <- 1.02 * attr(rightList,"height")
    ##standardise the merge height of pairs of leaves (this makes the plottd dendrogram look nicer)
    if (attr(leftList, "leaf") & attr(rightList, "leaf")) mergeHeight[i] <- 0.05*max(mergeHeight)
    ##define some required attributes for the new list
    nMembers[newIndex]            <- nMembers[leftMerge] + nMembers[rightMerge] 
    nMerges[newIndex]             <- 1 + max(nMerges[leftMerge], nMerges[rightMerge])
    attr(newMerge, "members")     <- nMembers[newIndex]
    attr(newMerge, "height")      <- mergeHeight[i]
    attr(newMerge, "logEvidence") <- mergeWeight[i]
    attr(newMerge, "leaf")        <- FALSE
    ##decide edge colour based on whether merger is accepted
    if (mergeWeight[i] < 0){
      ##mark nodes that lead to valid mergers
      if (attr(leftList, 'logEvidence')>=0 | attr(rightList, 'logEvidence')>=0)
        attr(newMerge, "nodePar")  <- list(lab.col = "red", lab.cex=.7, col="red", pch=16)
      ##label the edges
      attr(newMerge, "edgePar")  <- list(lab.col = "red", lab.cex=.7, col="red", pch=16, lty="dotted")
      attr(newMerge, "edgetext") <- labelArray[i]
      ##also want to re-label the child edges, so the plot looks sensible
      attr(newMerge[[1]], "edgePar")  <- list(lab.col = "red", lab.cex=.7, col="red", pch=16, lty="dotted")
      attr(newMerge[[2]], "edgePar")  <- list(lab.col = "red", lab.cex=.7, col="red", pch=16, lty="dotted")
    }
    else{
      attr(newMerge, "edgePar")  <- list(lab.col = "blue", lab.cex=.7, col="blue", pch=16)
    }
    ##set the midpoint so that the plot looks okay
    leftMidpoint  <- 0.5*nMembers[leftMerge]
    rightMidpoint <- 0.5*nMembers[rightMerge] + nMembers[leftMerge]
    attr(newMerge, "midpoint") <- mean(c(leftMidpoint,rightMidpoint)) - 0.5    ##this gives a reasonable-looking dendrogram
    ##store the new list in the tree
    tree[[newIndex]] <- newMerge
  }
  ##and the result we want should be in the last element of tree
  outputDendrogram                      <- newMerge
  class(outputDendrogram)               <- "dendrogram"
  attr(outputDendrogram, "logEvidence") <- out$logEvidence
                                                                              attr(outputDendrogram, "time") <- out$t
  ##return the output dendrogram
  outputDendrogram
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------







