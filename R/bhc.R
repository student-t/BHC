## (implemented March '08, Rich Savage)
## (C++ code originally written by Xang Yu; modified by Rich Savage)
##
## Bayesian Hierarchical Clustering (BHC) 
## --------------------------------------
## This is an R implementation of the BHC algorithm developed by Heller & Ghahramani (2005), set up to handle multinomial data.
## This R function links to C++ code originally authored by Yang Xu and modified to the multinomial case by Rich Savage.
## That code was based on original (Matlab) code written by Katherine Heller.
##
##The input data should be a data matrix with dimensions nDataItems * nFeatures.  The dataItems will then be clustered.
##
bhc <- function(data, itemLabels, nFeatureValues, verbose=FALSE){
  if (verbose) print("Running Bayesian Hierarchical Clustering....", quote=FALSE)
  ##-----------------------------------------------------------------------------
  ## IF NECESSARY, LOAD DYNAMICALLY THE BHC C++ CODE; INCLUDE NECESSARY R FILES -
  ##-----------------------------------------------------------------------------
  if(!is.loaded("bhcWrapper")){dyn.load("ScienceProjects/BayesianHierarchicalClustering/R_code/BHC/src/BHC.so")}
  ##----------------------------------------------------------------------
  ## FORMAT DATA FOR PASSING TO C++; FIND USEFUL VALUES ------------------
  ##----------------------------------------------------------------------
  data       <- as.matrix(data)
  data       <- aperm(data, c(2,1)); #Matrix transpose.  Do this to interface correctly with C++
  nDataItems <- (dim(data))[2];
  nFeatures  <- (dim(data))[1];
  ##----------------------------------------------------------------------
  ## IDEALLY, WE WANT SOME PRETTY ROBUST ASSERTIONS HERE TO MINIMISE -----
  ## PROBLEMS WITH THE LINKED CODE ---------------------------------------
  ##----------------------------------------------------------------------
  if (min(data) <  0)              stop("Error!  Some data have values <0")
  if (max(data) >= nFeatureValues) stop("Error!  Some data have values >= nFeatureValues (must be less than, as the C++ code starts counting at zero)")
  ##---------------------------------------------------------------------------------------
  ## OPTIMISE GLOBAL HYPERPARAMETER; RUN BHC ANALYSIS; CONSTRUCT OUTPUT DENDROGRAM OBJECT -
  ##---------------------------------------------------------------------------------------
  globalHyperParam <- FindOptimalHyperparameter(data, nDataItems, nFeatures, nFeatureValues, verbose=verbose)
  out              <- RunBhcWrapper(globalHyperParam, data, nDataItems, nFeatures, nFeatureValues, fullOutputSwitch=TRUE, verbose=verbose)  
  outputDendrogram <- ConstructDendrogramObject(out, nDataItems, nFeatures, itemLabels)
  ##----------------------------------------------------------------------
  ## PRINT THE LogEvidence TO SCREEN -------------------------------------
  ##----------------------------------------------------------------------
  if (verbose){
    print(paste("Hyperparameter:", globalHyperParam), quote=FALSE)
    print(paste("Lower bound on overall LogEvidence:", format(out$logEvidence, digits=5, trim=TRUE, scientific=TRUE)), quote=FALSE);
    print("*******************", quote=FALSE)
  }
  ##return the output dendrogram
  outputDendrogram
}
##*****************************************************************************
##*****************************************************************************
