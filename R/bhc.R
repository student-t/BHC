#
bhc <- function(data, itemLabels=NULL, nFeatureValues=0, timePoints=NULL, dataType="multinomial", noise=NULL, numReps=0, noiseMode=0, robust=0, numThreads=1, randomised=FALSE, m=2, verbose=FALSE){
  ##----------------------------------------------------------------------
  ## OPTIONALLY, PRINT USEFUL INFORMATION TO SCREEN ----------------------
  ##----------------------------------------------------------------------
  if (verbose){
    print("Running Bayesian Hierarchical Clustering....", quote=FALSE)
    print(paste("DataType:", dataType))
  }
  ##-------------------------------------------------------------------------
  ## CONVERT THE dataType TO AN INTEGER ID (EASIER TO PASS TO THE C++ CODE) - 
  ##-------------------------------------------------------------------------
  if (dataType=="multinomial")
    dataTypeID = 0
  else
    {if (dataType=="time-course")
       dataTypeID = 1
       else
       	{if (dataType=="cubicspline")
       dataTypeID = 2
    	else
          stop(paste("Error!  This isn't a valid dataType:", dataType))
      	}
   }
  ##----------------------------------------------------------------------
  ## FORMAT DATA FOR PASSING TO C++; FIND USEFUL VALUES ------------------
  ##----------------------------------------------------------------------
  data       <- as.matrix(data) ##enforce that it's a data.matrix
  data       <- aperm(data, c(2,1)) ##transpose so the C++ code understands it
  nFeatures  <- nrow(data)
  nDataItems <- ncol(data)
  nLabels    <- length(itemLabels)
  if (nLabels!=nDataItems)
    stop(paste("Error!  There are", nLabels, "labels, but", nDataItems, "data items."))
  if(m<2) m=2
  ##-----------------------------------------------------------------------
  ## MULTINOMIAL:  NORMALISE DATA; FIND NUMBER OF DISCRETE FEATURE VALUES -
  ##-----------------------------------------------------------------------
  if (dataType=="multinomial"){
    data           <- data - min(data) ##the C++ code wants values starting at zero 
    nFeatureValues <- length(unique(data[1:length(data)]))
  }
  ##----------------------------------------------------------------------
  ## IF REQUIRED, OPTIMISE THE GLOBAL HYPERPARAMETER ---------------------
  ##----------------------------------------------------------------------
  ##ultimately, we might want to move this inside the C++ code
  if (dataType=="multinomial")
    globalHyperParam <- FindOptimalHyperparameter(dataTypeID, data, timePoints, noise,
                                                  nDataItems, nFeatures, nFeatureValues,
                                                  verbose=verbose)
  else
    globalHyperParam <- 0 ##dummy value as we don't need a global Hyperparam in this case
  ##----------------------------------------------------------------------
  ## RUN BHC ANALYSIS; CONSTRUCT OUTPUT DENDROGRAM OBJECT ----------------
  ##----------------------------------------------------------------------
  out              <- RunBhcWrapper(globalHyperParam, dataTypeID, data, timePoints,
                                    nDataItems, nFeatures, nFeatureValues, noise,
                                    numReps, noiseMode, robust, fullOutputSwitch=TRUE,
                                    numThreads, randomised, m, verbose)
  outputDendrogram <- ConstructDendrogramObject(out, nDataItems, nFeatures, itemLabels)
  ##----------------------------------------------------------------------
  ## PRINT THE LogEvidence TO SCREEN -------------------------------------
  ##----------------------------------------------------------------------
  if (verbose){
    print(paste("Hyperparameter:", globalHyperParam), quote=FALSE)
    print(paste("Lower bound on overall LogEvidence:",
                format(out$logEvidence, digits=5, trim=TRUE, scientific=TRUE)),
          quote=FALSE);
    print("*******************", quote=FALSE)
  }
  ##return the output dendrogram
  outputDendrogram
}
##*****************************************************************************
##*****************************************************************************
