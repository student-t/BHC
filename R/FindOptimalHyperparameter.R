##Function to run a series of numerical optimsations on the overall log-Evidence, to find
##the optimal value for the global hyperparameter.
##This function assumes that the log-Evidence is well-behaved in the region of param=1
FindOptimalHyperparameter <- function(data, nDataItems, nFeatures, nFeatureValues, verbose=FALSE){
  if (verbose) print("Optimising global hyperparameter...", quote=FALSE)
  ##----------------------------------------------------------------------
  ## BASIC SCAN OF (1D) HYPERPARAMETER SPACE TO FIND ROUGH REGION WHERE --
  ## THE OPTIMAL VALUE LIVES ---------------------------------------------
  ##----------------------------------------------------------------------
  param.lower <- 1
  param.upper <- 1
  logEv.lower <- RunBhcWrapper(param.lower, data, nDataItems, nFeatures, nFeatureValues)  
  logEv.upper <- logEv.lower
  ##find a sensible lower bound
  repeat{
    param.lower <- 0.5 * param.lower
    logEv.new   <- RunBhcWrapper(param.lower, data, nDataItems, nFeatures, nFeatureValues)  
    if ((logEv.new - logEv.lower) < 0) break else logEv.lower = logEv.new
  }
  ##find a sensible upper bound
  repeat{
    param.upper <- 2 * param.upper
    logEv.new   <- RunBhcWrapper(param.upper, data, nDataItems, nFeatures, nFeatureValues)  
    if ((logEv.new - logEv.upper) < 0) break else logEv.upper = logEv.new
  }
  ##----------------------------------------------------------------------
  ## GIVEN A ROUGH RANGE, RUN A SMARTER OPTIMISATION ALGORITHM -----------
  ##----------------------------------------------------------------------
  ##the tolerance here is in log-Evidence units.  This means that 1 is a fairly sensible value, because
  ##values smaller than that descend into the analysis' uncertainty
  ##(1 -> probability ratio of 'e' -> approx a one-sigma significance (Z-score))
  optimalOutput <- optimise(RunBhcWrapper, interval =c(param.lower, param.upper) , data, nDataItems, nFeatures, nFeatureValues,
                            maximum = TRUE, tol=1., verbose=verbose)
  ##return the optimal hyperparameter value
  optimalOutput$maximum
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------

