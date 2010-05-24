##wrapper function to run the interfaced C++ code, ensuring correct formating etc
RunBhcWrapper <- function(globalHyperParameter, data, nDataItems, nFeatures, nFeatureValues, fullOutputSwitch=FALSE, verbose=FALSE){
  ##generate the output structure
  out <- .C("bhcWrapper",
            as.integer(data),
            as.integer(nDataItems),
            as.integer(nFeatures),
            as.double(globalHyperParameter),
            as.integer(nFeatureValues),
            logEvidence=as.double(0),
            node1=vector(mode='integer',length=nDataItems-1),
            node2=vector(mode='integer', length=nDataItems-1),
            mergeOrder=vector(mode='integer', length=nDataItems-1),
            mergeWeight=vector(mode='numeric', length=nDataItems-1),
            PACKAGE="BHC")
  ##optionally, return the logEvidence so we can use optimisation routines
  if (verbose)          print(c(globalHyperParameter, out$logEvidence), quote=FALSE)
  if (fullOutputSwitch) out else out$logEvidence
}
##*****************************************************************************
##*****************************************************************************
