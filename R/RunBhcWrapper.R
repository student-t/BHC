##Wrapper function to run the interfaced C++ code, ensuring correct formating etc
RunBhcWrapper <- function(globalHyperParameter, dataTypeID, data, timePoints, nDataItems, nFeatures, nFeatureValues, noise, numReps=0, noiseMode=0, robust=0, fullOutputSwitch=FALSE, numThreads=1, verbose=FALSE){
  ##generate the output structure
  if(dataTypeID==0){ # Multinomial case: use the 1.1.0 code
    out <- .C("bhcWrapper_multinomial",
            as.integer(data),
            as.integer(nDataItems),
            as.integer(nFeatures),
            as.double(globalHyperParameter),
            as.integer(nFeatureValues),
            logEvidence=as.double(123),
            node1=vector(mode='integer',length=nDataItems-1),
            node2=vector(mode='integer', length=nDataItems-1),
            mergeOrder=vector(mode='integer', length=nDataItems-1),
            mergeWeight=vector(mode='numeric', length=nDataItems-1),
            PACKAGE="BHC")
  }else
  {
    out <- .C("bhcWrapper",
              as.integer(dataTypeID),
              as.double(data), ##the data(matrix) is input as a vector read down the matrix columns; only the data is input not the row or column names
              as.double(timePoints),
              as.integer(nDataItems),
              as.integer(nFeatures),
              as.double(globalHyperParameter),
              as.double(noise),
              as.integer(numReps),
              as.integer(noiseMode),
              as.integer(robust),
              as.integer(nFeatureValues),
              logEvidence=as.double(123),
              node1=vector(mode='integer',length=nDataItems-1),
              node2=vector(mode='integer', length=nDataItems-1),
              mergeOrder=vector(mode='integer', length=nDataItems-1),
              mergeWeight=vector(mode='numeric', length=nDataItems-1),
              as.integer(numThreads))
  }
  ##PACKAGE="base")
  ##optionally, return the logEvidence so we can use optimisation routines
  if (verbose) print(c(globalHyperParameter, out$logEvidence), quote=FALSE)
  if (fullOutputSwitch) out else out$logEvidence
}
##*****************************************************************************
##*****************************************************************************
