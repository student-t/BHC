\name{BHC-package}
\alias{BHC-package}
\alias{BHC}
\docType{package}
\title{
Bayesian Hierarchical Clustering
}
\description{
The BHC method performs bottom-up hierarchical clustering, using a
Dirichlet Process (infinite mixture) to model uncertainty in the data
and Bayesian model selection to decide at each step which clusters to
merge.  This avoids several limitations of traditional methods, for
example how many clusters there should be and how to choose a
principled distance metric.  This implementation accepts multinomial
(i.e. discrete, with 2+ categories) data.
}
\details{
\tabular{ll}{
Package: \tab BHC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-02\cr
License: \tab GPL-3\cr
}
The Bayesian Hierarchical Clustering (BHC) algorithm is a black-box,
accessed via the bhc() function (see help for details on how to do
this).  All other functions in this package are accessed via bhc() and
the user should not need to access them directly.
}
\author{
Rich Savage (C++ code originally written for binomial case by Yang Xu)

Maintainer: Rich Savage <r.s.savage@warwick.ac.uk>
}
\references{\emph{Bayesian Hierarchical Clustering}, Heller +
  Ghahramani, Gatsby Unit Technical Report GCNU-TR 2005-002 (2005); also
  see shorter version in ICML-2005; \emph{R/BHC: fast Bayesian
  hierarchical clustering for microarray data}, Savage et al, BMC
  Bioinformatics 10:242 (2009)}
\keyword{ package }
\examples{
require(graphics)
require(BHC)
require(affydata)
require(gcrma)

data(Dilution)
ai        <- compute.affinities(cdfName(Dilution))
Dil.expr  <- gcrma(Dilution,affinity.info=ai,type="affinities")
testData  <- exprs(Dil.expr)
keep      <- sd(t(testData))>0
testData  <- testData[keep,]
testData  <- testData[1:100,]
geneNames <- row.names(testData)

nGenes         <- (dim(testData))[1];
nFeatures      <- (dim(testData))[2];
nFeatureValues <- 4
##NORMALISE EACH EXPERIMENT TO ZERO MEAN, UNIT VARIANCE
for (i in 1:nFeatures){
    newData      <- testData[,i]
    newData      <- (newData - mean(newData)) / sd(newData)
    testData[,i] <- newData
}
##DISCRETISE THE DATA ON A GENE-BY-GENE BASIS
##(defining the bins by equal quartiles)
for (i in 1:nGenes){
  newData      <- testData[i,]
  newData      <- rank(newData) - 1
  testData[i,] <- newData
}
##PERFORM THE CLUSTERING
hc <- bhc(testData, geneNames, nFeatureValues=nFeatureValues)
plot(hc, axes=FALSE)
}
