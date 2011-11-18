\name{bhc}
\alias{bhc}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{WriteOutClusterLabels}
\alias{DiscretiseData}
\alias{FindOptimalBinning}
\title{Function to perform Bayesian Hierarchical Clustering on a 2D
  array of discretised (i.e. multinomial) data}
\description{
The method performs bottom-up hierarchical clustering, using a Dirichlet
Process (infinite mixture) to model uncertainty in the data and Bayesian
model selection to decide at each step which clusters to merge.  This
avoids several limitations of traditional methods, for example how many
clusters there should be and how to choose a principled distance metric.
This implementation accepts multinomial (i.e. discrete, with 2+
categories) or time-series data. This version also includes a randomised
algorithm which is more efficient for larger data sets.
}
\usage{
bhc(data, itemLabels, nFeatureValues, timePoints, dataType,
    noise, numReps, noiseMode, robust, numThreads, randomised, m, verbose)
}
\details{
  Typical usage for the multinomial case:
  \preformatted{bhc(data, itemLabels).}
  To use the squared-exponential covariance:
  \preformatted{bhc(data, itemLabels, 0, timePoints, "time-course",
    noise, numReps, noiseMode),}
  and the cubic spline covariance:
  \preformatted{bhc(data, itemLabels, 0, timePoints, "cubicspline",
    noise, numReps, noiseMode).}
  To use the randomised algorithm, simply include the following two
  arguments:
  \preformatted{bhc(data, itemLabels, 0, timePoints, "time-course",
    noise, numReps, noiseMode, randomised=TRUE, m=10)}
    }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A 2D array containing discretised data. The dimensions of
    \code{data} should be \code{nDataItems * nFeatures}, and the
    algorithm will cluster the data items.
  }
  \item{itemLabels}{A character array containing \code{nDataItems} entries,
    one for each data item in the analysis. The leaf nodes of the
    output dendrogram will be labelled with these labels.}
  \item{nFeatureValues}{Deprecated. This is a legacy argument, retained for
    backwards compatibility. Any value passed to it will have no effect.}
  \item{timePoints}{An array of length \code{nFeatures}, containing the time
    points of the measurements.}
  \item{dataType}{A string specifying the data type. Either
    \code{``multinomial''}, \code{``time-course''}, or \code{``cubicspline''}.}
  \item{noise}{Noise term for each gene, required only if
    \env{noiseMode=2.} The noise term for each gene is calculated as
    \deqn{\frac{\sum(\mathrm{residuals}^2)}{(\mathrm{number\, of\,
	  observations\, for\, gene} - 1)(\mathrm{number\, of\,
	  replicates})},} where (number of observations for
      gene) is typically (number of time points * number of replicates).}
  \item{numReps}{Number of replicates per observation.}
  \item{noiseMode}{Noise mode. If 0 then fitted noise; 2 estimated
    noise from replicates.}
  \item{robust}{0 to use single Gaussian likelihood, 1 to use mixture
    likelihood.}
  \item{numThreads}{The BHC library has been parallelised using
    OpenMP (currently on UN*X systems only). Specify here the number of
    threads to use (the default value is 1).}
  \item{randomised}{Set to TRUE if you wish to use the randomised
    algorithm.}
  \item{m}{If randomised is set to TRUE, then this is the dimension of
    the randomly chosen subset D_m in the randomised algorithm.}
  \item{verbose}{If set to TRUE, the algorithm will output
    some information to screen as it runs.}
}
%\details{}
\value{
  A DENDROGRAM object (see the R stats package for details).
}
\references{\emph{Bayesian Hierarchical Clustering}, Heller + Ghahramani, Gatsby Unit Technical Report GCNU-TR
  2005-002 (2005); also see shorter version in ICML-2005;
\emph{R/BHC:fast Bayesian hierarchical clustering for microarray data}, Savage et
  al, BMC Bioinformatics 10:242 (2009); 
\emph{Bayesian hierarchical clustering for microarray time series data with replicates and outlier measurements}, Cooke et al, currently under review
}
\author{
  Rich Savage, Emma Cooke, Robert Darkins, and Yang Xu
}
%\note{}
\seealso{\code{\link{hclust}}}
\examples{
##BUILD SAMPLE DATA AND LABELS
data         <- matrix(0,15,10)
itemLabels   <- vector("character",15)
data[1:5,]   <- 1 ; itemLabels[1:5]   <- "a"
data[6:10,]  <- 2 ; itemLabels[6:10]  <- "b"
data[11:15,] <- 3 ; itemLabels[11:15] <- "c"
timePoints   <- 1:10 # for the time-course case

##DATA DIMENSIONS
nDataItems <- nrow(data)
nFeatures  <- ncol(data)

##RUN MULTINOMIAL CLUSTERING
hc1 <- bhc(data, itemLabels, verbose=TRUE)
plot(hc1, axes=FALSE)

##RUN TIME-COURSE CLUSTERING
hc2 <- bhc(data, itemLabels, 0, timePoints, "time-course",
          numReps=1, noiseMode=0, numThreads=2, verbose=TRUE)
plot(hc2, axes=FALSE)

##OUTPUT CLUSTER LABELS TO FILE
WriteOutClusterLabels(hc1, "labels.txt", verbose=TRUE)

##FOR THE MULTINOMIAL CASE, THE DATA CAN BE DISCRETISED
newData      <- data[] + rnorm(150, 0, 0.1);
percentiles  <- FindOptimalBinning(newData, itemLabels, transposeData=TRUE, verbose=TRUE)
discreteData <- DiscretiseData(t(newData), percentiles=percentiles)
discreteData <- t(discreteData)
hc3          <- bhc(discreteData, itemLabels, verbose=TRUE)
plot(hc3, axes=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
