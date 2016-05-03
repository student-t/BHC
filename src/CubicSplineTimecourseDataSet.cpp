/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <limits>

#include "CubicSplineTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

/* ----------------------------------------------------------------------
   Default constructor.
---------------------------------------------------------------------- */

CubicSplineTimecourseDataSet::CubicSplineTimecourseDataSet() {}

/* ---------------------------------------------------------------------- */

CubicSplineTimecourseDataSet::CubicSplineTimecourseDataSet(string dataFile)
{
  ReadInData(dataFile);
}

/* ---------------------------------------------------------------------- */

CubicSplineTimecourseDataSet::
CubicSplineTimecourseDataSet(const vector<vector<double> >& inputData)
{
  // Copy the data into the object
  data = inputData;

  // Find the data size
  nDataItems = data.size();
  nFeatures = data[0].size();
  nTimePoints = nFeatures;

  // Print out some useful information
  //cout << "----------" << endl;
  //cout << "nTimeSeries: " << nDataItems << endl;
  //cout << "nTimePoints: " << nTimePoints << endl;
  //cout << "----------" << endl;
}

/* ----------------------------------------------------------------------
   Compute the log-evidence for a single cluster containing the data items
   identified by itemIndex.
---------------------------------------------------------------------- */
 
double CubicSplineTimecourseDataSet::
SingleClusterLogEvidence(const vector<int>& itemIndex,
			 double& lengthScale,
			 double& noiseFreeScale,
			 double& noiseSigma,
			 double& mixtureComponent)
{
  // Declarations
  int i, j, index;
  const int nCurrentItems=itemIndex.size();
  double logEvidence=-numeric_limits<double>::infinity(), replicateNoise;
  vector<double> yValues=vector<double>(nCurrentItems*nTimePoints);

  // Extract the data points for this current cluster;
  // store the relevant data in row-major order in yValues
  for (i=0; i<nCurrentItems; i++)
  {
    index=i;
    for (j=0; j<nTimePoints; j++)
    {
      yValues[index]=data[itemIndex[i]][j];
      index+=nCurrentItems;
    }
  }

  // Optimise the hyperparameters (length-scale)
  if (noise_mode == 0)
  {
    OptimiseHyperparameters(yValues, noiseFreeScale, noiseSigma);
    logEvidence = ComputeMaximisedLogEvidence(yValues, noiseFreeScale, noiseSigma);
  }
  else if (noise_mode == 1)
  {
    noiseSigma = noiseData[0];
    OptimiseHyperparametersFixedNoise(yValues, noiseFreeScale, noiseSigma);
    logEvidence = ComputeMaximisedLogEvidence(yValues, noiseFreeScale, noiseSigma);
  }
  else if (noise_mode ==2)
  {
    replicateNoise = GetClusterSEMS(itemIndex);
    OptimiseHyperparametersEstimatedNoise(yValues, noiseFreeScale, noiseSigma, replicateNoise);
    logEvidence = ComputeMaximisedLogEvidence(yValues, noiseFreeScale, noiseSigma);
  }
  else
  {
    cout << "Error! noise_mode not recognised" <<endl;
  }

  return logEvidence;
}

/* ---------------------------------------------------------------------- */

double CubicSplineTimecourseDataSet::ComputeMaximisedLogEvidence(const vector<double>& yValues, double& noiseFreeScale, double& noiseSigma)
{
  // Declarations
  int blockSize;
  BlockCovarianceMatrix covarFunction;
  double logEvidence;

  // Find block size
  blockSize = yValues.size() / nTimePoints;

  // Calculate maximised log evidence
  covarFunction = CubicSplineCovarianceFunction(blockSize, noiseFreeScale);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
  logEvidence = ComputeLogEvidence(covarFunction, yValues);
  
  return logEvidence;
}

/* ----------------------------------------------------------------------
   Optimise the hyperparameters for the GP model, returning the optimal log-evidence.
   
   xValues should contain the distinct time points (*not* duplicates).
   This method uses a simple implementation of a gradient ascent method.
   Hyperparameters are:
      - length-scale (of the SE covariance function)
      - the noiseFreeScale (ie. amplitude of the SE covariance function term)
      - noiseSigma

   ! If this method needs optimising, then you could get it to use
   ! the DFPMaximise() function in TimecourseDataSet. See the
   ! SquaredExponentialTimecourseDataSet class for an example of how.
---------------------------------------------------------------------- */
 
void CubicSplineTimecourseDataSet::
OptimiseHyperparameters(const vector<double>& yValues,
			double& noiseFreeScale,
			double& noiseSigma)
{
  //DECLARATIONS
  int i;
  int blockSize;
  const int nIterations=50;
  double minNoiseFreeScale=0.1, minNoiseSigma=0.05, maxNoiseSigma=1.0;
  noiseFreeScale=1, noiseSigma=0.5;
  double gradScale = 0.01;
  double maxNoiseGrad, maxNoiseStep=0.01;//if the noise steps get too big, the fit can become unstable
  double grad_noiseFreeScale, grad_noiseSigma;
  BlockCovarianceMatrix covarFunction, covarDerivative_nf, covarCopy, workingMatrix;//derivative WRT length-scale
  vector<double> alpha, vec1, vec2;
  
  //FIND BLOCK SIZE
  blockSize = yValues.size() / nTimePoints;
  
  //USE GRADIENT ASCENT TO FIND THE OPTIMAL HYPERPARAMETER VALUES
  for (i=0; i<nIterations; i++)
    {
      //***IS IT WORTH STEPPING ONE PARAM AT A TIME??***
      
      //find the covariance function for this pass
      covarFunction = CubicSplineCovarianceFunction(blockSize, noiseFreeScale);
      covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
      
      //find useful intermediate values with which to construct the gradients
      covarFunction.InvertMatrix();
      alpha = covarFunction.VectorMultiply(yValues);
      covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
      
      //compute the gradients
      grad_noiseFreeScale = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
      
      grad_noiseSigma = ComputeNoiseGradient(covarFunction, alpha, noiseSigma);
      
      //enforce a maxNoiseStep to keep the optimisation stable
      maxNoiseGrad = maxNoiseStep / gradScale;
      if (fabs(grad_noiseSigma) > maxNoiseGrad)
	{
	  grad_noiseSigma /= fabs(grad_noiseSigma);
	  grad_noiseSigma *= maxNoiseGrad;
	}
      //compute the parameter steps and apply them
      noiseFreeScale += gradScale * grad_noiseFreeScale;
      noiseSigma     += gradScale * grad_noiseSigma;
      //enforce priors of minimum parameter values
      if (noiseFreeScale<minNoiseFreeScale) noiseFreeScale = minNoiseFreeScale;
      if (noiseSigma<minNoiseSigma)         noiseSigma = minNoiseSigma;
      if (noiseSigma>maxNoiseSigma)         noiseSigma = maxNoiseSigma;
    }
}

/* ---------------------------------------------------------------------- */

double CubicSplineTimecourseDataSet::GetMLIINoise(const vector<int>& itemIndex)
{
  // Declarations
  int i, j;
  int index, counter=0;
  int nCurrentItems=itemIndex.size();
  vector<double> extractedData, yValues;
  double fittedNoise;

  // Extract the data points for this current cluster;
  // data vector is nDataItems*nTimePoints in size
  for (i=0; i<nCurrentItems; i++)
  {
    for (j=0; j<nTimePoints; j++)
    {
      index = itemIndex[i];
      extractedData.push_back(data[index][j]);//store the corresponding data value
    }
  }

  // Re-construct the data vector so it has the correct ordering for
  // the BlockCovarianceMatrix form we'll be using
  yValues = extractedData;
  counter = 0;
  for (i=0; i<nCurrentItems; i++)
  {
    for (j=0; j<nTimePoints; j++)
    {
      index          = j*nCurrentItems + i;
      yValues[index] = extractedData[counter];
      counter++;
    }
  }
  //Optimise the hyperparameters and return the fitted noise
  fittedNoise = CalculateFittedNoiseHyperparameter(yValues);
  
  return fittedNoise;
}

/* ---------------------------------------------------------------------- */

double CubicSplineTimecourseDataSet::
CalculateFittedNoiseHyperparameter(const vector<double>& yValues)
{
  //DECLARATIONS
  int i;
  int blockSize;
  const int nIterations=50;
  double minNoiseFreeScale=0.1, minNoiseSigma=0.05, maxNoiseSigma=1.0;
  double noiseFreeScale=1, noiseSigma=0.5;
  double gradScale = 0.01;
  double maxNoiseGrad, maxNoiseStep=0.01;//if the noise steps get too big, the fit can become unstable
  double grad_noiseFreeScale, grad_noiseSigma;
  BlockCovarianceMatrix covarFunction, covarDerivative_nf, covarCopy, workingMatrix;//derivative WRT length-scale
  vector<double> alpha, vec1, vec2;
  //FIND BLOCK SIZE
  blockSize = yValues.size() / nTimePoints;

  //USE GRADIENT ASCENT TO FIND THE OPTIMSL HYPERPARAMETER VALUES
  for (i=0; i<nIterations; i++)
  {
    //***IS IT WORTH STEPPING ONE PARAM AT A TIME??***
    //find the covariance function for this pass
    covarFunction = CubicSplineCovarianceFunction(blockSize, noiseFreeScale);
    covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
    //find useful intermediate values with which to construct the gradients
    covarFunction.InvertMatrix();
    alpha              = covarFunction.VectorMultiply(yValues);
    covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
    //compute the gradients
    grad_noiseFreeScale = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
    grad_noiseSigma     = ComputeNoiseGradient(covarFunction, alpha, noiseSigma);
    //enforce a maxNoiseStep to keep the optimisation stable
    maxNoiseGrad = maxNoiseStep / gradScale;
    if (fabs(grad_noiseSigma) > maxNoiseGrad)
    {
      grad_noiseSigma /= fabs(grad_noiseSigma);
      grad_noiseSigma *= maxNoiseGrad;
    }
    //compute the parameter steps and apply them
    noiseFreeScale += gradScale * grad_noiseFreeScale;
    noiseSigma     += gradScale * grad_noiseSigma;
    //enforce priors of minimum parameter values
    if (noiseFreeScale<minNoiseFreeScale) noiseFreeScale = minNoiseFreeScale;
    if (noiseSigma<minNoiseSigma)         noiseSigma     = minNoiseSigma;
    if (noiseSigma>maxNoiseSigma)         noiseSigma     = maxNoiseSigma;
  }

  return(noiseSigma);
}

/* ---------------------------------------------------------------------- */

void CubicSplineTimecourseDataSet::OptimiseHyperparametersEstimatedNoise(const vector<double>& yValues, double& noiseFreeScale, double& noiseSigma, double replicateNoise)
{
  //DECLARATIONS
  int i;
  int blockSize;
  const int nIterations=50;
  double minNoiseFreeScale=0.1, minNoiseSigma=replicateNoise, maxNoiseSigma=1.2;
  noiseFreeScale=1, noiseSigma=0.5;
  double gradScale = 0.01;
  double maxNoiseGrad, maxNoiseStep=0.01;//if the noise steps get too big, the fit can become unstable
  double grad_noiseFreeScale, grad_noiseSigma;
  BlockCovarianceMatrix covarFunction, covarDerivative_nf, covarCopy, workingMatrix;//derivative WRT length-scale
  vector<double> alpha, vec1, vec2, params;
  
  //FIND BLOCK SIZE
  blockSize = yValues.size() / nTimePoints;
  
  //FIND THE GAMMA PARAMETERS
  if (replicateNoise < 1.0)
  {
    params = OptimiseGammaParams(replicateNoise);
  }
  
  //USE GRADIENT ASCENT TO FIND THE OPTIMSL HYPERPARAMETER VALUES
  for (i=0; i<nIterations; i++)
  {
    //***IS IT WORTH STEPPING ONE PARAM AT A TIME??***
    //find the covariance function for this pass
    covarFunction = CubicSplineCovarianceFunction(blockSize, noiseFreeScale);
    covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
    //find useful intermediate values with which to construct the gradients
    covarFunction.InvertMatrix();
    alpha              = covarFunction.VectorMultiply(yValues);
    covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
    //compute the gradients
    grad_noiseFreeScale = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
    if (replicateNoise < 1.0)
    {
      grad_noiseSigma     = ComputeNoiseGradientIncludingGammaHyperprior(covarFunction, alpha, noiseSigma, params);
    }
    else
    {
      grad_noiseSigma = ComputeNoiseGradient(covarFunction, alpha, noiseSigma);
    }
    //enforce a maxNoiseStep to keep the optimisation stable
    maxNoiseGrad = maxNoiseStep / gradScale;
    if (fabs(grad_noiseSigma) > maxNoiseGrad)
    {
      grad_noiseSigma /= fabs(grad_noiseSigma);
      grad_noiseSigma *= maxNoiseGrad;
    }
    //compute the parameter steps and apply them
    noiseFreeScale += gradScale * grad_noiseFreeScale;
    noiseSigma     += gradScale * grad_noiseSigma;
    //enforce priors of minimum parameter values
    if (noiseFreeScale<minNoiseFreeScale) noiseFreeScale = minNoiseFreeScale;
    if (noiseSigma>maxNoiseSigma)         noiseSigma     = maxNoiseSigma;
    if (noiseSigma<minNoiseSigma)         noiseSigma     = minNoiseSigma;
  }
}

/* ---------------------------------------------------------------------- */

void CubicSplineTimecourseDataSet::OptimiseHyperparametersFixedNoise(const vector<double>& yValues, double& noiseFreeScale, double& noiseSigma)
{
  //DECLARATIONS
  int i;
  int blockSize;
  const int nIterations=50;
  double minNoiseFreeScale=0.1;
  noiseFreeScale=1;
  double gradScale = 0.01;
  double grad_noiseFreeScale;
  BlockCovarianceMatrix covarFunction, covarDerivative_nf, covarCopy, workingMatrix;//derivative WRT length-scale
  vector<double> alpha, vec1, vec2;
  
  //FIND BLOCK SIZE
  blockSize = yValues.size() / nTimePoints;

  //USE GRADIENT ASCENT TO FIND THE OPTIMSL HYPERPARAMETER VALUES
  for (i=0; i<nIterations; i++)
  {
    //***IS IT WORTH STEPPING ONE PARAM AT A TIME??***
    //find the covariance function for this pass
    covarFunction = CubicSplineCovarianceFunction(blockSize, noiseFreeScale);
    covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
    //find useful intermediate values with which to construct the gradients
    covarFunction.InvertMatrix();
    alpha              = covarFunction.VectorMultiply(yValues);
    covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
    //compute the gradients
    grad_noiseFreeScale = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
    //compute the parameter steps and apply them
    noiseFreeScale += gradScale * grad_noiseFreeScale;
    //enforce priors of minimum parameter values
    if (noiseFreeScale<minNoiseFreeScale) noiseFreeScale = minNoiseFreeScale;
  }
}

/* ---------------------------------------------------------------------- 
   Compute a noise-less square exponential (SE) covariance function
   changes for cubic spline.
---------------------------------------------------------------------- */

BlockCovarianceMatrix CubicSplineTimecourseDataSet::
CubicSplineCovarianceFunction(int blockSize, double noiseFreeScale)
{
  // Declarations
  int i, j;
  double covarElement, mintimePoint, partcovarElement;
  BlockCovarianceMatrix  blockMatrix;

  // Initialise the block matrix
  blockMatrix.nRank = nTimePoints;
  blockMatrix.blockSize = blockSize;

  // Initialise the covariance function
  for (i=0; i<nTimePoints; i++)
  {
    blockMatrix.noiseFreeCoeff.push_back(vector<double>(nTimePoints, 0));
    blockMatrix.noisyCoeff.push_back(0);
  }
  
  // Compute each element of the covariance function
  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      if (timePoints[i] > timePoints[j])
        mintimePoint = timePoints[j];
      else
        mintimePoint = timePoints[i];
      covarElement = fabs(timePoints[i] - timePoints[j]);
      covarElement *= 0.5*pow(mintimePoint, 2);
      partcovarElement = pow(mintimePoint, 3);
      partcovarElement /= 3.0;
      covarElement += partcovarElement;
      covarElement *= noiseFreeScale;
      
      // and store in 2 elements (covariance is symmetric)
      // this duplicates effort for i==j; not a big deal, computationally :-)
      blockMatrix.noiseFreeCoeff[i][j] = covarElement;
      blockMatrix.noiseFreeCoeff[j][i] = covarElement;
    }
  }
  return blockMatrix;
}

/* ---------------------------------------------------------------------- */

double CubicSplineTimecourseDataSet::
ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
				      const int blockSize,
				      const vector<double>& params)
{
  BlockCovarianceMatrix covarFunction;
  covarFunction = CubicSplineCovarianceFunction(blockSize, params[0]);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, params[1]);
  return ComputeLogEvidence(covarFunction, yValues);
}

/* ---------------------------------------------------------------------- */

// params[0] = noiseFreeScale, params[1] = noiseSigma
void CubicSplineTimecourseDataSet::
ComputeGradientsFromHyperparameters(const vector<double>& yValues,
				    const int blockSize,
				    const vector<double>& params,
				    vector<double>& grad)
{
  BlockCovarianceMatrix covarFunction, covarDerivative_nf;
  vector<double> alpha;
  
  covarFunction = CubicSplineCovarianceFunction(blockSize, params[0]);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, params[1]);
  
  //find useful intermediate values with which to construct the gradients
  covarFunction.InvertMatrix();
  alpha = covarFunction.VectorMultiply(yValues);
  covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
  
  //compute the gradients
  grad[0] = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
  grad[1] = ComputeNoiseGradient(covarFunction, alpha, params[1]);
}

/* ---------------------------------------------------------------------- */

void CubicSplineTimecourseDataSet::
ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						  const int blockSize,
						  const vector<double>& params,
						  double& logEv,
						  vector<double>& grad)
{
  BlockCovarianceMatrix covarFunction, covarDerivative_nf;
  vector<double> alpha;
  
  covarFunction = CubicSplineCovarianceFunction(blockSize, params[0]);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, params[1]);

  // Compute the log-evidence
  logEv = ComputeLogEvidence(covarFunction, yValues);
  
  //find useful intermediate values with which to construct the gradients
  covarFunction.InvertMatrix();
  alpha = covarFunction.VectorMultiply(yValues);
  covarDerivative_nf = CubicSplineCovarianceFunction(blockSize, 1);
  
  //compute the gradients
  grad[0] = ComputeGradient(covarFunction, covarDerivative_nf, alpha);
  grad[1] = ComputeNoiseGradient(covarFunction, alpha, params[1]);
}

/* ---------------------------------------------------------------------- */

void CubicSplineTimecourseDataSet::
ImposeConstraintsOnHyperparameters(vector<double>& params) {}

/* ---------------------------------------------------------------------- */
