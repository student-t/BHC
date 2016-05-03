/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <limits>

#include "SquaredExponentialTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

/* ---------------------------------------------------------------------- */

SquaredExponentialTimecourseDataSet::SquaredExponentialTimecourseDataSet() {}

/* ---------------------------------------------------------------------- */

SquaredExponentialTimecourseDataSet::SquaredExponentialTimecourseDataSet(string dataFile)
{
  ReadInData(dataFile);
}

/* ---------------------------------------------------------------------- */

SquaredExponentialTimecourseDataSet::
SquaredExponentialTimecourseDataSet(const vector<vector<double> >& inputData)
{
  data = inputData;

  // Find the data size
  nDataItems  = data.size();
  nFeatures   = data[0].size();
  nTimePoints = nFeatures;

  //PRINT OUT SOME USEFUL INFORMATION
  //cout << "----------" << endl;
  //cout << "nTimeSeries: " << nDataItems << endl;
  //cout << "nTimePoints: " << nTimePoints << endl;
  //cout << "----------" << endl;
}

/* ----------------------------------------------------------------------
   Compute the log-evidence for a single cluster containing the data
   items identified by 'itemIndex'.
------------------------------------------------------------------------- */

double SquaredExponentialTimecourseDataSet::
SingleClusterLogEvidence(const vector<int>& itemIndex,
			 double& lengthScale,
			 double& noiseFreeScale,
			 double& noiseSigma,
			 double& mixtureComponent)
{
  // Declarations
  int i, j, index;
  const int nCurrentItems=itemIndex.size();
  double replicateNoise,logEvidence=-numeric_limits<double>::infinity();
  vector<double> yValues=vector<double>(nCurrentItems*nTimePoints);

  // Extract the data points for this current cluster;
  // store the relevant data in yValues in row-major order
  for (i=0; i<nCurrentItems; i++)
  {
    index=i;
    for (j=0; j<nTimePoints; j++)
    {
      yValues[index]=data[itemIndex[i]][j];
      index+=nCurrentItems;
    }
  }

  // Optimise the hyperparameters (length-scale, noise-free-scale, noise-sigma)
  if (noise_mode == 0 )
  {
    OptimiseHyperparameters(yValues, lengthScale, noiseFreeScale, noiseSigma);
    logEvidence = ComputeMaximisedLogEvidence(yValues, lengthScale, noiseFreeScale, noiseSigma);
  }
  else if (noise_mode == 1)
  {
    noiseSigma = noiseData[0];
    OptimiseHyperparametersFixedNoise(yValues, lengthScale, noiseFreeScale, noiseSigma);
    logEvidence = ComputeMaximisedLogEvidence(yValues, lengthScale, noiseFreeScale, noiseSigma);
  }
  else if (noise_mode == 2)
  {
    replicateNoise = GetClusterSEMS(itemIndex);
    OptimiseHyperparametersEstimatedNoise(yValues, lengthScale, noiseFreeScale, noiseSigma, replicateNoise);
    logEvidence = ComputeMaximisedLogEvidence(yValues, lengthScale, noiseFreeScale, noiseSigma);
  }
  else
  {
    cout << "Error: noise_mode not recogised" <<endl;
  }
  return logEvidence;
}

/* ---------------------------------------------------------------------- */

double  SquaredExponentialTimecourseDataSet::
ComputeMaximisedLogEvidence(vector<double> yValues,
			    double& lengthScale,
			    double& noiseFreeScale,
			    double& noiseSigma)
{
  // Declarations
  int blockSize;
  BlockCovarianceMatrix covarFunction;
  double logEvidence;
  
  // Find block size
  blockSize = yValues.size() / nTimePoints;

  // Calculate maximised log-evidence
  covarFunction = SquareExponentialCovarianceFunction(lengthScale,
						      blockSize, 
						      noiseFreeScale);
  covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
  logEvidence = ComputeLogEvidence(covarFunction, yValues);

  return logEvidence;
}

/* ----------------------------------------------------------------------
   Optimise the hyperparameters for the GP model, returning the optimal
   log-evidence. The yValues should contain distinct time points.
   This uses a quasi-Newton method (implemented in TimecourseDataSet)
   to optimise the hyperparameters:
      - lengthScale (of the SE covariance function)
      - noiseFreeScale (the amplitude of the SE covariance function term)
      - noiseSigma
------------------------------------------------------------------------- */

void SquaredExponentialTimecourseDataSet::
OptimiseHyperparameters(const vector<double>& yValues,
			double& lengthScale,
			double& noiseFreeScale,
			double& noiseSigma)
{
  int i;
  int blockSize = yValues.size() / nTimePoints;
  double bestLengthScale;
  double bestLogEv, trialLogEv;
  vector<double> params(3);
  vector<int> fix(3);
  
  // Guess a starting state using a coarse-grained method
  bestLengthScale=2.0;
  params[1]=1.0;
  params[2]=0.5;
  bestLogEv=-numeric_limits<double>::infinity();
  for(i=2; i<=10; i+=2)
    {
      params[0] = static_cast<double>(i);
      trialLogEv = ComputeLogEvidenceFromHyperparameters(yValues,
							 blockSize,
							 params);
      if(trialLogEv > bestLogEv)
	{
	  bestLengthScale = params[0];
	  bestLogEv = trialLogEv;
	}
    }
  params[0]=bestLengthScale;

  // Now do the actual maximisation
  fix[0]=fix[1]=fix[2]=0; // do not fix any params
  double fret=0;
  double gtol=fast_switch ? (1.0e-1) : (1.0e-2); // the convergence tolerance
  DFPMaximise(params,fix,gtol,fret,blockSize,yValues);
  
  // Return the result
  lengthScale=params[0];
  noiseFreeScale=params[1];
  noiseSigma=params[2];
}

/* ---------------------------------------------------------------------- */

void SquaredExponentialTimecourseDataSet::
ImposeConstraintsOnHyperparameters(vector<double>& params)
{
  params[0] = MAX(params[0], 0.2); // lengthScale
  params[1] = MAX(params[1], 0.2); // noiseFreeScale
  params[2] = MAX(params[2], 0.05); // noiseSigma
  params[2] = MIN(params[2], 1.0); // noiseSigma
}

/* ---------------------------------------------------------------------- */

double SquaredExponentialTimecourseDataSet::
ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
				      const int blockSize,
				      const vector<double>& params)
{
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  return ComputeLogEvidence(covarFunc, yValues);
}

/* ---------------------------------------------------------------------- */

void SquaredExponentialTimecourseDataSet::
ComputeGradientsFromHyperparameters(const vector<double>& yValues,
				    const int blockSize,
				    const vector<double>& params,
				    vector<double>& grad)
{
  // Build the covariance function for these hyperparameters
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  // Useful intermediate values for computing the gradient
  covarFunc.InvertMatrix(); // K = K^-1
  vector<double> alpha = covarFunc.VectorMultiply(yValues); // K^-1 * y
  BlockCovarianceMatrix covarDeriv_ls =
    SquareExponential_lengthDerivative(params[0],
				       blockSize,
				       params[1]);
  BlockCovarianceMatrix covarDeriv_nf = SquareExponentialCovarianceFunction(params[0],
									    blockSize,
									    1);
  // Compute the gradient at this point in hyperspace
  grad[0] = ComputeGradient(covarFunc, covarDeriv_ls, alpha);
  grad[1] = ComputeGradient(covarFunc, covarDeriv_nf, alpha);
  grad[2] = ComputeNoiseGradient(covarFunc, alpha, params[2]);
}

/* ---------------------------------------------------------------------- */

void SquaredExponentialTimecourseDataSet::
ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						  const int blockSize,
						  const vector<double>& params,
						  double& logEv,
						  vector<double>& grad)
{
  // Build the covariance function for these hyperparameters
  BlockCovarianceMatrix covarFunc = 
    AddNoiseToCovarianceFunction(SquareExponentialCovarianceFunction(params[0],
								     blockSize,
								     params[1]),
				 params[2]);
  // Compute the log-evidence
  logEv = ComputeLogEvidence(covarFunc, yValues);
  
  // Useful intermediate values for computing the gradient
  covarFunc.InvertMatrix(); // K = K^-1
  vector<double> alpha = covarFunc.VectorMultiply(yValues); // K^-1 * y
  BlockCovarianceMatrix covarDeriv_ls =
    SquareExponential_lengthDerivative(params[0],
				       blockSize,
				       params[1]);
  BlockCovarianceMatrix covarDeriv_nf = SquareExponentialCovarianceFunction(params[0],
									    blockSize,
									    1);
  // Compute the gradient at this point in hyperspace
  grad[0] = ComputeGradient(covarFunc, covarDeriv_ls, alpha);
  grad[1] = ComputeGradient(covarFunc, covarDeriv_nf, alpha);
  grad[2] = ComputeNoiseGradient(covarFunc, alpha, params[2]);
}

/* ---------------------------------------------------------------------- */

double SquaredExponentialTimecourseDataSet::GetMLIINoise(vector<int> itemIndex)
{
  //DECLARATIONS
  int                      i, j;
  int                      index, counter=0;
  int                      nCurrentItems=itemIndex.size();
  vector<double>           extractedData, yValues;
  double				   fittedNoise;
  //EXTRACT THE DATA POINTS FOR THIS CURRENT CLUSTER
  //data vector is nDataItems*nTimePoints
  for (i=0; i<nCurrentItems; i++)
  {
    for (j=0; j<nTimePoints; j++)
    {
      index = itemIndex[i];
      extractedData.push_back(data[index][j]);//store the corresponding data value
    }
  }
  //RE-CONSTRUCT THE DATA VECTOR SO IT HAS THE CORRECT ORDERING FOR THE BlockCovarianceMatrix FORM WE'LL BE USING
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

double SquaredExponentialTimecourseDataSet::CalculateFittedNoiseHyperparameter(const vector<double>& yValues)
{
  double lengthScale, noiseFreeScale, noiseSigma;
  OptimiseHyperparameters(yValues,lengthScale,noiseFreeScale,noiseSigma);
  return noiseSigma;
}

/* ---------------------------------------------------------------------- */

void  SquaredExponentialTimecourseDataSet::
OptimiseHyperparametersEstimatedNoise(vector<double> yValues,
				      double& lengthScale,
				      double& noiseFreeScale,
				      double& noiseSigma,
				      double replicateNoise)
{
  // Declarations
  int i;
  int blockSize;
  int nIterations=50;
  double minLengthScale = 0.2,
    minNoiseFreeScale=0.1,
    minNoiseSigma=replicateNoise,
    maxNoiseSigma=1.2;
  double gradScale = 0.01;
  double maxNoiseGrad, maxNoiseStep=0.01;//if the noise steps get too big, the fit can become unstable
  double currentLengthScale;
  double currentLogEvidence, gridLogEvidence,
    logEvidence=-numeric_limits<double>::infinity();
  double grad_lengthScale, grad_noiseFreeScale, grad_noiseSigma;
  BlockCovarianceMatrix covarFunction, covarDerivative_ls, covarDerivative_nf, covarCopy, workingMatrix;//derivative WRT length-scale
  vector<double> alpha, vec1, vec2, params;
  lengthScale=2.0;
  noiseFreeScale=1.0;
  noiseSigma=replicateNoise;
  
  // Find block size
  blockSize = yValues.size() / nTimePoints;
  
  // Find gamma parameters
  if (replicateNoise < 1.0)
  {
    params = OptimiseGammaParams(replicateNoise);
  }

  //USE A COARSE GRID TO FIND A REASONABLE STARTING POSITION FOR lengthScale
  //this is useful because the length-scale is slow to converge
  for (i=0; i<10; i++)
  {
    currentLengthScale    = 1*(i+1);
    covarFunction         = SquareExponentialCovarianceFunction(currentLengthScale, blockSize, noiseFreeScale);
    covarFunction         = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
    currentLogEvidence    = ComputeLogEvidence(covarFunction, yValues);
    //  cout << "NoiseFreeScale: " << currentNoiseFreeScale << "\tLogEv: " << currentLogEvidence << endl;
    //check whether this is the best value so far
    if (currentLogEvidence>logEvidence)
    {
      logEvidence    = currentLogEvidence;
      lengthScale    = currentLengthScale;
    }
  }
  gridLogEvidence = logEvidence;//store a value for debugging
  //USE GRADIENT ASCENT TO FIND THE OPTIMAL HYPERPARAMETER VALUES
  for (i=0; i<nIterations; i++)
  {
    //find the covariance function for this pass
    covarFunction = SquareExponentialCovarianceFunction(lengthScale, blockSize, noiseFreeScale);
    covarFunction = AddNoiseToCovarianceFunction(covarFunction, noiseSigma);
    //find useful intermediate values with which to construct the gradients
    covarFunction.InvertMatrix(); //changing the covarFunction
    alpha              = covarFunction.VectorMultiply(yValues);
    covarDerivative_ls = SquareExponential_lengthDerivative(lengthScale, blockSize, noiseFreeScale);
    covarDerivative_nf = SquareExponentialCovarianceFunction(lengthScale, blockSize, 1);
    //compute the gradients
    grad_lengthScale    = ComputeGradient(covarFunction, covarDerivative_ls, alpha);
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
    lengthScale    += 10 * gradScale * grad_lengthScale;//speed-up the change of length scale here with a hardwired factor of 10
    noiseFreeScale += gradScale * grad_noiseFreeScale;
    noiseSigma     += gradScale * grad_noiseSigma;
    //enforce priors of minimum parameter values
    if (lengthScale<minLengthScale)       lengthScale    = minLengthScale;
    if (noiseFreeScale<minNoiseFreeScale) noiseFreeScale = minNoiseFreeScale;
    if (noiseSigma>maxNoiseSigma)         noiseSigma     = maxNoiseSigma;
    if (noiseSigma<minNoiseSigma)         noiseSigma     = minNoiseSigma;
  }
}

/* ----------------------------------------------------------------------
   Optimise the hyperparameters for the GP model, returning the optimal
   log-evidence in the case where noise is fixed for whole clustering.
---------------------------------------------------------------------- */

void SquaredExponentialTimecourseDataSet::
OptimiseHyperparametersFixedNoise(vector<double> yValues,
				  double& lengthScale,
				  double& noiseFreeScale,
				  double& noiseSigma)
{
  int i;
  int blockSize = yValues.size() / nTimePoints;
  double bestLengthScale;
  double bestLogEv, trialLogEv;
  vector<double> params(3);
  vector<int> fix(3);

  // Guess a starting state using a coarse-grained method
  bestLengthScale=2.0;
  params[1]=1.0;
  params[2]=0.5;
  bestLogEv=-numeric_limits<double>::infinity();
  for(i=2; i<=10; i+=2)
    {
      params[0] = static_cast<double>(i);
      trialLogEv = ComputeLogEvidenceFromHyperparameters(yValues,
							 blockSize,
							 params);
      if(trialLogEv > bestLogEv)
	{
	  bestLengthScale = params[0];
	  bestLogEv = trialLogEv;
	}
    }
  params[0]=bestLengthScale;

  // Do the maximisation
  fix[0]=fix[1]=0;
  fix[2]=1; // fix noiseSigma
  double fret=0;
  double gtol=fast_switch ? (1.0e-1) : (1.0e-2); // the convergence tolerance
  DFPMaximise(params,fix,gtol,fret,blockSize,yValues);
  
  // Return the result
  lengthScale=params[0];
  noiseFreeScale=params[1];
  noiseSigma=params[2];
}

/* ----------------------------------------------------------------------
   Compute a noise-less square exponential (SE) covariance function.
---------------------------------------------------------------------- */

BlockCovarianceMatrix SquaredExponentialTimecourseDataSet::
SquareExponentialCovarianceFunction(double lengthScale,
				    int blockSize,
				    double noiseFreeScale)
{
  // Declarations
  int                    i, j;
  double                 covarElement;
  BlockCovarianceMatrix  blockMatrix;

  // Initialise the block matrix
  blockMatrix.nRank     = nTimePoints;
  blockMatrix.blockSize = blockSize;

  // Initialise the covariance function
  blockMatrix.noisyCoeff=vector<double>(nTimePoints, 0.0);
  blockMatrix.noiseFreeCoeff=
    vector<vector<double> >(nTimePoints, vector<double>(nTimePoints, 0.0));

  // Compute each element of the covariance function
  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      covarElement = fabs(timePoints[i] - timePoints[j]);
      covarElement = covarElement*covarElement;
      covarElement /= 2*lengthScale*lengthScale;
      covarElement *= -1;
      covarElement  = exp(covarElement);
      covarElement *= noiseFreeScale;
      // and store in 2 elements (covariance is symmetric)
      // this duplicates effort for i==j; not a big deal, computationally :-)
      blockMatrix.noiseFreeCoeff[i][j] = covarElement;
      blockMatrix.noiseFreeCoeff[j][i] = covarElement;
    }
  }
  
  return blockMatrix;
}

/* ----------------------------------------------------------------------
   Compute the partial derivative w.r.t. length-scale of a square
   exponential (SE) covariance function. Note that we've hard-coded
   noise=0 here (a bit ugly).
---------------------------------------------------------------------- */

BlockCovarianceMatrix SquaredExponentialTimecourseDataSet::
SquareExponential_lengthDerivative(double lengthScale,
				   int blockSize,
				   double noiseFreeScale)
{
  // Declarations
  int i, j;
  double covarElement, deltaTime;
  BlockCovarianceMatrix  blockMatrix;

  assert(lengthScale > 0);
  assert(blockSize > 0);

  // Initialise the block matrix
  blockMatrix.nRank = nTimePoints;
  blockMatrix.blockSize = blockSize;

  // Initialise the covariance function
  blockMatrix.noiseFreeCoeff=
    vector<vector<double> >(nTimePoints, vector<double>(nTimePoints, 0.0));
  blockMatrix.noisyCoeff = vector<double>(nTimePoints, 0.0);

  // Compute each element of the covariance function
  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      deltaTime     = fabs(timePoints[i] - timePoints[j]);
      covarElement  = deltaTime*deltaTime;
      covarElement /= 2*lengthScale*lengthScale;
      covarElement *= -1;
      covarElement  = exp(covarElement);
      covarElement *= deltaTime*deltaTime;
      covarElement /= lengthScale*lengthScale*lengthScale;
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
