/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include "TimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"
#include <limits>
#include <numeric>

/* ---------------------------------------------------------------------- */

TimecourseDataSet::TimecourseDataSet()
{}

/* ---------------------------------------------------------------------- */

TimecourseDataSet::TimecourseDataSet(string dataFile)
{
  ReadInData(dataFile);
}

/* ----------------------------------------------------------------------
   Constructor: read in data from input vector; need this in order
   to link to R.
---------------------------------------------------------------------- */

TimecourseDataSet::TimecourseDataSet(const vector<vector<double> >& inputData)
{
  int i;
  
  //COPY THE DATA INTO THE OBJECT
  data = inputData;
  cout << "Need to deal with the time points here test!" << endl;

  //FIND THE DATA SIZE
  nDataItems  = data.size();
  nFeatures   = data[0].size();
  nTimePoints = nFeatures;

  for (i=0; i<nTimePoints; i++)
  {
    timePoints.push_back(i); // what is this all about?
  }

  //PRINT OUT SOME USEFUL INFORMATION
  ///cout << "----------" << endl;
  ///cout << "nTimeSeries: " << nDataItems << endl;
  ///cout << "nTimePoints: " << nTimePoints << endl;
  ///cout << "----------" << endl;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::ReadInData(string dataFile)
{
  //DECLARATIONS
  int     i, j;
  double  inputValue;
  fstream file;
  
  //FIND THE SIZE OF THE INPUT DATA; WE WILL ACTUALLY USE nTimePoints
  //INSTEAD OF nFeatures HERE (this makes the code easier to read)
  FindDataSize(dataFile);
  nTimePoints = nFeatures;
  nDataItems--;//we must subtract one here to account for the header
  
  //OPEN DATA FILE (WITH ERROR-CHECKING)
  file.open (dataFile.c_str(), fstream::in);//dataFile is a string object, so need to convert to a char array here!
  if(file.fail())
  {
    cout<<"Failed to open file "<<dataFile<<"."<<endl;
    if(!system("PAUSE"))/*didn't work*/;
  }
  
  //EITHER READ IN THE TIME-POINTS, OR GENERATE A UNIFORM SET
  for (i=0; i<nTimePoints; i++)
  {
    file >> inputValue;
    timePoints.push_back(inputValue);//store the latest time point read in from file
    //cout << timePoints[i] << endl;
  }

  //NORMALISE THE TIME POINTS TO SOMETHING SENSIBLE
  //want to make the min gap a unit spacing
  //this just makes the hyperparameter fitting more straightforward
  //actually, might be easier to handle this in the fitting routine (change the lengthScale)

  //READ THE DATA INTO OUR DATA VECTOR
  for (i=0; i<nDataItems; i++)
  {
    data.push_back(vector<double>(nTimePoints, 0.0));//add a vector for this item
    for (j=0; j<nTimePoints; j++)
    {
      file >> inputValue;
      data[i][j] = inputValue;
      //cout << "data[i][j]" << i << j << " " << data[i][j] << " " << endl;
    }
  }
  file.close();
  //PRINT OUT SOME USEFUL INFORMATION
  //cout << "----------" << endl;
  //cout << dataFile << endl;
  //cout << "nTimeSeries: " << nDataItems << endl;
  //cout << "nTimePoints: " << nTimePoints << endl;
  //cout << "----------" << endl;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::ReadInTimePoints(vector<double> timePoints_copy)
{
  //DECLARATIONS
  int i;
  
  for (i=0; i<nTimePoints; i++)
  {
    timePoints.push_back(timePoints_copy[i]);
  }
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::ReadInNoise(string dataFile)
{
  //DECLARATIONS
  int     i;
  double  inputValue;
  fstream file;
  
  //OPEN DATA FILE (WITH ERROR-CHECKING)
  file.open (dataFile.c_str(), fstream::in);//dataFile is a string object, so need to convert to a char array here!
  if(file.fail())
  {
    cout<<"Failed to open file "<<dataFile<<"."<<endl;
    if(!system("PAUSE"))/*didn't work*/;
  }
  i=0;
  while ( file >>inputValue )
  {
    i++;
    noiseData.push_back(inputValue);
  }
  file.close();
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::ReadInNoise(vector<double> noise_copy)
{
  noiseData = noise_copy;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::SetNoiseMode(int mode)
{
  noise_mode = mode;
  ///cout << "Noise mode = " << noise_mode << endl;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::SetReps(int num_reps)
{
  reps = num_reps;
  ///cout << "Num reps = " << reps << endl;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::SetDataType(string type)
{
  dataType = type;
  ///cout << "dataType:" << dataType << endl;
}

/* ---------------------------------------------------------------------- */

void TimecourseDataSet::SetRobustMode(int mode) // also sets the dataRange if Robust = 1
{
  //DECLARATIONS
  int i, j;
  vector<double> inputvals;
  vector<double>::iterator itmin, itmax;
  double min, max;

  robust_mode = mode;
  ///cout << "Robust mode = " << mode << endl;
  if (robust_mode == 1)
  {
    //extract the data values into vector
    for (i=0; i<nDataItems; i++)
    {
      for (j=0; j<nTimePoints; j++)
      {
        inputvals.push_back(data[i][j]);//store the corresponding data value
      }
    }
    itmin = min_element( inputvals.begin(), inputvals.end() );
    itmax = max_element( inputvals.begin(), inputvals.end() );
    min = *itmin;
    max = *itmax;
    dataRange = max-min; //This dataRange is constant for the dataSet
  }
}

/* ---------------------------------------------------------------------- */

double TimecourseDataSet::GetClusterNoise(int nodeID)
{
  if (noise_mode == 0)
    {
      //cout << "GCN: " << 0 << endl;
      return 0;
    }
  else if (noise_mode == 1)
    {
      //cout << "GCN: Fixed = " << noiseData[0] << endl;
      return noiseData[0];
    }
  else if (noise_mode == 2)
    {
      //cout << "GCN: Node[" << nodeID << "] = g" << noiseData[nodeID] << endl;
      return noiseData[nodeID];
    }
  else
    {
      // shouldn't ever get here
      return 0;
    }
}

/* ---------------------------------------------------------------------- */

double TimecourseDataSet::GetMergedClusterNoise(vector<int> mergedclusterindex)
{
  double newnoise;	// New calculated noise

  if (noise_mode == 0)
    {
      //cout << "GMCN: " << 0 << endl;
      return 0;
    }
  else if (noise_mode == 1)
    {
      //cout << "GMCN: Fixed = " << noiseData[0] << endl;
      return noiseData[0];
    }
  else if (noise_mode == 2)
    {
      newnoise = GetClusterSEMS(mergedclusterindex);
      return newnoise;
    }
  else
    {
      // shouldn't ever get here
      return 0;
    }
}

/* ----------------------------------------------------------------------
   Standard deviation of all residuals in a cluster.
---------------------------------------------------------------------- */

double TimecourseDataSet::GetClusterSE(vector<int> itemIndex)
{
  //DECLARATIONS
  unsigned int ui;
  int i,index;
  double clusterSE;
  vector<double> timepointmeans;
  vector< vector<double> > residuals;
  
  //INITIALISE timepointmeans as a vector of zeros and residuals as a
  //nTimepoints*itemIndex.size() matrix of zeros
  timepointmeans =  vector<double>(nTimePoints, 0);
  for (ui=0; ui<itemIndex.size(); ui++)
  {
    residuals.push_back(vector<double>(nTimePoints, 0));
  }
  clusterSE = 0.0;

  //Find the cluster mean	and residuals:(yValue(t) – pointwise mean of all
  //genes in cluster(t))
  for (i=0; i<nTimePoints; i++)
  {
    for (ui=0; ui<itemIndex.size(); ui++)
    {
      index = itemIndex[ui]; //Geneindex
      timepointmeans[i] += data[index][i];
    }
    timepointmeans[i] /= itemIndex.size();
    //cout << "timepointmeans[i]" << timepointmeans[i] << endl;
    for (ui=0; ui<itemIndex.size(); ui++)
    {
      index = itemIndex[ui]; //Geneindex
      residuals[ui][i] = data[index][i]-timepointmeans[i];
      //cout << "residuals" << ui <<i << residuals[ui][i] << endl;
    }
  }

  //Calculate the standard error
  for (i=0; i<nTimePoints; i++)
  {
    for (ui=0; ui<itemIndex.size(); ui++)
    {
      clusterSE += pow(residuals[ui][i], 2);
    }
  }
  clusterSE /= (itemIndex.size()*nTimePoints - 1.0);
  //clusterSE = sqrt(clusterSE);
  return clusterSE;
}

/* ----------------------------------------------------------------------
   SEMS - standard error of the mean squared.
---------------------------------------------------------------------- */

double TimecourseDataSet::GetClusterSEMS(vector<int> itemIndex)
{
  unsigned int i;
  int index, numobs;
  double clusterSEMS;
  
  clusterSEMS = 0.0;
  for (i=0; i< itemIndex.size(); i++)
  {
    index = itemIndex[i];
    clusterSEMS += noiseData[index]; //noiseData should be sigma^2
    //cout << "clusterSEM" << clusterSEM <<endl;
  }
  //cout << "reps:" << reps << endl;
  numobs = reps*nTimePoints;
  ////cout << "itemsize:" <<itemIndex.size() << endl;
  //cout <<"the constant:" << 95.0/((itemIndex.size()*96.0) - 1.0) << endl;
  clusterSEMS *= (numobs-1)/((itemIndex.size()*numobs) - 1.0); //96 = 4 replicates * 24 timepoints
  //cout << "constant * clusterSEM" << clusterSEM <<endl;
  //clusterSEM = sqrt(clusterSEM); //only ever deal with sigma^2
  return clusterSEMS;
}

/* ---------------------------------------------------------------------- */

vector<double> TimecourseDataSet::GetDataForCluster(vector<int> itemIndex)
{
  //DECLARATIONS
  int                      i, j;
  int                      index, counter=0;
  int                      nCurrentItems=itemIndex.size();
  vector<double>           extractedData, yValues;

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
  
  //RE-CONSTRUCT THE DATA VECTOR SO IT HAS THE CORRECT ORDERING FOR THE
  //BlockCovarianceMatrix FORM WE'LL BE USING
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
  return yValues;
}

/* ----------------------------------------------------------------------
   Method to compute the determinant of a covariance function matrix.
---------------------------------------------------------------------- */

double TimecourseDataSet::ComputeLogDeterminant(double* choleskyMatrix, int nVariables)
{
  //DECLARATIONS
  int    i;
  double logDeterminant = 0;
  //COMPUTE THE LOG-DETERMINANT, GIVEN THE CHOLESKY DECOMPOSITION MATRIX
  for (i=0; i<nVariables; i++)
    logDeterminant += choleskyMatrix[i + nVariables*i];//diagonal values in 1D format
  //also need a factor of two because DET(A) = DET(L)^2
  logDeterminant *= 2;
  return(logDeterminant);
}

/* ----------------------------------------------------------------------
   Compute the log-evidence for a GP model, given a covariance function
   and data. Currently we're computing two LU decompositions; can
   we pre-compute it?
---------------------------------------------------------------------- */

double TimecourseDataSet::ComputeLogEvidence(BlockCovarianceMatrix blockMatrix,
					     vector<double> data)
{
  //DECLARATIONS
  double logEvidence;
  const double PI=3.14159265358979324;

  logEvidence  = -0.5 * blockMatrix.ComputeMatrixDeterminant();
  logEvidence -= 0.5 * nTimePoints * blockMatrix.blockSize * log(2*PI);

  //invert the BlockMatrix
  blockMatrix.InvertMatrix();
  //compute the likelihood term
  logEvidence -= 0.5*blockMatrix.ComputeLogLikelihoodProduct(data);

  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)
  //ARE THERE ANY DANGERS TO DOING THIS??
  if (logEvidence==numeric_limits<double>::infinity())
    logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}

/* ----------------------------------------------------------------------
   Compute the part of the log-evidence for a robust case with a single
   missing observation from yValues.
---------------------------------------------------------------------- */

double TimecourseDataSet::
ComputePartRobustLogEvidenceMissingSingleObservation(double logDetCovarFunctionk,
						     BlockCovarianceMatrix invblockMatrix,
						     vector<double> yValsMissingSingleObservation)
{
  //DECLARATIONS
  double logEvidence;
  const double PI=4.0*atan(1.0);
  //CONSTRUCT THE LOG-EVIDENCE USING THE BlockMatrix CLASS
  logEvidence  = -0.5 * logDetCovarFunctionk;
  logEvidence -= 0.5 * log(2*PI) * ((nTimePoints * invblockMatrix.blockSize)-1);

  //compute the likelihood term
  logEvidence -= 0.5*invblockMatrix.ComputeLogLikelihoodProductMissingSingleObservation(yValsMissingSingleObservation);

  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)
  //ARE THERE ANY DANGERS TO DOING THIS??
  //if (logEvidence==numeric_limits<double>::infinity())
  //  logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}

/* ----------------------------------------------------------------------
   Compute the partial derivative for a hyperparameter (given implicitly
   in covarianceDerivative).

   See page 114 of Rasmussen and Williams for the relevant equation.
---------------------------------------------------------------------- */

double TimecourseDataSet::
ComputeGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
		const BlockCovarianceMatrix& covarianceDerivative,
		const vector<double>& alpha)
{
  // At first glance, it would appear to be quicker to explicitly build the
  // rows as you go along; but doing it this way cannot exploit the STL
  // and so it turns out quicker to just use GetRow() as below - R.D.

  // Declarations
  size_t k;
  int j;
  double gradient=0.0;
  const int alpha_size=alpha.size();
  vector<double> vec1=vector<double>(alpha_size);
  vector<double> vec2=vector<double>(alpha_size);
  double alpha_dot_vec2, vec1_dot_vec2;
  int block_j1, block_j2, block_j1_counter, block_j2_counter;
  vector<double>::const_iterator inIt1, inIt2, ub;
  vector<double>::iterator vecIt1, vecIt2, outItEnd1, outItEnd2;
  
  // Compute the gradient
  block_j1_counter = block_j2_counter = block_j1 = block_j2 = 0;
  for (j=0; j<alpha_size; j++)
  {
    //vec1 = inverseCovarianceFunction.GetRow(j);
    //vec2 = covarianceDerivative.GetRow(j);
    // ^ we merge these two operations into one (convoluted but quick):
    ///////////////////////////////////////////////////////////////////////
    inIt1 = inverseCovarianceFunction.noiseFreeCoeff[block_j1].begin();
    vecIt1 = vec1.begin();
    inIt2 = covarianceDerivative.noiseFreeCoeff[block_j2].begin();
    vecIt2 = vec2.begin();
    ub = inverseCovarianceFunction.noiseFreeCoeff[block_j1].end();
    while (inIt1 != ub)
      {
	outItEnd1 = vecIt1 + inverseCovarianceFunction.blockSize;
	outItEnd2 = vecIt2 + covarianceDerivative.blockSize;
	fill(vecIt1, outItEnd1, *inIt1++);
	fill(vecIt2, outItEnd2, *inIt2++);
	vecIt1 = outItEnd1;
	vecIt2 = outItEnd2;
      }
    vec1[j] *= 1.0 + inverseCovarianceFunction.noisyCoeff[block_j1];
    vec2[j] *= 1.0 + covarianceDerivative.noisyCoeff[block_j2];
    // block_j1 = j / inverseCovarianceFunction.blockSize
    // block_j2 = j / covarianceDerivative.blockSize
    // ^ we can do this without doing the expensive divisions:
    if(++block_j1_counter==inverseCovarianceFunction.blockSize)
      {
	block_j1_counter=0; ++block_j1;
      }
    if(++block_j2_counter==covarianceDerivative.blockSize)
      {
	block_j2_counter=0; ++block_j2;
      }
    /////////////////////////////////////////////////////////////////////
    
    // Compute alpha.dK/dtheta and K^-1.dK/dtheta
    alpha_dot_vec2 = vec1_dot_vec2 = 0;
    for (k=0; k<alpha.size(); k++)
    {
      alpha_dot_vec2 += alpha[k] * vec2[k];
      vec1_dot_vec2 += vec1[k] * vec2[k];
    }
    
    // This element's contribution to the trace
    gradient += alpha[j] * alpha_dot_vec2 - vec1_dot_vec2;
  }

  return gradient * 0.5;
}

/* ---------------------------------------------------------------------- */

double TimecourseDataSet::
ComputeNoiseGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
		     const vector<double>& alpha,
		     double noiseSigma)
{
  //DECLARATIONS
  unsigned int i;
  double gradient=0, currentElement;
  
  //COMPUTE THE GRADIENT
  for (i=0; i<alpha.size(); i++)
  {
    currentElement = alpha[i] * alpha[i];
    currentElement -= inverseCovarianceFunction.GetElement(i,i);
    gradient += currentElement;
  }
  //was: gradient * noiseSigma;
  return gradient * 0.5;
}

/* ---------------------------------------------------------------------- */

double TimecourseDataSet::
ComputeNoiseGradientIncludingGammaHyperprior(BlockCovarianceMatrix inverseCovarianceFunction, vector<double> alpha, double noiseSigma, vector<double> gammaParams)
{
  //DECLARATIONS
  unsigned int i;
  double gradient=0, currentElement, gammaCont;
  vector<double> vec1;
  //COMPUTE THE GRADIENT
  for (i=0; i<alpha.size(); i++)
  {
    currentElement  = alpha[i] * alpha[i];
    currentElement -= inverseCovarianceFunction.GetElement(i,i);
    gradient       += currentElement;
  }
  gradient *= 0.5;
  gammaCont = gammaParams[0] -1.0; //gammaParams[0] = alpha
  gammaCont /= noiseSigma;
  gammaCont -= gammaParams[1]; //gammaParams[1] = beta
  gradient += gammaCont;
  return(gradient);
}

/* ----------------------------------------------------------------------
   Add Gaussian white noise to the input covariance function.
---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
AddNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, double noiseSigma)
{
  // Declarations
  int i;
  double sigmaSquared;
  
  // Add sigma^2 (variange) to the diagonal elements
  sigmaSquared           = pow(noiseSigma, 2);
  for (i=0; i<blockMatrix.nRank; i++)
    // normalisation because of the way we construct a block matrix
    blockMatrix.noisyCoeff[i] = sigmaSquared / blockMatrix.noiseFreeCoeff[i][i];

  return(blockMatrix);
}

/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
AddFixedNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix,
				  double noise_std_error)
{
  // Declarations
  int i;
  double sigmaSquared;

  // Add sigma^2 (variance) to the diagonal elements
  sigmaSquared           = pow(noise_std_error, 2);
  for (i=0; i<blockMatrix.nRank; i++)
    //normalisation because of the way we construct a block matrix
    blockMatrix.noisyCoeff[i] = sigmaSquared / blockMatrix.noiseFreeCoeff[i][i];
  
  return(blockMatrix);
}

/* ----------------------------------------------------------------------
   Compute a square exponential (SE) covariance function for a K matrix
   missing an observation at time k from the full K matrix.
   This has same values as the full K matrix for: noisyCoeff, noiseFreeCoeff,
   blockSize, nRank; except the timePoint corresponding to missing obs
   comes first, other timePoints are in order.
---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
CovarianceFunctionMissingSingleObservation(vector<vector<double> >KnFC,
					   vector <double> KnC,
					   int KblockSize,
					   int KnRank, 
					   int timePoint)
{
  //DECLARATIONS
  unsigned int i;
  vector<double>::iterator it; //it is an iterator for a vector of doubles
  vector<double> copyRow;
  BlockCovarianceMatrix  blockMatrix;
  //INITIALISE THE BLOCK MATRIX
  blockMatrix.nRank     = KnRank;
  blockMatrix.blockSize = KblockSize;
  //INITIALISE THE COVARIANCE FUNCTION
  blockMatrix.noiseFreeCoeff = KnFC;
  blockMatrix.noisyCoeff = KnC;
  //timepoints start at 0
  //Delete value for kth timepoint from noisyCoeff and insert it instead as the first element in the vector
  it = blockMatrix.noisyCoeff.begin() + timePoint;
  blockMatrix.noisyCoeff.erase(it);
  it = blockMatrix.noisyCoeff.begin();
  blockMatrix.noisyCoeff.insert(it, KnC[timePoint]);
  //Delete the column for the kth timepoint and insert it as the first column
  for (i=0; i<KnC.size(); i++)  //going down the rows of the noiseFreeCoeff
  {
    blockMatrix.noiseFreeCoeff[i].erase(blockMatrix.noiseFreeCoeff[i].begin()+timePoint);
    blockMatrix.noiseFreeCoeff[i].insert(blockMatrix.noiseFreeCoeff[i].begin(), KnFC[i][timePoint]);
  }
  //Swap over the rows
  copyRow = blockMatrix.noiseFreeCoeff[timePoint]; //Copy the row for time timePoint before deleting it
  blockMatrix.noiseFreeCoeff.erase(blockMatrix.noiseFreeCoeff.begin()+timePoint); //Delete the row for the kth timepoint from noiseFreeCoeff

  blockMatrix.noiseFreeCoeff.insert(blockMatrix.noiseFreeCoeff.begin(), copyRow);  //...and insert it instead as the first row

  return(blockMatrix);
}

/* ----------------------------------------------------------------------
   Compute the partial derivative w.r.t. length-scale of a square
   exponential (SE) covariance function. Note that we've hard-coded
   noise=0 here (ugly).
---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
SquareExponential_lengthDerivative(double lengthScale,
				   int blockSize,
				   double noiseFreeScale)
{
  //DECLARATIONS
  int                    i, j;
  double                 covarElement, deltaTime;
  BlockCovarianceMatrix  blockMatrix;
  //ASSERTIONS ABOUT INPUTS
  assert(lengthScale > 0);
  assert(blockSize > 0);
  //INITIALISE THE BLOCK MATRIX
  blockMatrix.nRank     = nTimePoints;
  blockMatrix.blockSize = blockSize;
  //INITIALISE THE COVARIANCE FUNCTION (do this so we can assign to 2 elements each time)
  for (i=0; i<nTimePoints; i++)
  {
    blockMatrix.noiseFreeCoeff.push_back(vector<double>(nTimePoints, 0));
    blockMatrix.noisyCoeff.push_back(0);
  }
  //COMPUTE EACH ELEMENT OF THE COVARIANCE FUNCTION
  for (i=0; i<nTimePoints; i++)
  {
    for (j=i; j<nTimePoints; j++)
    {
      deltaTime     = fabs(timePoints[i] - timePoints[j]);
      covarElement  = pow(deltaTime, 2);
      covarElement /= 2*pow(lengthScale, 2);
      covarElement *= -1;
      covarElement  = exp(covarElement);
      covarElement *= pow(deltaTime, 2);
      covarElement /= pow(lengthScale, 3);
      covarElement *= noiseFreeScale;
      //and store in 2 elements (covariance is symmetric)
      //this duplicates effort for i==j; not a big deal, computationally :-)
      blockMatrix.noiseFreeCoeff[i][j] = covarElement;
      blockMatrix.noiseFreeCoeff[j][i] = covarElement;
    }
  }
  return(blockMatrix);
}

/* ----------------------------------------------------------------------
   Optimise the noise hyperparameter given the measurement error
   (SEM calculated from replicates).
---------------------------------------------------------------------- */

vector<double> TimecourseDataSet::OptimiseGammaParams(double clusterSEMS)
{
  //DECLARATIONS
  double gap, alpha, beta, c, d, maxx, diff, bestAlpha, bestBeta;
  vector<double> bestParams;
  //INITIALISE
  gap = 100.0;
  //iterate over beta to find alpha
  for (beta =2.0; beta < 150.0; beta = beta+0.5)
  {
    alpha = GetAlpha(beta, clusterSEMS);
    c = GammaDistribution(clusterSEMS, alpha, beta);
    maxx = (alpha -1.0)/beta;
    d = GammaDistribution(maxx, alpha, beta)/100.0;
    diff = fabs(c - d);
    //cout << "diff:" << diff<< endl;
    if (diff < gap)
    {
      gap = diff;
      bestAlpha = alpha;
      bestBeta = beta;
    }
    //cout << "bestAlpha:" << bestAlpha << "bestBeta:" << bestBeta << endl;
  }
  bestParams.push_back(bestAlpha);
  bestParams.push_back(bestBeta);
  return(bestParams);

}

/* ----------------------------------------------------------------------
  Get the alphap parameter for Gamma distribution, given the beta parameter
  and sigmaSEM.
---------------------------------------------------------------------- */

double TimecourseDataSet::GetAlpha(double beta, double clusterSEMS)
{
  //DECLARATIONS
  double alpha, c, d, e;
  //Calculate Alpha from beta and sigma
  c = beta*(clusterSEMS - 1.0);
  d = log(clusterSEMS);
  //cout << "d:" << d << endl;
  e = c/d;
  alpha = e +1.0;
  //cout << "beta:" << beta << "alpha:" << alpha << endl;
  return(alpha);
}

/* ----------------------------------------------------------------------
   Get the density for the Gamma Distribution at the input value, for params
   alpha and beta.
---------------------------------------------------------------------- */

double TimecourseDataSet::GammaDistribution(double input, double alpha, double beta)
{
  //DECLARATIONS
  double c, d, e, f, gammadist;
  //Calculate Gamma distribution
  c = pow(beta, alpha);
  d = tgamma(alpha);
  e = pow(input, (alpha-1.0));
  f = exp(-(beta*input));
  gammadist = c*e*f;
  gammadist /= d;
  return(gammadist);
}

/* ---------------------------------------------------------------------- */

int TimecourseDataSet::GetNoiseMode()
{
  return(noise_mode);
}

/* ---------------------------------------------------------------------- */

int TimecourseDataSet::GetRobustMode()
{
  return(robust_mode);
}

/* ----------------------------------------------------------------------
   Compute a noise-less square exponential (SE) covariance function
   which misses all observation from one timepoint.
---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
SquareExponentialCovarianceFunctionMissOneTimePoint(double lengthScale,
						    int blockSize,
						    double noiseFreeScale,
						    int timePoint)
{
  //DECLARATIONS
  int                    i, j, k, l; //i,j are timepoints and j,k are matrix positions
  double                 covarElement;
  BlockCovarianceMatrix  blockMatrix;
  //INITIALISE THE BLOCK MATRIX
  blockMatrix.nRank     = nTimePoints-1;
  blockMatrix.blockSize = blockSize;
  //INITIALISE THE COVARIANCE FUNCTION (do this so we can assign to 2 elements each time)
  for (i=0; i<nTimePoints-1; i++)
  {
    blockMatrix.noiseFreeCoeff.push_back(vector<double>(nTimePoints-1, 0));
    blockMatrix.noisyCoeff.push_back(0);
  }
  k=0;
  l=0;
  //COMPUTE EACH ELEMENT OF THE COVARIANCE FUNCTION
  for (i=0; i<nTimePoints; i++)
  {
    if (i != timePoint)
    {
      for (j=i; j<nTimePoints; j++)
      {
        if (j != timePoint)
        {
          covarElement = fabs(timePoints[i] - timePoints[j]);
          covarElement = pow(covarElement, 2);
          covarElement /= 2*pow(lengthScale, 2);
          covarElement *= -1;
          covarElement  = exp(covarElement);
          covarElement *= noiseFreeScale;
          //and store in 2 elements (covariance is symmetric)
          //this duplicates effort for i==j; not a big deal, computationally :-)
          blockMatrix.noiseFreeCoeff[k][l] = covarElement;
          blockMatrix.noiseFreeCoeff[l][k] = covarElement;
          l++;
        }
      }
      k++;
      l=k;
    }
  }
  return(blockMatrix);
}

/* ----------------------------------------------------------------------
   Perform a line seach. This is used by DFPMaximise() below.
   Adapted from Numerical Recipes in C++.
---------------------------------------------------------------------- */

void TimecourseDataSet::
LineSearch(vector<double>& xold, const double fold, vector<double>& g,
	   vector<double>& p, vector<double>& x, double& f, const double stpmax,
	   bool& check, const int blockSize, const vector<double>& yValues)
{
  const double ALF=1.0e-3, TOLX=numeric_limits<double>::epsilon();
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double rhs1,rhs2,slope,sum,temp,test,tmplam;

  int n=xold.size();
  check=false;
  sum=0.0;
  for(i=0; i<n; i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    {
      temp=stpmax/sum;
      for(i=0; i<n; i++) p[i] *= temp;
    }
  slope=0.0;
  for(i=0; i<n; i++) slope += g[i]*p[i];
  //if(slope >= 0.0) cout << "Roundoff problem in line_search: " << slope << endl;
  test=0.0;
  for(i=0; i<n; i++)
    {
      temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
      if(temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for(;;)
    {
      for(i=0; i<n; i++) x[i]=xold[i] + alam*p[i];
      ImposeConstraintsOnHyperparameters(x);
      f=ComputeLogEvidenceFromHyperparameters(yValues,blockSize,x);
      f=-f;
      if(alam < alamin)
	{
	  for(i=0; i<n; i++) x[i]=xold[i];
	  check=true;
	  return;
	}else if(f <= fold + ALF*alam*slope) return;
      else
	{
	  if(alam == 1.0)
	    tmplam = -slope/(2.0*(f-fold*slope));
	  else
	    {
	      rhs1=f-fold-alam*slope;
	      rhs2=f2-fold-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if(a == 0.0) tmplam = -slope/(2.0*b);
	      else
		{
		  disc=b*b-3.0*a*slope;
		  if(disc < 0.0) tmplam=0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		  else tmplam=-slope/(b+sqrt(disc));
		}
	      if(tmplam>0.5*alam)
		tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2=f;
      alam=MAX(tmplam,0.1*alam);
    }
}

/* ----------------------------------------------------------------------
   Uses the Broyden-Fletcher-Goldfarb-Shanno (BFGS) variant of the
   Davidon-Fletcher-Powell (DFP) optimisation method (quasi-Newton) to
   perform function maximisation.

   This function can only be called from a class that inherits this class,
   and which must define the functions:
      - ComputeLogEvidenceAndGradientsFromHyperparameters()
      - ComputeGradientsFromHyperparameters()
      - ComputeLogEvidenceFromHyperparameters()
      - ImposeConstraintsOnHyperparameters()

   Adapted from Numerical Recipes in C++.
------------------------------------------------------------------------- */

void TimecourseDataSet::
DFPMaximise(vector<double>& p, const vector<int>& fix, const double gtol, double& fret, const int blockSize, const vector<double>& yValues)
{
  const int ITMAX=100;
  const double EPS=numeric_limits<double>::epsilon();
  const double TOLX=4*EPS, STPMX=100.0;
  bool check, fixedparam=false;
  int i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

  int n=p.size();
  vector<double> dg(n), g(n), hdg(n), pnew(n), xi(n);
  vector<vector<double> > hessin=vector<vector<double> >(n,vector<double>(n,0.0));

  for(i=0; i<n; i++) if(fix[i]) {fixedparam=true; break;}
  
  ComputeLogEvidenceAndGradientsFromHyperparameters(yValues,blockSize,
						    p,fp,g);
  fp=-fp;
  for(i=0; i<n; i++) g[i]=-g[i];
  if(fixedparam) for(i=0; i<n; i++) if(fix[i]) g[i]=0;
  for(i=0; i<n; i++)
    {
      hessin[i][i]=1.0; // identity
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }
  stpmax=STPMX*MAX(sqrt(sum),double(n));
  for(its=0; its<ITMAX; its++)
    {
      LineSearch(p,fp,g,xi,pnew,fret,stpmax,check,blockSize,yValues);

      fp=fret;
      for(i=0; i<n; i++)
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}
      test=0.0;
      for(i=0; i<n; i++)
	{
	  temp=fabs(xi[i])/MAX(fabs(p[i]),1.0);
	  if(temp>test)test=temp;
	}
      if(test < TOLX) return;
      for(i=0; i<n; i++) dg[i]=g[i];
      ComputeGradientsFromHyperparameters(yValues,blockSize,
					  p,g);
      for(i=0; i<n; i++) g[i]=-g[i];
      if(fixedparam) for(i=0; i<n; i++) if(fix[i]) g[i]=0;
      test=0.0;
      den=MAX(fret,1.0);
      for(i=0; i<n; i++)
	{
	  temp=fabs(g[i])*MAX(fabs(p[i]),1.0)/den;
	  if(temp>test)test=temp;
	}
      if(test<gtol) return;
      for(i=0; i<n; i++) dg[i]=g[i]-dg[i];
      for(i=0; i<n; i++)
	{
	  hdg[i]=0.0;
	  for(j=0; j<n; j++) hdg[i] += hessin[i][j]*dg[j];
	}
      fac=fae=sumdg=sumxi=0.0;
      for(i=0; i<n; i++)
	{
	  fac += dg[i]*xi[i];
	  fae += dg[i]*hdg[i];
	  sumdg += dg[i]*dg[i];
	  sumxi += xi[i]*xi[i];
	}
      if(fac > sqrt(EPS*sumdg*sumxi))
	{
	  fac=1.0/fac;
	  fad=1.0/fae;
	  for(i=0; i<n; i++) dg[i]=fac*xi[i]-fad*hdg[i];
	  for(i=0; i<n; i++)
	    {
	      for(j=i; j<n; j++)
		{
		  hessin[i][j] += fac*xi[i]*xi[j]
		    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
		  hessin[j][i]=hessin[i][j];
		}
	    }
	}
      for(i=0; i<n; i++)
	{
	  xi[i]=0.0;
	  for(j=0; j<n; j++) xi[i] -= hessin[i][j]*g[j];
	}
    }
  // If we get here, then we reached max number of iterations before convergence
}

/* ---------------------------------------------------------------------- */

BlockCovarianceMatrix TimecourseDataSet::
SquareExponentialCovarianceFunctionMissingSingleObservation(vector<vector<double> > KnFC,
							    vector<double> KnC,
							    int KblockSize,
							    int KnRank,
							    int timePoint)
{
  BlockCovarianceMatrix bcm;
  cout << "SquareExponentialCovarianceFunctionMissingSingleObservation() has not been implemented in TimecourseDataSet!" << endl;
  return bcm;
}

/* ---------------------------------------------------------------------- */

double TimecourseDataSet::GetMLIINoise(vector<int> itemIndex)
{
  cout << "GetMLIINoise() has not been implemented in TimecourseDataSet!" << endl;
  return -1;
}

/* ---------------------------------------------------------------------- */
