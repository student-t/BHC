#include "TimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"
#include <limits>
#include <numeric>


//CONSTRUCTOR
TimecourseDataSet::TimecourseDataSet()
{}



//CONSTRUCTOR - READ IN DATA FROM FILE
TimecourseDataSet::TimecourseDataSet(string dataFile)
{
  //READ IN THE DATA FROM FILE
  ReadInData(dataFile);
}



//CONSTRUCTOR - READ IN DATA FROM INPUT VECTOR
//need this in order to link to R
TimecourseDataSet::TimecourseDataSet(const vector<vector<double> >& inputData)
{
  //DECLARATIONS
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



//READ DATA IN FROM A FILE
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



void TimecourseDataSet::ReadInTimePoints(vector<double> timePoints_copy)
{
  //DECLARATIONS
  int	i;
  
  for (i=0; i<nTimePoints; i++)
  {
    timePoints.push_back(timePoints_copy[i]);
  }
}



//READ DATA IN FROM A FILE
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



//READ DATA IN FROM A FILE WHEN RUNNING FROM R
void TimecourseDataSet::ReadInNoise(vector<double> noise_copy)
{
  noiseData = noise_copy;
}



void TimecourseDataSet::SetNoiseMode(int mode)
{
  noise_mode = mode;
  ///cout << "Noise mode = " << noise_mode << endl;
}



void TimecourseDataSet::SetReps(int num_reps)
{
  reps = num_reps;
  ///cout << "Num reps = " << reps << endl;
}



void TimecourseDataSet::SetDataType(string type)
{
  dataType = type;
  ///cout << "dataType:" << dataType << endl;
}



void TimecourseDataSet::SetRobustMode(int mode)  //also sets the dataRange if Robust = 1
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



//double TimecourseDataSet::GetMergedClusterNoise(double node1Noise, double node2Noise){
//	double newnoise;	// New calculated noise

//	if (noise_mode == 0){
//cout << "GMCN: " << 0 << endl;
//		return 0;
//	}
//	else if (noise_mode == 1){
//cout << "GMCN: Fixed = " << noiseData[0] << endl;
//		return noiseData[0];
//	}
//	else if (noise_mode == 2){
//		newnoise = pow((95.0/191.0)*(pow(node1Noise,2) + pow(node2Noise, 2)),0.5); //the constant is not always e.g. 95/191
//cout << "GMCN: Propagating = " << newnoise << endl;
//		return newnoise;
//	}
//}



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



double TimecourseDataSet::GetClusterSE(vector<int> itemIndex)  //standard deviation of all residuals in a cluster
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



//SEMS Standard Error of the Mean Squared
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



/*
double TimecourseDataSet::GetMLIINoise(vector<int> itemIndex){
  //DECLARATIONS
  int                      i, j;
  int                      index, counter=0;
  int                      nCurrentItems=itemIndex.size();
  vector<double>           extractedData, yValues;
  double				   fittedNoise;
  //EXTRACT THE DATA POINTS FOR THIS CURRENT CLUSTER
  //data vector is nDataItems*nTimePoints
  for (i=0; i<nCurrentItems; i++){
    for (j=0; j<nTimePoints; j++){
      index = itemIndex[i];
      extractedData.push_back(data[index][j]);//store the corresponding data value
    }
  }
  //RE-CONSTRUCT THE DATA VECTOR SO IT HAS THE CORRECT ORDERING FOR THE BlockCovarianceMatrix FORM WE'LL BE USING
  yValues = extractedData;
  counter = 0;
  for (i=0; i<nCurrentItems; i++){
    for (j=0; j<nTimePoints; j++){
      index          = j*nCurrentItems + i;
      yValues[index] = extractedData[counter];
      counter++;
    }
  }
  //Optimise the hyperparameters and return the fitted noise
  fittedNoise = CalculateFittedNoiseHyperparameter(yValues);
  return(fittedNoise);
}
*/



//METHOD TO COMPUTE THE ROBUST LOGEVIDENCE MISSING OUT A TIMEPOINT EACH TIME
/*
double TimecourseDataSet::ComputeRobustLogEvidence(vector<double> extractedData, vector<double> yValues, int nCurrentItems, double* p_lengthScale, double* p_NoiseFree, double* p_Noise){
	//DECLARATIONS
	int		k,i,j, index, counter;
	double  lft, lst, cft, cst, denom, diff, a, b, B, logEvidence, min, max;
	long double ft, evidence, pst, sst;
	vector<double>	yValuesCopy;
	vector<double>::iterator itmin, itmax;
	BlockCovarianceMatrix    covarFunction;
	//INITIALISE
	b = 0.01; //prob a value is an outlier
	a = 1-b; //prob a value is a regular value
	//B = 0.83; //the constant value of the outlier PDF, imagine stdev is 3 times greater for outlier, 1/0.4*3=0.83
	sst = 0;
	//if (nCurrentItems = 1){ //could still miss out a timepoint for single genes
	//	return ComputeMaximisedLogEvidence(yValues, p_lengthScale, p_NoiseFree, p_Noise);
	//}
	//Set B as 1/range data
	itmin = min_element( yValues.begin(), yValues.end() );
	itmax = max_element( yValues.begin(), yValues.end() );
	min = *itmin;
	max = *itmax;
	diff = max-min;
	B = 1.0/diff;
	//cout << "min:" << min << "max" << max << "B:" << B << endl;
	//COMPUTE FIRST TERM
	lft = ComputeMaximisedLogEvidence(yValues, p_lengthScale, p_NoiseFree, p_Noise);
	ft = exp((long double)lft);
	//COMPUTE SECOND TERM
	//miss out a timepoint each time
	for (k=0; k<nTimePoints;k++){
		yValuesCopy = yValues;
		for (i=nTimePoints*(nCurrentItems-1) + k; i >=  0; i=i-nTimePoints){
			yValuesCopy.erase (yValuesCopy.begin()+i);
		}
		covarFunction = SquareExponentialCovarianceFunctionMissOneTimePoint(*p_lengthScale, nCurrentItems, *p_NoiseFree, k); //misses out timepoint k
		covarFunction = AddNoiseToCovarianceFunction(covarFunction, *p_Noise);
		lst = ComputeLogEvidence(covarFunction, yValuesCopy);
		pst = exp((long double)lst);
		sst += pst;
	}
	//ADD TERMS AND TAKE LOG
	cft = pow(a, nCurrentItems);
	cst = pow(a, (nCurrentItems-1));
	cst *=pow(b, 1);
	//denom = cft+cst;
	//cst /= denom;
	//cft /= denom;
	cst *= B;
	evidence = cft*ft + cst*sst;
	logEvidence = log(evidence);
	return(logEvidence);
}
*/



//METHOD TO COMPUTE THE DETERMINANT OF A COVARIANCE FUNCTION MATRIX
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



// COMPUTE THE LOG-EVIDENCE FOR A GP MODEL, GIVEN A COVARIANCE FUNCTION AND DATA
//currently, we're computing two LU decompositions :-(
//can we pre-compute it, and therefore speed up the det and inv functions??
double TimecourseDataSet::ComputeLogEvidence(BlockCovarianceMatrix blockMatrix, vector<double> data)
{
  //DECLARATIONS
  double logEvidence;
  const double PI=4.0*atan(1.0); // ugly
  //CONSTRUCT THE LOG-EVIDENCE USING THE BlockMatrix CLASS
  //cout << "blockMatrix.noiseFreeCoeff[0]:" << blockMatrix.noiseFreeCoeff[0][0] << " " << blockMatrix.noiseFreeCoeff[0][1]<< " " << blockMatrix.noiseFreeCoeff[0][2] << " " << blockMatrix.noiseFreeCoeff[0][3] << endl;
  //cout << "blockMatrix.noiseFreeCoeff[1]:" << blockMatrix.noiseFreeCoeff[1][0] << " " << blockMatrix.noiseFreeCoeff[1][1]<< " " << blockMatrix.noiseFreeCoeff[1][2] << " " << blockMatrix.noiseFreeCoeff[1][3]<< endl;
  //cout << "blockMatrix.noiseFreeCoeff[2]:" << blockMatrix.noiseFreeCoeff[2][0] << " " << blockMatrix.noiseFreeCoeff[2][1]<< " " << blockMatrix.noiseFreeCoeff[2][2] << " " << blockMatrix.noiseFreeCoeff[2][3]<< endl;
  //cout << "blockMatrix.noiseFreeCoeff[3]:" << blockMatrix.noiseFreeCoeff[3][0] << " " << blockMatrix.noiseFreeCoeff[3][1]<< " " << blockMatrix.noiseFreeCoeff[3][2] << " " << blockMatrix.noiseFreeCoeff[3][3]<< endl;
  //cout << "blockMatrix.noisyCoeff:" << blockMatrix.noisyCoeff[0] << " " << blockMatrix.noisyCoeff[1] << " " << blockMatrix.noisyCoeff[2] << " " << blockMatrix.noisyCoeff[3] << endl;

  logEvidence  = -0.5 * blockMatrix.ComputeMatrixDeterminant();
  //cout << "part1:" << logEvidence << endl;
  logEvidence -= 0.5 * nTimePoints * blockMatrix.blockSize * log(2*PI);
  //cout << "part2:" << -0.5 * nTimePoints * blockMatrix.blockSize * log(2*PI) << endl;
  //invert the BlockMatrix
  blockMatrix.InvertMatrix();
  //compute the likelihood term
  logEvidence -= 0.5*blockMatrix.ComputeLogLikelihoodProduct(data);
  //cout << "part3:" << -0.5*blockMatrix.ComputeLogLikelihoodProduct(data) << endl;

  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)
  //ARE THERE ANY DANGERS TO DOING THIS??
  if (logEvidence==numeric_limits<double>::infinity())
    logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}



// COMPUTE THE PART OF THE LOG-EVIDENCE FOR A ROBUST CASE WITH A SINGLE MISSING OBS FROM Y VALUES
double TimecourseDataSet::ComputePartRobustLogEvidenceMissingSingleObservation(double logDetCovarFunctionk, BlockCovarianceMatrix invblockMatrix, vector<double> yValsMissingSingleObservation)
{
  //DECLARATIONS
  double logEvidence;
  const double PI=4.0*atan(1.0);
  //CONSTRUCT THE LOG-EVIDENCE USING THE BlockMatrix CLASS
  logEvidence  = -0.5 * logDetCovarFunctionk;
  //cout << "partdet:" << logEvidence << endl;
  logEvidence -= 0.5 * log(2*PI) * ((nTimePoints * invblockMatrix.blockSize)-1);
  //cout << "partconst:" << -(0.5 * log(2*PI) * ((nTimePoints * invblockMatrix.blockSize)-1)) << endl;
  //compute the likelihood term
  logEvidence -= 0.5*invblockMatrix.ComputeLogLikelihoodProductMissingSingleObservation(yValsMissingSingleObservation);
  //cout << "partinv:" << -0.5*invblockMatrix.ComputeLogLikelihoodProductMissingSingleObservation(yValsMissingSingleObservation);
  //cout << "y.size:" << yValsMissingSingleObservation.size() << " " << yValsMissingSingleObservation[0] << " " << yValsMissingSingleObservation[1]<< " "  << yValsMissingSingleObservation[2]<< " "  << yValsMissingSingleObservation[3]<< " "  << yValsMissingSingleObservation[4];
  //cout << "invMatrixk.noiseFreeCoeff:" << invblockMatrix.noiseFreeCoeff[0][0] << " " << invblockMatrix.noiseFreeCoeff[0][1] << " " << invblockMatrix.noiseFreeCoeff[1][0] << " " << invblockMatrix.noiseFreeCoeff[1][1];
  //cout << "invMatrixk.noisyCoeff:" <<  invblockMatrix.noisyCoeff[0] << " " << invblockMatrix.noisyCoeff[1] << endl;
  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)
  //ARE THERE ANY DANGERS TO DOING THIS??
  //if (logEvidence==numeric_limits<double>::infinity())
  //  logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}



// OPTIMISE THE HYPERPARAMETERS FOR THE GP MODEL, RETURNING THE OPTIMAL LOG-EVIDENCE
//xValues should contain the distinct time points (*not* duplicates)
//This method uses a simple implementation of a gradient ascent method.
//Hyperparameters are:
//   - length-scale (of the SE covariance function)
//   - the noiseFreeScale (ie. amplitude of the SE covariance function term)
//   - noiseSigma
//


//double TimecourseDataSet::CalculateFittedNoiseHyperparameter(vector<double> yValues){
// cout << "CFNHP: This method should be replaced by one for the subclass being used" << endl;
// }



// Compute the partial derivative for a hyperparameter (which is given
// implicitly in covarianceDerivative).
//
// See page 114 of Rasmussen and Williams.
double TimecourseDataSet::ComputeGradient(const BlockCovarianceMatrix& inverseCovarianceFunction, const BlockCovarianceMatrix& covarianceDerivative, const vector<double>& alpha)
{
  // I have found that using GetRow(), as below, is the quickest
  // way of doing this job, since explicitly building the rows whilst
  // performing the below operation cannot make use of the STL (whereas
  // GetRow() can).
  
  //DECLARATIONS
  unsigned int k;
  int j;
  double gradient=0.0;
  const int alpha_size=alpha.size();
  vector<double> vec1=vector<double>(alpha_size);
  vector<double> vec2=vector<double>(alpha_size);
  double alpha_dot_vec2, vec1_dot_vec2;
  int block_j1, block_j2, block_j1_counter, block_j2_counter;
  vector<double>::const_iterator inIt1, inIt2, ub;
  vector<double>::iterator vecIt1, vecIt2, outItEnd1, outItEnd2;
  
  //COMPUTE THE GRADIENT
  block_j1_counter = block_j2_counter = block_j1 = block_j2 = 0;
  for (j=0; j<alpha_size; j++)
  {
    //vec1 = inverseCovarianceFunction.GetRow(j);
    //vec2 = covarianceDerivative.GetRow(j);
    // ^ we merge these two operations into one (quite convoluted):
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
    // Book-keeping to avoid the following divisions:
    //  block_j1 = j / inverseCovarianceFunction.blockSize
    //  block_j2 = j / covarianceDerivative.blockSize
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



//COMPUTE THE PARTIAL GRADIENT FOR A GIVEN HYPERPARAMETER
// This function is pretty negligible at the moment, but might
// benefit later from using the STL.
double TimecourseDataSet::ComputeNoiseGradient(const BlockCovarianceMatrix& inverseCovarianceFunction, const vector<double>& alpha, double noiseSigma)
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



//COMPUTE THE PARTIAL GRADIENT FOR A GIVEN HYPERPARAMETER
double TimecourseDataSet::ComputeNoiseGradientIncludingGammaHyperprior(BlockCovarianceMatrix inverseCovarianceFunction, vector<double> alpha, double noiseSigma, vector<double> gammaParams)
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



//ADD GAUSSIAN WHITE NOISE TO THE INPUT COVARIANCE FUNCTION
BlockCovarianceMatrix TimecourseDataSet::AddNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, double noiseSigma)
{
  //DECLARATIONS
  int i;
  double sigmaSquared;
  //ADD SIGMA^2 (VARIANCE) TO THE DIAGONAL ELEMENTS
  sigmaSquared           = pow(noiseSigma, 2);
  for (i=0; i<blockMatrix.nRank; i++)
    //cout << "nRank" << blockMatrix.nRank << endl; why does this break it??
    blockMatrix.noisyCoeff[i] = sigmaSquared / blockMatrix.noiseFreeCoeff[i][i];//normalisation because of the way we construct a block matrix
  //RETURN THE NOISY COVARIANCE FUNCTION
  return(blockMatrix);
}



//ADD GAUSSIAN WHITE NOISE TO THE INPUT COVARIANCE FUNCTION
BlockCovarianceMatrix TimecourseDataSet::AddFixedNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, double noise_std_error)
{
  //DECLARATIONS
  int i;
  double sigmaSquared;
  //ADD SIGMA^2 (VARIANCE) TO THE DIAGONAL ELEMENTS
  sigmaSquared           = pow(noise_std_error, 2);
  for (i=0; i<blockMatrix.nRank; i++)
    //cout << "nRank" << blockMatrix.nRank << endl; why does this break it??
    blockMatrix.noisyCoeff[i] = sigmaSquared / blockMatrix.noiseFreeCoeff[i][i];//normalisation because of the way we construct a block matrix
  //RETURN THE NOISY COVARIANCE FUNCTION
  return(blockMatrix);
}



// COMPUTE A SQUARE EXPONENTIAL (SE) COVARIANCE FUNCTION FOR A K MATRIX MISSING OBSERVATION AT TIME k, FROM FULL K MATRIX
// THIS HAS SAME VALUES AS THE FULL K MATRIX FOR: NOISYCOEFF NOISEFREECOEFF, BLOCKSIZE, NRANK; EXCEPT THE TIMEPOINT CORRES. TO MISSING OBS COMES FIRST, OTHER TIMEPOINTS ARE IN ORDER
BlockCovarianceMatrix TimecourseDataSet::CovarianceFunctionMissingSingleObservation(vector< vector<double> > KnFC, vector <double> KnC, int KblockSize, int KnRank, int timePoint)
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



// COMPUTE THE PARTIAL DERIVATIVE WRT LENGTH-SCALE OF A SQUARE EXPONENTIAL (SE) COVARIANCE FUNCTION
//NOTE THAT WE'VE HARD-WIRED NOISE=0 HERE (THIS IS A BIT UGLY)
BlockCovarianceMatrix TimecourseDataSet::SquareExponential_lengthDerivative(double lengthScale, int blockSize, double noiseFreeScale)
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



//OPTIMISE THE NOISE HYPERPARAMETER GIVEN THE MEASUREMENT ERROR (SEM CALCULATED FROM REPLICATES)
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



//Get the alpha parameter for Gamma distribution, given the beta parameter and sigmaSEM
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



//Get the density for the Gamma Distribution at the input value, for params alpha and beta
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



int TimecourseDataSet::GetNoiseMode()
{
  return(noise_mode);
}
int TimecourseDataSet::GetRobustMode()
{
  return(robust_mode);
}



// COMPUTE A NOISE-LESS SQUARE EXPONENTIAL (SE) COVARIANCE FUNCTION, WHICH MISSES ALL OBSERVATIONS FROM ONE TIMEPOINT
BlockCovarianceMatrix TimecourseDataSet::SquareExponentialCovarianceFunctionMissOneTimePoint(double lengthScale, int blockSize, double noiseFreeScale, int timePoint)
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



//COMPARE LOG-LIKELIHOOD CALCULATIONS

/*
if (( newLogLike - oldLogLike) > 1e-6){
  cout << "*******************" << endl;
  cout  << "BlockSize: " << blockSize << "\tNew logLike: " << newLogLike << "\tOld logLike: " << oldLogLike
  << "\tRatio: " << newLogLike / oldLogLike << endl;
  //	<< "\tDifference: " <<  newLogLike - oldLogLike   << "\tRatio: " << newLogLike / oldLogLike << endl;
}
*/

/*
cout << "-------------(block)" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    dum = blockMatrix.noiseFreeCoeff[i][j];
    if (i==j) dum *= blockMatrix.noisyCoeff[i] + 1;

    cout << dum  << "\t " ;
  }
  cout << endl;
}

cout << "-------------(IT++ inversion)" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    cout <<  workingMatrix.get(i,j) << "\t " ;
  }
  cout << endl;
}
cout << "-------------(difference)" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    dum = blockMatrix.noiseFreeCoeff[i][j];
    if (i==j) dum *= blockMatrix.noisyCoeff[i] + 1;
    cout <<  dum - workingMatrix.get(i,j) << "\t " ;
  }
  cout << endl;
}
cout << "-------------" << endl;
*/


// ***ALSO WANT TO CHECK THE DATA VECTOR FOR nData==1***
//for (i=0; i<data.size(); i++)
//cout << data[i] << "\t" << newData[i] << "\t" << dataVec.get(i) << endl;


//  cout << "NewLogEv: " << newLogEv << "\tOldLogEv: " << logEvidence << endl;








//DEBUG- THIS ONLY WORKS FOR blockSize==1
/*
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){

    blockMatrix.noiseFreeCoeff[i][j]=1;
    if (i!=j) {
workingMatrix.set(i,j,0);
blockMatrix.noiseFreeCoeff[i][j]=0;
    }

  }
}
*/


//COMPARE LOG-DETERMINANT CALCULATIONS
//    cout << "BlockSize: " << blockSize << endl;


/*
if (( blockMatrix.ComputeMatrixLogDeterminant() -  log(det(workingMatrix)))>1e-6){
  cout << "*******************" << endl;
  cout  << "BlockSize: " << blockSize << "\tNew log(det): " << blockMatrix.ComputeMatrixLogDeterminant() << "\tOld log(det): " << log(det(workingMatrix))
  << "\tDifference: " <<  blockMatrix.ComputeMatrixLogDeterminant() -  log(det(workingMatrix))   << endl;
}
*/

/*
cout << "-------------" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    cout <<  workingMatrix.get(i,j) << " " ;
  }
  cout << endl;
}
cout << "-------------" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    cout <<  blockMatrix.noiseFreeCoeff[i][j] << " " ;
  }
  cout << endl;
}
cout << "-------------" << endl;
for (i=0; i<nTimePoints; i++){
  for (j=0; j<nTimePoints; j++){
    if (i==j)
cout <<  blockMatrix.noisyCoeff[i] << " " ;
    else
cout << 0 << " ";
  }
  cout << endl;
}
cout << "-------------" << endl;
*/



/*
// COMPUTE THE LOG-EVIDENCE FOR A GP MODEL, GIVEN A COVARIANCE FUNCTION AND DATA
//currently, we're computing two LU decompositions :-(
//can we pre-compute it, and therefore speed up the det and inv functions??
double TimecourseDataSet::ComputeLogEvidence(vector< vector<double> > covarianceFunction, vector<double> data){
  //DECLARATIONS
  int          i, j;
  double       logEvidence;
  int          nVariables=covarianceFunction.size(), nData=data.size();
  mat          workingMatrix(nVariables, nVariables);
  vec          workingVector(nData), dataVec(nData);
  const double PI=4.0*atan(1.0);

  BlockCovarianceMatrix    blockMatrix;
  vector< vector<double> > noiseFreeCoeff;//'a' values that define the noiseless component of the blocks
  vector<double>           noisyCoeff;    //'b' values that define the noisy, diagonal component
  int                      blockSize;
  double                   newLogLike, oldLogLike, dum, newLogEv;
  vector<double>           newData;
  int                      counter=0, index;


  //PUT THE DATA IN THE CORRECT FOR TO USE THE IT++ LIBRARY
  for (i=0; i<nData; i++)
    dataVec.set(i, data[i]);
  //PUT THE COVARIANCE FUNCTION IN THE CORRECT FORM TO USE THE IT++ LIBRARY
  for (i=0; i<nVariables; i++){
    for (j=0; j<nVariables; j++){
      workingMatrix.set(i, j, covarianceFunction[i][j]);
    }
  }
  //CONSTRUCT a BlockCovarianceMatrix CLASS VERSION OF THE COVARIANCE MATRIX
  blockSize   = nData / nTimePoints;
  blockMatrix = BlockCovarianceMatrix(workingMatrix, nTimePoints, blockSize);



  //LOGEV TERM: DETERMINANT
  logEvidence   = -0.5 * log(det(workingMatrix));
  //LOGEV TERM: CONSTANT
  logEvidence  -= 0.5 * nVariables * log(2*PI);



  //RE-CONSTRUCT THE DATA VECTOR SO IT HAS THE CORRECT ORDERING FOR THE BlockMatrix
  newData = data;
  counter = 0;
  for (i=0; i<blockSize; i++){
    for (j=0; j<nTimePoints; j++){
      index      = j*blockSize + i;

      //cout << i << "\t"<< counter << "\t" << index << endl;

      newData[index] = data[counter];
      counter++;
    }
  }



  //LOGEV TERM: MATRIX INVERSE, DATA
  workingMatrix = inv(workingMatrix);
  workingVector = workingMatrix * dataVec;//this uses a method from IT++
  logEvidence  -= 0.5 * (dataVec * workingVector);//this is a dot product



  //CONSTRUCT THE LOG-EVIDENCE USING THE BlockMatrix CLASS
  newLogEv  = -0.5 * blockMatrix.ComputeMatrixLogDeterminant();
  newLogEv -= 0.5 * nVariables * log(2*PI);


  //INVERT THE BlockMatrix
  blockMatrix.InvertMatrix();




  //COMPUTE THE LOG-LIKE USING OUR NEW METHOD
  oldLogLike = dataVec * workingVector;
  newLogLike = blockMatrix.ComputeLogLikelihoodProduct(newData);


  newLogEv -= 0.5*blockMatrix.ComputeLogLikelihoodProduct(newData);

  logEvidence = newLogEv;

  //IT MIGHT BE SENSIBLE TO FORBID logEv=inf HERE (RESULT OF A SINGULAR MATRIX)
  //ARE THERE ANY DANGERS TO DOING THIS??
  if (logEvidence==numeric_limits<double>::infinity())
    logEvidence = -numeric_limits<double>::infinity();//-ve inf gives us ev=0, which will be rejected

  return(logEvidence);
}
*/


/////////////////////



// Perform a line search
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



// This method does the BFGS/DFP minimisation.
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
  // If we get here, then reached max number of iterations before convergence
}

