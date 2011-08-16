#ifndef TIMECOURSEDATASET_H
#define TIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>
#include "DataSet.h"
#include "BlockCovarianceMatrix.h"

// This class inherits the superclass DataSet.
// NOTE: in TimecourseDataSet, nFeatures is deprecated (in favour of nTimePoints)
class TimecourseDataSet : public DataSet
{
 public:
  //CONSTRUCTORS
  TimecourseDataSet();
  TimecourseDataSet(string dataFile);
  TimecourseDataSet(const vector<vector<double> >& inputData);
  
  //OTHER METHODS
  virtual double SingleClusterLogEvidence(const vector<int>& itemIndex,
					double& lengthScale,
					double& noiseFreeScale,
					double& noiseSigma,
					double& mixtureComponent) = 0;
  void ReadInData(string dataFile);
  void ReadInNoise(string dataFile);
  void ReadInNoise(vector<double> noise);
  void ReadInTimePoints(vector<double> timePoints);
  double GetClusterNoise(int nodeID);
  double GetMergedClusterNoise(vector<int>);
  double GetClusterSE(vector<int> ItemIndex);
  double GetClusterSEMS(vector<int> ItemIndex);
  //double GetMLIINoise(vector<int> itemIndex);
  void SetNoiseMode(int mode);
  int GetNoiseMode();
  void SetRobustMode(int mode);
  int GetRobustMode();
  void SetReps(int num_reps);
  vector<double> GetDataForCluster(vector<int> itemIndex);
  void SetDataType(string type);

  //TAGS
  int nTimePoints;
  double globalNoiseSigma; // estimate

 protected:
  //METHODS
  virtual double ComputeLogDeterminant(double* choleskyMatrix, int nVariables);
  virtual double ComputeLogEvidence(BlockCovarianceMatrix blockMatrix, vector<double> data);
  double ComputeMaximisedLogEvidence(vector<double> yValues,
				   double *lengthScale,
				   double *NoiseFree,
				   double *Noise);
  virtual double ComputePartRobustLogEvidenceMissingSingleObservation(double detCovarFunctionk,
								    BlockCovarianceMatrix invblockMatrix,
								    vector<double> yValsMissingSingleObservation);
  void OptimiseHyperparameters(const vector<double>& yValues,
			       double& lengthScale,
			       double& noiseFreeScale,
			       double& noiseSigma);
  void OptimiseHyperparametersFixedNoise(vector<double> yValues,
					 double& lengthScale,
					 double& noiseFreeScale,
					 double& noiseSigma);
  void OptimiseHyperparametersEstimatedNoise(vector<double> yValues,
					     double& lengthScale,
					     double& noiseFreeScale,
					     double& noiseSigma,
					     double replicateNoise);
  BlockCovarianceMatrix
    SquareExponentialCovarianceFunctionMissOneTimePoint(double lengthScale,
							int blockSize,
							double noiseFreeScale,
							int timePoint);
  BlockCovarianceMatrix SquareExponential_lengthDerivative(double lengthScale,
							   int blockSize,
							   double noiseFreeScale);
  BlockCovarianceMatrix AddNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix, 
						     double noiseSigma);
  BlockCovarianceMatrix AddFixedNoiseToCovarianceFunction(BlockCovarianceMatrix blockMatrix,
							  double noise_std_error);
  double ComputeGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
		       const BlockCovarianceMatrix& covarianceDerivative,
		       const vector<double>& alpha);
  double ComputeNoiseGradient(const BlockCovarianceMatrix& inverseCovarianceFunction,
			    const vector<double>& alpha,
			    double noiseSigma);
  double GetAlpha(double beta, double sigmaSEM);
  double GammaDistribution(double input, double alpha, double beta);
  double ComputeNoiseGradientIncludingGammaHyperprior(BlockCovarianceMatrix inverseCovarianceFunction,
						    vector<double> alpha,
						    double noiseSigma,
						    vector<double> gammaParams);
  vector<double> OptimiseGammaParams(double sigmaSEM);

  BlockCovarianceMatrix CovarianceFunctionMissingSingleObservation(vector<vector<double> > KnFC,
								   vector <double> KnC,
								   int KblockSize,
								   int KnRank,
								   int timePoint);

  virtual double ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
						     const int blockSize,
						     const vector<double>& params) = 0;
  virtual void ComputeGradientsFromHyperparameters(const vector<double>& yValues,
						   const int blockSize,
						   const vector<double>& params,
						   vector<double>& grad) = 0;
  virtual void
    ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						      const int blockSize,
						      const vector<double>& params,
						      double& logEv,
						      vector<double>& grad) = 0;

  virtual void ImposeConstraintsOnHyperparameters(vector<double>& params) = 0;
  
  void LineSearch(vector<double>& xold, const double fold, vector<double>& g,
		  vector<double>& p, vector<double>& x, double& f, const double stpmax,
		  bool& check, const int blockSize, const vector<double>& yValues);

  void DFPMaximise(vector<double>& p, const vector<int>& fix, const double gtol, double& fret,
		   const int blockSize, const vector<double>& yValues);
  
  //TAGS
 protected:
  //array to store the time-course data; this will be nDataItems * nTimePoints
  vector<vector<double> > data;
  //vector to store the time for each point in the time-courses
  vector<double> timePoints;
  //the range of the y data ymax-ymin
  double dataRange;
};

#endif // TIMECOURSEDATASET_H
