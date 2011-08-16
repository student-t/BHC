#ifndef SQUAREDEXPONENTIALTIMECOURSEDATASET_H
#define SQUAREDEXPONENTIALTIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>
#include "header.h"
#include "TimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

// This class inherits the superclass DataSet.
// NOTE: in TimecourseDataSet, nFeatures is deprecated (in favour of nTimePoints)
class SquaredExponentialTimecourseDataSet : public TimecourseDataSet
{

 public:
  //CONSTRUCTORS
  SquaredExponentialTimecourseDataSet();
  SquaredExponentialTimecourseDataSet(string dataFile);
  SquaredExponentialTimecourseDataSet(const vector<vector<double> >& inputData);
  double GetMLIINoise(vector<int> itemIndex);

  //OTHER METHODS
  double SingleClusterLogEvidence(const vector<int>& itemIndex, double& lengthScale, double& NoiseFree, double& Noise, double& mixtureComponent);  // returns the log-evidence for a cluster that contains data items defined by itemIndex
 protected:
  //METHODS
  double ComputeMaximisedLogEvidence(vector<double> yValues, double& lengthScale, double& noiseFreeScale, double& noiseSigma);
  double CalculateFittedNoiseHyperparameter(const vector<double>& yValues);
  void OptimiseHyperparameters(const vector<double>& yValues, double& lengthScale, double& NoiseFree, double& Noise );
  void OptimiseHyperparametersFixedNoise(vector<double> yValues, double& lengthScale, double& NoiseFree, double& Noise);
  void OptimiseHyperparametersEstimatedNoise(vector<double> yValues, double& p_lengthScale, double& p_NoiseFree, double& p_Noise, double replicateNoise);
  BlockCovarianceMatrix SquareExponentialCovarianceFunction(double lengthScale, int blockSize, double noiseFreeScale);

  BlockCovarianceMatrix SquareExponential_lengthDerivative(double lengthScale, int blockSize, double noiseFreeScale);

  double ComputeLogEvidenceFromHyperparameters(const vector<double>& yValues,
					     const int blockSize,
					     const vector<double>& params);

  void ComputeGradientsFromHyperparameters(const vector<double>& yValues,
					   const int blockSize,
					   const vector<double>& params,
					   vector<double>& grad);
  void
    ComputeLogEvidenceAndGradientsFromHyperparameters(const vector<double>& yValues,
						      const int blockSize,
						      const vector<double>& params,
						      double& logEv,
						      vector<double>& grad);

  void ImposeConstraintsOnHyperparameters(vector<double>& params);
};

#endif // SQUAREDEXPONENTIALTIMECOURSEDATASET_H
