/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef SQUAREDEXPONENTIALTIMECOURSEDATASET_H
#define SQUAREDEXPONENTIALTIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>

#include "header.h"
#include "TimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

class SquaredExponentialTimecourseDataSet : public TimecourseDataSet
{
 public:
  SquaredExponentialTimecourseDataSet();
  SquaredExponentialTimecourseDataSet(string dataFile);
  SquaredExponentialTimecourseDataSet(const vector<vector<double> >& inputData);
  
  double GetMLIINoise(vector<int> itemIndex);
  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				  double& lengthScale,
				  double& NoiseFree,
				  double& Noise,
				  double& mixtureComponent);
 protected:
  double ComputeMaximisedLogEvidence(vector<double> yValues,
				     double& lengthScale,
				     double& noiseFreeScale,
				     double& noiseSigma);
  double CalculateFittedNoiseHyperparameter(const vector<double>& yValues);
  void OptimiseHyperparameters(const vector<double>& yValues,
			       double& lengthScale,
			       double& NoiseFree,
			       double& Noise);
  void OptimiseHyperparametersFixedNoise(vector<double> yValues,
					 double& lengthScale,
					 double& NoiseFree,
					 double& Noise);
  void OptimiseHyperparametersEstimatedNoise(vector<double> yValues,
					     double& p_lengthScale,
					     double& p_NoiseFree,
					     double& p_Noise,
					     double replicateNoise);
  BlockCovarianceMatrix SquareExponentialCovarianceFunction(double lengthScale,
							    int blockSize,
							    double noiseFreeScale);
  BlockCovarianceMatrix SquareExponential_lengthDerivative(double lengthScale,
							   int blockSize,
							   double noiseFreeScale);
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
