/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef CUBICSPLINETIMECOURSEDATASET_H
#define CUBICSPLINETIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>

#include "TimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

class CubicSplineTimecourseDataSet : public TimecourseDataSet
{
 public:
  CubicSplineTimecourseDataSet();
  CubicSplineTimecourseDataSet(string dataFile);
  CubicSplineTimecourseDataSet(const vector<vector<double> >& inputData);
  double GetMLIINoise(const vector<int>& itemIndex);
  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);
  
 protected:
  double ComputeMaximisedLogEvidence(const vector<double>& yValues,
				   double& noiseFreeScale,
				   double& noiseSigma);
  double CalculateFittedNoiseHyperparameter(const vector<double>& yValues);
  void OptimiseHyperparameters(const vector<double>& yValues,
			       double& noiseFreeScale,
			       double& noiseSigma);
  void OptimiseHyperparametersFixedNoise(const vector<double>& yValues,
					 double& noiseFreeScale,
					 double& noiseSigma);
  void OptimiseHyperparametersEstimatedNoise(const vector<double>& yValues,
					     double& noiseFreeScale,
					     double& noiseSigma,
					     double replicateNoise );
  BlockCovarianceMatrix CubicSplineCovarianceFunction(int blockSize,
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
#endif
