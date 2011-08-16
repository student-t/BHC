//Need this code to make sure we only use this once during compilation
#ifndef ROBUSTCUBICSPLINETIMECOURSEDATASET_H
#define ROBUSTCUBICSPLINETIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>
#include "CubicSplineTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

// This class inherits the superclass DataSet.
// NOTE: in SplineTimecourseDataSet, nFeatures is deprecated (in favour of nTimePoints)
class RobustCubicSplineTimecourseDataSet : public CubicSplineTimecourseDataSet
{
 public:
  //CONSTRUCTORS
  RobustCubicSplineTimecourseDataSet();
  RobustCubicSplineTimecourseDataSet(string dataFile);
  RobustCubicSplineTimecourseDataSet(const vector<vector<double> >& inputData);


  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);

 protected:
  //METHODS
  double ComputeRobustLogEvidence(const vector<double>& yValues,
				int nCurrentItems,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);
};
#endif
