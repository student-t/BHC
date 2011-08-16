//Need this code to make sure we only use this once during compilation
#ifndef ROBUSTSQUAREDEXPONENTIALTIMECOURSEDATASET_H
#define ROBUSTSQUAREDEXPONENTIALTIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>
#include "SquaredExponentialTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"
/*
  This class inherits the superclass DataSet.
  We do this so that BHC can run with different types of data.
  (this is like using a Java interface)

  NOTE: in TimecourseDataSet, nFeatures is deprecated (in favour of nTimePoints)

*/
class RobustSquaredExponentialTimecourseDataSet : public SquaredExponentialTimecourseDataSet
{
 public:
  //CONSTRUCTORS
  RobustSquaredExponentialTimecourseDataSet();
  RobustSquaredExponentialTimecourseDataSet(string dataFile);
  RobustSquaredExponentialTimecourseDataSet(const vector<vector<double> >& inputData);
  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);

 protected:
  //METHODS

  double ComputeRobustLogEvidence(vector<double> yValues,
				int nCurrentItems,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);

};
#endif
