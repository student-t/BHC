/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef ROBUSTSQUAREDEXPONENTIALTIMECOURSEDATASET_H
#define ROBUSTSQUAREDEXPONENTIALTIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>

#include "SquaredExponentialTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

class RobustSquaredExponentialTimecourseDataSet : public SquaredExponentialTimecourseDataSet
{
 public:
  RobustSquaredExponentialTimecourseDataSet();
  RobustSquaredExponentialTimecourseDataSet(string dataFile);
  RobustSquaredExponentialTimecourseDataSet(const vector<vector<double> >& inputData);
  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);

 protected:
  double ComputeRobustLogEvidence(vector<double> yValues,
				int nCurrentItems,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);
};
#endif
