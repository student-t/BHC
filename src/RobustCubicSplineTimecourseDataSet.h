/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef ROBUSTCUBICSPLINETIMECOURSEDATASET_H
#define ROBUSTCUBICSPLINETIMECOURSEDATASET_H

#include <iostream>
#include <string>
#include <vector>

#include "CubicSplineTimecourseDataSet.h"
#include "BlockCovarianceMatrix.h"

class RobustCubicSplineTimecourseDataSet : public CubicSplineTimecourseDataSet
{
 public:
  RobustCubicSplineTimecourseDataSet();
  RobustCubicSplineTimecourseDataSet(string dataFile);
  RobustCubicSplineTimecourseDataSet(const vector<vector<double> >& inputData);

  double SingleClusterLogEvidence(const vector<int>& itemIndex,
				double& lengthScale,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);
 protected:
  double ComputeRobustLogEvidence(const vector<double>& yValues,
				int nCurrentItems,
				double& noiseFreeScale,
				double& noiseSigma,
				double& mixtureComponent);
};
#endif
