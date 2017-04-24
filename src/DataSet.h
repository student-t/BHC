/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef DATASET_H
#define DATASET_H

#include "header.h"
#include "BlockCovarianceMatrix.h"

/* ----------------------------------------------------------------------
   This is an abstract base class which allows BHC to handle different
   types of data. Note that since the virtual functions have not been
   implemented here, then DataSet cannot be instantiated directly.
---------------------------------------------------------------------- */

class DataSet
{
 public:
  virtual ~DataSet() {}
  virtual double SingleClusterLogEvidence(const vector<int>& itemIndex,
					  double& lengthScale,
					  double& noiseFreeScale,
					  double& noiseSigma,
					  double& mixtureComponent) = 0;
  virtual void ReadInData(string dataFile) = 0;
  virtual void ReadInNoise(string dataFile) = 0;
  virtual void ReadInNoise(vector<double> noise) = 0;
  virtual void ReadInTimePoints(vector<double> timePoints) = 0;
  virtual double GetClusterNoise(int nodeID) = 0;
  virtual double GetMergedClusterNoise(vector<int>) = 0;
  virtual double GetClusterSE(vector<int>) = 0;
  virtual double GetClusterSEMS(vector<int>) = 0;
  virtual void SetNoiseMode(int mode) = 0;
  virtual void SetReps(int num_reps) = 0;
  virtual int GetNoiseMode() = 0;
  virtual void SetRobustMode(int mode) = 0;
  virtual int GetRobustMode() = 0;
  virtual void SetDataType(string type) = 0;
  virtual double GetMLIINoise(vector<int> itemIndex) = 0;
  virtual BlockCovarianceMatrix SquareExponentialCovarianceFunctionMissingSingleObservation(vector<vector<double> > KnFC, vector<double> KnC, int KblockSize, int KnRank, int timePoint) = 0;
  virtual vector<double> GetDataForCluster(vector<int> itemIndex) = 0;

  void FindDataSize(string dataFile);
  void ReadInDataVector(vector<int> inputData, int nDataItems, int nFeatures);
  int Get_nDataItems();
  int Get_nFeatures();
  
  
  int nDataItems; // the number of data items stored by DataSet
  int nFeatures; // the number of features possessed by each data item
  int noise_mode; // indicates whether there is precalculated fixed noise
  int robust_mode;
  int reps;//the number of replicates per observation
  string dataType;
  vector<double>  noiseData; // noise values (computed from estimators)
};
#endif

