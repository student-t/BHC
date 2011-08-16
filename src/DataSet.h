//Need this code to make sure we only use this once during compilation
#ifndef DATASET_H
#define DATASET_H
//some standard libraries we'd like
#include "header.h"
#include "BlockCovarianceMatrix.h"
//define the class itself

// This is an abstract base-class that allows us to define different specific DataSet
// classes, so that BHC can handle different types of data.
// It is an example of polymorphism. Since the virtual functions have not been
// implemented here, then DataSet cannot be instantiated directly.
class DataSet
{
public:
  //VIRTUAL DESTRUCTOR (need this because we're using virtual functions)
  virtual ~DataSet() {}
  //PLACE-HOLDERS!  Any subclass must provide these method (correct for the type of data)
  virtual double SingleClusterLogEvidence(const vector<int>& itemIndex,
					double& lengthScale,
					double& noiseFreeScale,
					double& noiseSigma,
					double& mixtureComponent);
  virtual void   ReadInData(string dataFile);
  virtual void   ReadInNoise(string dataFile);
  virtual void   ReadInNoise(vector<double> noise);
  virtual void   ReadInTimePoints(vector<double> timePoints);
  virtual double   GetClusterNoise(int nodeID);
  virtual double   GetMergedClusterNoise(vector<int>);
  virtual double   GetClusterSE(vector<int>);
  virtual double   GetClusterSEMS(vector<int>);
  virtual void  SetNoiseMode(int mode);
  virtual void SetReps(int num_reps);
  virtual int GetNoiseMode();
  virtual void  SetRobustMode(int mode);
  virtual int GetRobustMode();
  virtual void  SetDataType(string type);
  virtual double GetMLIINoise(vector<int> itemIndex);
  virtual BlockCovarianceMatrix SquareExponentialCovarianceFunctionMissingSingleObservation(vector< vector<double> > KnFC, vector <double> KnC, int KblockSize, int KnRank, int timePoint);
  virtual vector<double> GetDataForCluster(vector<int> itemIndex);


  //OTHER METHODS
  void   FindDataSize(string dataFile);
  void   ReadInDataVector(vector<int> inputData, int nDataItems, int nFeatures);
  //GETS AND SETS
  int Get_nDataItems();
  int Get_nFeatures();
  //TAGS
  int nDataItems;                                         //the number of data items stored by DataSet
  int nFeatures;                                          //the number of features possessed by each data item
  int	 noise_mode;				   //indicates whether there is precalculated fixed noise
  int robust_mode;
  int reps;													//the number of replicates per observation
  string dataType;
  vector<double>		   noiseData;						//vector to store the noise values for the data (computed from estimators)

private:
};
#endif

