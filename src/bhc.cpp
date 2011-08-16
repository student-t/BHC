#include "header.h"
#include <limits>
#include "MultinomialDataSet.h"
#include "SquaredExponentialTimecourseDataSet.h"
#include "RobustSquaredExponentialTimecourseDataSet.h"
#include "CubicSplineTimecourseDataSet.h"
#include "RobustCubicSplineTimecourseDataSet.h"
#include "Node.h"
#include "DirichletProcessMixture.h"
#include "BlockCovarianceMatrix.h"
// 1.1.0 includes
#include "multinomial_header.h"


//
// --------------------------------------
// Bayesian Hierarchical Clustering (BHC)
// --------------------------------------
//
// This is a C++ implementation of the BHC algorithm developed by Heller & Ghahramani
// (2005), and set up to handle multinomial data. The bhcWrapper() function below
// acts as the interface between the C++ code and R. Input and output data are
// exchanged via the function arguments.
//
// The C++ code was originally authored by Yang Xu and modified to the multinomial
// case by Rich Savage (Oct 2008) -- that code being based on original (Matlab) code
// written by Katherine Heller. It has been further developed by Emma Cooke, and
// optimised by Robert Darkins (2011).
//
// This version includes a version of the BHC algorithm that runs on time-series data.
//
// DataTypeID should be interpreted as follows:
//    0 - multinomial
//    1 - time-series (using squared-exponential covariance)
//    2 - cubicspline (using cubic spline covariance)
//
// The inputData is ordered as follows:
//    item1_feature1  item1_feature2 ... item2_feature1 item2_feature2 ... etc
//
// The DataSet classes expect a 2D array of data (nItems * nFeatures).
//

extern "C" { //required for linking properly to R
  void bhcWrapper(int* dataTypeID_input, double* inputData, double* timePoints,
		  int* nDataItems_input, int* nFeatures_input, double* ghpInput,
		  double* noise_input, int* num_reps, int* set_noise_input,
		  int* robust_input, int* nFeatureValues_input, double* logEvidence,
		  int* node1, int* node2, int* mergeOrder, double* mergeWeight,
		  int* numThreads)
  {
    //----------------------------------------------------------------------
    // DECLARATIONS --------------------------------------------------------
    //----------------------------------------------------------------------
    int i, j, counter1=0;
    string dataFile, dataType, outputFile;
    vector<Node> treeNode; //array of Node objects that we use for the BHC greedy algorithm
    DataSet* dataSet=NULL;  //this is a super-class, so we can use different types of data here
    DirichletProcessMixture bhc;
    //do the following assignment since our inputs/outputs are pointers
    int nDataItems=*nDataItems_input, 
      nFeatures=*nFeatures_input,
      dataTypeID=*dataTypeID_input,
      num_reps_copy = *num_reps,
      set_noise_copy = *set_noise_input,
      robust_copy = *robust_input;
    double globalHyperParameter= *ghpInput;
    vector<double> noise, timePoints_copy;
    vector<vector<int> > data_int;//for ease, we define both these 2D arrays; we only need one
    vector<vector<double> > data_double;//

    //----------------------------------------------------------------------
    // PARALLELISE IF REQUESTED --------------------------------------------
    //----------------------------------------------------------------------
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(MAX(*numThreads,1));
#endif
    
    //----------------------------------------------------------------------
    // DECIDE ON THE TYPE OF DATE WE'RE ANALYSING --------------------------
    //----------------------------------------------------------------------
    //this is doublely just so we have a more informative variable than just an integer ID
    if      (dataTypeID==0) dataType = "multinomial";
    else if (dataTypeID==1) dataType = "time-course";
    else if (dataTypeID==2) dataType = "cubicspline";
    else
    {
      cout<<"Oops! invalid dataType! tag one"<<endl;
      return;
    }

    //----------------------------------------------------------------------
    // CONVERT THE DATA BACK TO A 2D ARRAY OF THE RELEVANT TYPE ------------
    //----------------------------------------------------------------------
    //we choose to do this here so that the BHC code doesn't have to know about the
    //R interface.  We could also do this via class constructors.
    if (dataType=="multinomial")
    {
      for (i=0; i<nDataItems; i++)
      {
        data_int.push_back(vector<int>(nFeatures, 0));
        for (j=0; j<nFeatures; j++)
        {
          data_int[i][j] = inputData[counter1++];
        }
      }

    }

    //the only difference between here and the above code is we're copying data
    //into an array of doubles, not ints
    else if (dataType=="time-course" || dataType=="cubicspline")
    {
      for (i=0; i<nDataItems; i++)
      {
        data_double.push_back(vector<double>(nFeatures, 0));
        for (j=0; j<nFeatures; j++)
        {
          data_double[i][j] = inputData[counter1++];
        }
      }
    }

    //Read in the timePoints
    if (dataType=="time-course" || dataType=="cubicspline")
    {
      for (i=0; i<nFeatures; i++)
      {
        timePoints_copy.push_back(timePoints[i]);
      }
    }
    
    //----------------------------------------------------------------------
    // COPYING OVER THE NOISE VECTOR ---------------------------------------
    //----------------------------------------------------------------------
    if (set_noise_copy ==2)
    {
      for (i=0; i<nDataItems; i++)
      {
        noise.push_back(noise_input[i]);
      }
    }
    else if (set_noise_copy ==1)
    {
      noise.push_back(noise_input[0]);
    }
    else if (set_noise_copy ==0)
    {
      noise.push_back(0.0);
    }

    //----------------------------------------------------------------------
    // INSTANTIATE THE REQUIRED TYPE OF DATASET OBJECT ---------------------
    //----------------------------------------------------------------------
    if (dataType=="multinomial")
      dataSet = new MultinomialDataSet(data_int, globalHyperParameter);
    else if (dataType=="time-course")
    {
      if (robust_copy == 0)
      {
        //cout << "time-course" << endl;
        dataSet = new SquaredExponentialTimecourseDataSet(data_double);
        dataSet->ReadInNoise(noise);
      }
      else if (robust_copy == 1)
      {
        //cout << "Robust time-course" << endl;
        dataSet = new RobustSquaredExponentialTimecourseDataSet(data_double);
        dataSet->ReadInNoise(noise);
      }
    }
    else if (dataType=="cubicspline")
    {
      if (robust_copy == 0)
      {
        dataSet = new CubicSplineTimecourseDataSet(data_double);
        dataSet->ReadInNoise(noise);
      }
      if (robust_copy == 1)
      {
        dataSet = new RobustCubicSplineTimecourseDataSet(data_double);
        dataSet->ReadInNoise(noise);
      }
    }
    else
    {
      cout<<"Oops! invalid dataType! tag two"<<endl;
      return;
    }

    dataSet->SetNoiseMode(set_noise_copy);
    if (set_noise_copy ==2)
    {
      dataSet->SetReps(num_reps_copy);
    }
    dataSet->SetRobustMode(robust_copy);
    dataSet->SetDataType(dataType);
    if(dataTypeID!=0)
      dataSet->ReadInTimePoints(timePoints_copy); // no time points for multinomial case

    //----------------------------------------------------------------------
    // RUN CLUSTERING ANALYSIS ---------------------------------------------
    //----------------------------------------------------------------------
    treeNode = bhc.GreedyClustering(dataSet,true);

    //----------------------------------------------------------------------
    // PUT DENDROGRAM DATA IN A FORM TO PASS BACK TO R ---------------------
    //----------------------------------------------------------------------
    for (i=0; i<nDataItems-1; i++)
    {
      node1[i] = treeNode[i+nDataItems].GetLeftChildIndex()  + 1;
      node2[i] = treeNode[i+nDataItems].GetRightChildIndex() + 1;
      mergeOrder[i] = i + 1;
      mergeWeight[i] = treeNode[i+nDataItems].GetClusterLogEvidence();
    }
    *logEvidence = treeNode[treeNode.size()-1].GetGlobalLogEvidence();//ugly way to access the last element!
  }


  void bhcWrapper_multinomial(int* inputData, int* nObservations, int* nFeatures,
			      double* ghpInput, int* nFeatureValues_input,
			      double* logEvidence, int* node1, int* node2,
			      int* mergeOrder, double* mergeWeight)
  {
    //----------------------------------------------------------------------
    // DECLARATIONS --------------------------------------------------------
    //----------------------------------------------------------------------
    NODE*  tr_node;
    int    i=0, j=0;
    double alp=0.001, minWeight= -numeric_limits<double>::infinity();
    int    obs=*nObservations, dim=*nFeatures, nFeatureValues=*nFeatureValues_input;
    double globalHyperParameter= *ghpInput; //do this as our inputs/outputs are pointers
    
    //----------------------------------------------------------------------
    // READ IN DATA; RUN CLUSTERING ANALYSIS -------------------------------
    //----------------------------------------------------------------------
    tr_node      = ReadInData(dim, obs, minWeight, nFeatureValues, inputData);
    *logEvidence = bayeslink_binf(tr_node,dim,obs,globalHyperParameter, alp,
				  minWeight, nFeatureValues);
    
    //----------------------------------------------------------------------
    // PUT DENDROGRAM DATA IN A FORM TO PASS BACK TO R ---------------------
    //----------------------------------------------------------------------
    for (i=obs;i<obs*2-1;i++)
      {
	node1[i-obs]       = tr_node[i].pleft+1;
	node2[i-obs]       = tr_node[i].pright+1;
	mergeOrder[i-obs]  = i-obs+1;
	mergeWeight[i-obs] = tr_node[i].wt[i];
      }
    
    //----------------------------------------------------------------------
    // FREE ALLOCATED MEMORY -----------------------------------------------
    //----------------------------------------------------------------------
    for (i=0;i<2*obs;i++) {
      for (j=0;j<nFeatureValues;j++)
	delete[] tr_node[i].dat[j];
      delete[] tr_node[i].num1;
      delete[] tr_node[i].num2;
      delete[] tr_node[i].wt;
    }
    delete[] tr_node;
  }
}
