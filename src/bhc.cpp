#include "header.h"
/*
(implemented March '08, Rich Savage)
(C++ code originally written by Xang Yu; modified by Rich Savage)

Bayesian Hierarchical Clustering (BHC) 
--------------------------------------
This is an R implementation of the BHC algorithm developed by Heller & Ghahramani (2005), set up to handle multinomial data.
This R function links to C++ code originally authored by Yang Xu and modified to the multinomial case by Rich Savage.
That code was based on original (Matlab) code written by Katherine Heller.

This is a wrapper function for the BHC C++ code.  It allows us to link properly to R and handles things like the passing of 
variables between C++ and R.
 */
extern "C" { //required for linking properly to R
  void bhcWrapper(int* inputData, int* nObservations, int* nFeatures, double* ghpInput, int* nFeatureValues_input, 
		  double* logEvidence, int* node1, int* node2, int* mergeOrder, double* mergeWeight) {
    //----------------------------------------------------------------------
    // DECLARATIONS --------------------------------------------------------
    //----------------------------------------------------------------------
    NODE*  tr_node;
    int    i=0, j=0;
    double alp=0.001, minWeight= -100000.;
    int    obs=*nObservations, dim=*nFeatures, nFeatureValues=*nFeatureValues_input;
    double globalHyperParameter= *ghpInput; //do this as our inputs/outputs are pointers (interface)
    //----------------------------------------------------------------------
    // READ IN DATA; RUN CLUSTERING ANALYSIS -------------------------------
    //----------------------------------------------------------------------
    tr_node      = ReadInData(dim, obs, minWeight, nFeatureValues, inputData);
    *logEvidence = bayeslink_binf(tr_node,dim,obs,globalHyperParameter, alp, minWeight, nFeatureValues);
    //----------------------------------------------------------------------
    // PUT DENDROGRAM DATA IN A FORM TO PASS BACK TO R ---------------------
    //----------------------------------------------------------------------
    for (i=obs;i<obs*2-1;i++){
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
//*****************************************************************************
//*****************************************************************************
//----------------------------------------------------------------------
// ----------------------------------------
//----------------------------------------------------------------------
