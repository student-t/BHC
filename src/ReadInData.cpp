#include "header.h"
//Routine to read in data for the BHC code
//For the R-interfaced version of this code, the data are passed to this function
//in a 1D array.  They're then copied into the NODE structures.
NODE* ReadInData(int dim, int obs, double min_wt, int nFeatureValues, int* inputData){
  //----------------------------------------------------------------------
  // DECLARATIONS --------------------------------------------------------
  //----------------------------------------------------------------------	
  NODE*   tr_node;
  int     i, j, k;
  int     inputValue, counter=0;
  //----------------------------------------------------------------------
  // ALLOCATE MEMORY, INITIALISE VALUES FOR THE DENDROGRAM NODES ---------
  //----------------------------------------------------------------------
  tr_node = new NODE[2*obs];
  for (i=0;i<2*obs;i++){
      tr_node[i].dat       = new double* [dim];
      for (j=0;j<dim;j++) 
	tr_node[i].dat[j]  = new double [nFeatureValues];
      tr_node[i].num1      = new double[2*obs];
      tr_node[i].num2      = new double[2*obs];
      tr_node[i].wt        = new double[2*obs];
      for (k=0;k<2*obs;k++) 
	(tr_node[i]).wt[k] = min_wt;
      tr_node[i].pleft     = -1;
      tr_node[i].pright    = -1;
      tr_node[i].flag      = 0;
  }
  //----------------------------------------------------------------------
  // INITIALISE THE DATA VALUES ------------------------------------------
  //----------------------------------------------------------------------
  for (i=0;i<obs;i++)
    for (j=0;j<dim;j++)
      for (k=0;k<nFeatureValues;k++)
	tr_node[i].dat[j][k] = 0;
  //----------------------------------------------------------------------
  // READ IN THE DATA VALUES FROM THE INPUT DATA ARRAY --------------------
  //----------------------------------------------------------------------
  for (i=0;i<obs;i++){
      for (j=0;j<dim;j++){
	inputValue = inputData[counter++];
	//check that this value is legal
	if (inputValue <0 | inputValue >= nFeatureValues)  
	  cout << "Error! Input value outside range " << i << " " << j << " " << inputValue << endl; 		     
	//if so, add one to the relevant counter
	//actually, we will only have one value here for each feature for each data item
	//therefore, we should just set the relevant value to unity
	tr_node[i].dat[j][int(inputValue)] = 1 ;
	}
      tr_node[i].vec_no = 1;
    }
  //----------------------------------------------------------------------
  // RETURN DATA NODES ---------------------------------------------------
  //----------------------------------------------------------------------
  return tr_node;
}
//*****************************************************************************
//*****************************************************************************
//----------------------------------------------------------------------
// ----------------------------------------
//----------------------------------------------------------------------
