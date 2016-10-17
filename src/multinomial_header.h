#include <iostream>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include "gammaln.h"

//Structure: NODE.
//Represents a node in the dendrogram.
struct NODE {
		double** dat;
		double   ck;
		double   nk;
		double*  num1;
		double*  num2;
		double*  wt;
		double   den;
		int      pleft;
		int      pright;
		int      flag;
		int      vec_no;
};

double bayeslink_binf(NODE* tr_node, int dim, int obs, 
		      double cc, double alp, double min_wt, int nFeatureValues);

double binevidence(NODE* tr_node, int dim, double** hyperParameters, int index1, int index2, int nFeatureValues);

NODE* ReadInData(int dim, int obs, double min_wt, int nFeatureValues, int* inputData);

double** CalculateHyperparameters(NODE* tr_node, int dim, int obs, int nFeatureValues, double cc);

void OutputDendrogramInformation(NODE* tr_node, int dim, int obs);
