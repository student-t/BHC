#ifndef DIRICHLETPROCESSMIXTURE_H
#define DIRICHLETPROCESSMIXTURE_H

#include "header.h"
#include "Node.h"

/*
  This class defines A Dirichlet Process (infinite mixture) model
*/
class DirichletProcessMixture
{
 public:
  DirichletProcessMixture();
  vector<Node> GreedyClustering(DataSet& dataSet,
				const bool verbose=false);
  vector<Node> RandomisedClustering(DataSet& dataSet,
				    const int m,
				    const int t,
				    const bool verbose);
  
 private:
  void RandomisedClustering(vector<Node>& D_sub,
			    vector<Node>& D,
			    DataSet& dataSet,
			    const int m,
			    const int t,
			    int& depth,
			    const bool verbose);
  void GreedyClustering(vector<Node>& treeNode,
			DataSet& dataSet,
			const bool verbose);
  vector<int> FindMostLikelyMerger(const vector<Node>& treeNode);
};

#endif // DIRICHLETPROCESSMIXTURE_H
