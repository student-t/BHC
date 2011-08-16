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
  vector<Node> GreedyClustering(DataSet* dataSet, bool verbose=false);
  
 private:
  vector<int> FindMostLikelyMerger(vector<Node>& treeNode);
};
#endif // DIRICHLETPROCESSMIXTURE_H
