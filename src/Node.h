/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef NODE_H
#define NODE_H

#include "header.h"
#include "DataSet.h"

/* ----------------------------------------------------------------------
   This class describes a node within our binary tree and represents
   a gene or a cluster thereof.

   It was rewritten by R. Darkins (2011). A key difference between this
   and the old implementation is that this version does not store information
   regarding its relationship with other nodes in the tree (for instance,
   it does not store hyperparameters). This was because storing such info
   made the randomised algorithm almost impossible to implement with
   any efficiency.
------------------------------------------------------------------------- */

class Node
{
 public:
  Node(); 
  static Node CreateDataNode(DataSet& dataSet, const int arg_nodeID);
  static Node CreateMergerNode(DataSet& dataSet,
			       const Node& node1,
			       const Node& node2,
			       const int arg_nodeID=0);
  
  int GetDataID() const;
  int GetNodeID() const;
  int GetLeftChildNodeID() const;
  int GetRightChildNodeID() const;
  double GetClusterLogEvidence() const;
  double GetLowerBoundLogEvidence() const;
  int GetNumItems() const;
  bool AllowedToMerge() const;
  vector<int> GetChildDataIDs() const;
  vector<int> GetChildNodeIDs() const;
  void SetNodeID(int newID);
  void SetLeftChildNodeID(int newID);
  void SetRightChildNodeID(int newID);
  void FlagAsMerged();
  void UnflagAsMerged();
  void OutputResultsToFile(const vector<Node>& treeNode, string outputFileName);

 private:
  static const double dirichletProcessParameter;
  int nodeID;
  int dataID;
  int leftChildNodeID;
  int rightChildNodeID;
  vector<int> childNodeIDs;
  vector<int> childDataIDs;
  double log_d_k;
  bool mergeFlag;
  double clusterLogEvidence;
  double lowerBoundLogEvidence;
  double mergePrior;
};
#endif // NODE_H
