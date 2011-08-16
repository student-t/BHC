//Need this code to make sure we only use this once during compilation
#ifndef NODE_H
#define NODE_H
//some standard libraries we'd like
#include "header.h"
//include any object classes we need to access
#include "DataSet.h"
//define the class itself
class Node
{
public:
  //CONSTRUCTORS
  Node();                                        //default
  Node(DataSet* dataSet, vector<Node> treeNode); //initialize the node with one item
  Node(DataSet *dataSet, vector<Node> treeNode, vector<int> index);//merge two nodes into a new node
  //GETS AND SETS
  double      GetMergeProbability(int index);
  double      GetGlobalLogEvidence();
  double      GetClusterLogEvidence();
  double	  GetClusterLengthScale(int index);
  double	  GetClusterNoiseFree(int index);
  double	  GetClusterNoise(int index);
  double	  GetClusterMixtureComponent(int index);
  double      GetLogEvidenceForNewCluster(int index);
  vector<int> GetItemIndex();
  int         GetLeftChildIndex();
  int         GetRightChildIndex();
  void        SetNodeID(int inputID);
  void        SetMergeFlag(int inputFlag);
  //OTHER METHODS
  void        ComputeMergeLogProbability(DataSet *dataSet, Node newNode);
  void        OutputResultsToFile(vector<Node> treeNode, string outputFileName);
  int         AllowedToMerge();
  void        FlagAsMerged();

private:
  int             nodeID;              //a unique numerical identifier for the node
  vector<int>     itemIndex;           //indices for the data items in this cluster
  double          d_k;                 //used in the construction of the merging prior (d_k in the paper) (formerly ck)
  double          mergePrior;          //pi_k: the prior used when deciding mergers (see paper) (formerly nk)
  vector<double>  num1;                //logEV = num1 - num2
  vector<double>  num2;                //logEV = num1 - num2
  vector<double>  mergeLogProbability; //array.  posterior probability of merger with each other node. (formerly wt)
  double          clusterLogEvidence;  //scalar. the logEvidence for this single mixture component (i.e. favouring the merger over partitioned data)
  double          logEvidence;         //scalar.  Lower bound of the logEvidence of the DP mixture model. (formerly den)
  int             leftChildIndex;      //index of the left child node
  int             rightChildIndex;     //index of the right child node
  int             mergeFlag;           //has this node been merged yet?
  int             nItems;              //this counts the number of items in the cluster (formerly vec_no)

  double          dirichletProcessParameter; //this was called 'alp' in the old codes
  double	  Noise; //noise for this node
  double	  lengthScale;
  double	  NoiseFree;
  double	  mixtureComponent;
  vector<double>  LogEvidenceForNewCluster; //the SingleClusterLogEvidence of merger with each other node
  vector<double>  clusterlengthScale; //lengthScale for this node merged with indexed nodes
  vector<double>  clusterNoiseFree; //NoiseFree for this node merged with indexed nodes
  vector<double>  clusterMixtureComponent; //mixture component for the Robust case, default is 1 if robust_mode = 0

protected:
  vector<double>  clusterNoise; //Noise for this node merged with indexed nodes
};
#endif

