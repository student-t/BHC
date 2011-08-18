#include "DirichletProcessMixture.h"
#include <limits>



// DEFAULT CONSTRUCTOR
DirichletProcessMixture::DirichletProcessMixture() {}



// RUN THE GREEDY CLUSTERING ALGORITHM
// verbose has a default value of 'false', set in the header file
vector<Node> DirichletProcessMixture::GreedyClustering(DataSet* dataSet, bool verbose)
{
  //DECLARATIONS
  vector<Node> treeNode; //vector of Node objects that we use for the BHC greedy algorithm
  vector<int> mergeIndex; //holds the indices of Nodes to be merged
  Node mergedNode;
  int i;

  // We want 4 decimal places of output for the log ev's
  cout.precision(4);
  
  //INSTANTIATE EACH STARTING (SINGLE ITEM) CLUSTER
  for (i=0; i<dataSet->nDataItems; i++)
  {
    treeNode.push_back(Node(dataSet, treeNode));
  }

  if(verbose)
    {
      cout << "BHC START ====================================" << endl;
      cout << "Number of data items:\t" << dataSet->nDataItems << endl;
      cout << "Number of features:\t" << dataSet->nFeatures << endl;
      cout << "Algorithm:\t\tGreedy clustering" << endl;
      cout << "\nNode\tMerges\t\tLog Ev" << endl;
      cout << "----\t------\t\t------" << endl;
    }

  // w is the max width (number of digits) for the node indices,
  // this is used to keep the output below tidy and aligned
  int w=1;
  int mult=10;
  while(2*(dataSet->nDataItems)-1 >= mult)
    {
      ++w;
      mult*=10;
    }
  
  //PERFORM MERGERS TO FORM HIERARCHICAL CLUSTERING STRUCTURE
  for (i=dataSet->nDataItems; i<2*(dataSet->nDataItems)-1; i++)
  {
    mergeIndex = FindMostLikelyMerger(treeNode);
    
    // Construct the merged node
    mergedNode = Node(dataSet, treeNode, mergeIndex);
    treeNode.push_back(mergedNode);

    // Flag the other nodes as merged
    treeNode[mergeIndex[0]].FlagAsMerged();
    treeNode[mergeIndex[1]].FlagAsMerged();

    //PRINT THE MERGER RESULTS TO SCREEN
    if(verbose)
      {
	cout << std::setw(w) << (i+1)
	     << "\t(" << std::setw(w) << mergeIndex[0]+1 << ","
	     << std::setw(w) << mergeIndex[1]+1 << ")\t\t"
	     << fixed << treeNode[mergeIndex[0]].GetMergeProbability(mergeIndex[1])
	     << endl;
      }
    
    //THE OLD OUTPUT:
    /*
    cout << "Node number:    " << i+1
	 << "\t MergeIndex1: " << mergeIndex[0]+1
	 << "\t MergeIndex2: " << mergeIndex[1]+1
	 << "\t LogEv: " << treeNode[mergeIndex[0]].GetMergeProbability(mergeIndex[1]);
    if (dataSet->dataType =="time-course")
    {
      cout << "\t LS: "	<< treeNode[mergeIndex[0]].GetClusterLengthScale(mergeIndex[1]);
    }
    cout << "\t NF: " << treeNode[mergeIndex[0]].GetClusterNoiseFree(mergeIndex[1])
	 << "\t N: "  << treeNode[mergeIndex[0]].GetClusterNoise(mergeIndex[1]);
    if (dataSet->GetNoiseMode() == 2)
    {
      cout << "\t Fitted N:" << dataSet->GetMLIINoise(mergedNode.GetItemIndex());
    }
    if (dataSet->GetRobustMode() ==1)
    {
      cout << "\t Mixture Component:"
	   << treeNode[mergeIndex[0]].GetClusterMixtureComponent(mergeIndex[1]);
    }
    if (dataSet->noise_mode == 2)
    {
      cout << "\t SEMS: " << dataSet->GetClusterSEMS(mergedNode.GetItemIndex());
    }
    cout << "\t SES: " << dataSet->GetClusterSE(mergedNode.GetItemIndex());
    cout << "\t SingleClusLogEv: "
	 << treeNode[mergeIndex[0]].GetLogEvidenceForNewCluster(mergeIndex[1])
	 << "\t ClusSize: " << mergedNode.GetItemIndex().size() << "\t" ;
    //FIND THE GENES IN THE MERGED CLUSTER
    for (k=0; k<mergedNode.GetItemIndex().size(); k++)
    {
      cout << mergedNode.GetItemIndex()[k]+1 << "," ;
    }
    cout << endl;
    */
  }

  if(verbose)
    {
      //PRINT OUT THE LOWER BOUND ON THE OVERALL MODEL LOG EVIDENCE
      cout << "\n(Lower bound) model Log Ev:\t"
	   << fixed << mergedNode.GetGlobalLogEvidence() << endl;
      cout.precision(2);
      cout << "BHC END ======================================" << endl;
    }
  
  return treeNode;
}



//FIND THE MOST LIKELY MERGER, GIVEN THE CURRENT SET OF NODES
vector<int> DirichletProcessMixture::FindMostLikelyMerger(vector<Node>& treeNode)
{
  //DECLARATIONS
  int i, j;
  double currentLogEvidence, maxLogEvidence= -numeric_limits<double>::infinity();
  vector<int> mergeIndex(2,-1);//initialise to disallowed values

  //SCAN EACH NODE IN TURN, LOOKING FOR THE MOST LIKELY MERGER
  //REMEMBERING TO NOT INCLUDE ALREADY-MERGED NODES
  for (i=1; i<int(treeNode.size()); i++)
    for (j=0; j<i; j++)
    {
      if (treeNode[i].AllowedToMerge() && treeNode[j].AllowedToMerge())
      {
        currentLogEvidence = treeNode[i].GetMergeProbability(j);
        if (currentLogEvidence >= maxLogEvidence)
        {
          mergeIndex[0] = i;
          mergeIndex[1] = j;
          maxLogEvidence = currentLogEvidence;
        }
      }
    }
  
  return mergeIndex;
}
