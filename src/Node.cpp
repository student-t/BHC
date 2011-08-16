#include "Node.h"
#include <limits>


//CONSTRUCTOR - NO DATA ITEMS
Node::Node() {}



//CONSTRUCTOR - ONE DATA ITEM (HENCE THE INPUT ISN'T A VECTOR)
//a node should only need to know about mergers with lower-ID other nodes
//Therefore, we pass this method a vector of the lower-ID nodes (if any)
Node::Node(DataSet *dataSet, vector<Node> treeNode)
{
  //DECLARATIONS
  int    i;
  vector<double> clusterData;
  //INITIALISE THE OBJECT VARIABLES
  nodeID                    = treeNode.size();//this makes it the n+1th node
  itemIndex.push_back(nodeID);//associated with the n+1th data item
  dirichletProcessParameter = 0.001;//for now, hardwire this to the value used in Heller, Xu's codes
  d_k                       = log(dirichletProcessParameter);//this is actually log(d_k) from the paper
  mergePrior                = 1;
  Noise			    = dataSet->GetClusterNoise(nodeID);
  logEvidence               = dataSet->SingleClusterLogEvidence(itemIndex, lengthScale, NoiseFree, Noise, mixtureComponent);//no point computing logevidence for a single gene
  leftChildIndex            = -1;
  rightChildIndex           = -1;
  mergeFlag                 = 0;
  nItems                    = 1;
  clusterlengthScale	= vector<double>(nodeID, 0);
  clusterNoiseFree = vector<double>(nodeID, 0);
  clusterNoise = vector<double>(nodeID, 0);
  clusterMixtureComponent = vector<double>(nodeID, 0);
  mergeFlag=0;
  //IF THIS ISN'T THE FIRST NODE, NEED SOME OTHER INITIALISATIONS
  if (nodeID>0)
  {
    //INITIALISE VECTORS TO HOLD MERGER LOG-PROBABILITIES PLUS USEFUL INTERMEDIATES
    num1                = vector<double>(nodeID, 0);
    num2                = vector<double>(nodeID, 0);
    mergeLogProbability = vector<double>(nodeID, -numeric_limits<double>::infinity());
    LogEvidenceForNewCluster = vector<double>(nodeID, 0);
    //CALCULATE THE MERGER PROBABILITIES
#pragma omp parallel for default(shared) schedule(dynamic,1)
    for (i=0; i<nodeID; i++)
      {
	ComputeMergeLogProbability(dataSet, treeNode[i]);
      }
  }
}



// Culprit
//CONSTRUCTOR - MERGE TWO EXISTING NODES TOGETHER
//when merging, we wish to construct a new node based on two child nodes
Node::Node(DataSet *dataSet, vector<Node> treeNode, vector<int> index)
{
  //DECLARATIONS
  int i;
  double tr1, tr2, a, b;//intermediates for computing weights
  Node node1, node2;
  //EXTRACT THE NODES THAT WE'RE GOING TO MERGE
  node1 = treeNode[index[0]];
  node2 = treeNode[index[1]];

  //INITIALISE/UPDATE SOME OBJECT TAGS
  dirichletProcessParameter = 0.001;//for now, hardwire this to the value used in Heller, Xu's codes
  mergeFlag                 = 0;
  nItems                    = node1.nItems + node2.nItems;
  mergePrior                = node1.mergePrior + node2.mergePrior;
  Noise			    = node1.clusterNoise[index[1]];  // Depends on node1.nodeID > node2.nodeID
  //STORE THE LOG-EVIDENCE FOR THIS MIXTURE COMPONENT
  clusterLogEvidence = node1.mergeLogProbability[node2.nodeID];
  //UPDATE THE VARIOUS nodeIDs
  leftChildIndex  = node1.nodeID;
  rightChildIndex = node2.nodeID;
  nodeID          = treeNode.size();//this makes it the n+1th node
  
  //UPDATE THE VECTOR OF DATA ITEM INDICES
  itemIndex.insert(itemIndex.end(), node1.itemIndex.begin(), node1.itemIndex.end());
  itemIndex.insert(itemIndex.end(), node2.itemIndex.begin(), node2.itemIndex.end());

  //COMPUTE USEFUL INTERMEDIATE VALUES
  tr1  = log(dirichletProcessParameter) + gammaln(node1.mergePrior + node2.mergePrior);
  tr2  = node1.d_k + node2.d_k;
  a    = max(tr1,tr2);
  b    = min(tr1,tr2);
  d_k  = a + log(1+exp(b-a));
  //UPDATE OVERALL LOG-EVIDENCE ESTIMATE
  a           = max(node1.num1[node2.nodeID], node1.num2[node2.nodeID]);
  b           = min(node1.num1[node2.nodeID], node1.num2[node2.nodeID]);
  logEvidence = a + log(1+exp(b-a));
  //INITIALISE VECTORS TO HOLD MERGER LOG-PROBABILITIES PLUS USEFUL INTERMEDIATES
  num1                = vector<double>(nodeID, 0);
  num2                = vector<double>(nodeID, 0);
  mergeLogProbability = vector<double>(nodeID, -numeric_limits<double>::infinity());
  LogEvidenceForNewCluster = vector<double>(nodeID, 0);
  clusterlengthScale	= vector<double>(nodeID, 0);
  clusterNoiseFree = vector<double>(nodeID, 0);
  clusterNoise = vector<double>(nodeID, 0);
  clusterMixtureComponent = vector<double>(nodeID, 0);
  //CALCULATE THE MERGER PROBABILITIES
#pragma omp parallel for default(shared) schedule(dynamic,1)
  for (i=0; i<nodeID; i++)
  {
    ComputeMergeLogProbability(dataSet, treeNode[i]);
  }
}



//COMPUTE THE MERGER PROBABILITY BETWEEN THIS NODE AND ANOTHER NODE
void Node::ComputeMergeLogProbability(DataSet *dataSet, Node newNode)
{
  //DECLARATIONS
  double tr1, tr2, a, b, ckt, pk, gell;//intermediates for computing weights
  vector<int> mergeIndex=itemIndex;
  vector<int> newIndex  =newNode.itemIndex;
  double newlengthScale, newNoiseFree, newNoise, newMixtureComponent; //the hyperparameters for the cluster
  vector<double> clusterData;
  //CHECK TO SEE IF THIS MERGE IS ALLOWED
  if (mergeFlag==1 || newNode.mergeFlag==1)
    mergeLogProbability[newNode.nodeID] = -numeric_limits<double>::infinity();//if not allowed, skip the calculation and set P(merge)=0
  else
  {
    //ADD THE NEW INDEX TO THE OBJECT
    mergeIndex.insert(mergeIndex.end(), newIndex.begin(), newIndex.end());
    //COMPUTE USEFUL INTERMEDIATE VALUES
    tr1  = log(dirichletProcessParameter) + gammaln(mergePrior + newNode.mergePrior);
    tr2  = d_k + newNode.d_k;
    a    = max(tr1,tr2);
    b    = min(tr1,tr2);
    ckt  = a + log(1+exp(b-a));
    pk   = tr1 - ckt;

    //cout << "mergeIndex:" ;
    //for (i=0; i<mergeIndex.size(); i++){
    //cout << mergeIndex[i] << " ";
    //}
    //cout << endl;
    gell = dataSet->SingleClusterLogEvidence(mergeIndex, newlengthScale, newNoiseFree, newNoise, newMixtureComponent);
    //cout << "gell in Node::ComputeMergeLogProbability:" << gell << endl;
    //UPDATE THE OBJECT TAGS
    num1[newNode.nodeID] = pk + gell;
    num2[newNode.nodeID] = tr2 - ckt + logEvidence + newNode.logEvidence; // backtrace

    mergeLogProbability[newNode.nodeID] = num1[newNode.nodeID] - num2[newNode.nodeID];

    LogEvidenceForNewCluster[newNode.nodeID] = gell;
    clusterlengthScale[newNode.nodeID]	= newlengthScale;
    clusterNoiseFree[newNode.nodeID] = newNoiseFree;
    clusterNoise[newNode.nodeID] = newNoise;
    if (dataSet->GetRobustMode() == 1)
    {
      clusterMixtureComponent[newNode.nodeID] = newMixtureComponent;
    }
  }
}



//METHOD TO WRITE OUT THE RESULTS FROM A VECTOR OF NODES IN SUCH A WAY THAT
//THEY CAN EASILY BE READ INTO R
void Node::OutputResultsToFile(vector<Node> treeNode, string outputFileName)
{
  //DECLARATIONS
  int           i;
  int          nItems;
  ofstream     outputFile(outputFileName.c_str());
  //FIND SOME USEFUL VALUES
  nItems = (treeNode.size() + 1)/2; //there are always n-1 mergers in the BHC algorithm
  //WRITE OUT RESULTS TO FILE
  if (outputFile.is_open())
  {
    for (i=nItems; i<2*nItems-1; i++)
    {
      outputFile << treeNode[i].GetRightChildIndex()+1 << " "
                 << treeNode[i].GetLeftChildIndex()+1  << " "
                 << i-nItems+1                         << " "
                 << treeNode[i].GetClusterLogEvidence() << endl;
    }
    outputFile.close();
  }
  else cout << "Unable to open file";
}



//METHODS TO HANDLE THE MERGER LOGIC
int  Node::AllowedToMerge()
{
  return (mergeFlag==0);
}
void Node::FlagAsMerged()
{
  mergeFlag=1;
}



//'GET' METHODS
double      Node::GetMergeProbability(int index)
{
  return (mergeLogProbability[index]);
}
double      Node::GetClusterLogEvidence()
{
  return (clusterLogEvidence);
}
double      Node::GetGlobalLogEvidence()
{
  return (logEvidence);
}
vector<int> Node::GetItemIndex()
{
  return (itemIndex);
}
int         Node::GetLeftChildIndex()
{
  return (leftChildIndex);
}
int         Node::GetRightChildIndex()
{
  return (rightChildIndex);
}
double	    Node::GetClusterLengthScale(int index)
{
  return (clusterlengthScale[index]);
}
double	    Node::GetClusterNoiseFree(int index)
{
  return (clusterNoiseFree[index]);
}
double	    Node::GetClusterNoise(int index)
{
  return (clusterNoise[index]);
}
double		Node::GetLogEvidenceForNewCluster(int index)
{
  return (LogEvidenceForNewCluster[index]);
}
double		Node::GetClusterMixtureComponent(int index)
{
  return (clusterMixtureComponent[index]);
}



//'SET' METHODS
void Node::SetNodeID(int inputID)
{
  nodeID = inputID;
}
void Node::SetMergeFlag(int inputFlag)
{
  mergeFlag = inputFlag;
}

