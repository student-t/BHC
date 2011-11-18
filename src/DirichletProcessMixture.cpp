/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <limits>

#include "DirichletProcessMixture.h"

/* ----------------------------------------------------------------------
   Default constructor.
---------------------------------------------------------------------- */
   
DirichletProcessMixture::DirichletProcessMixture() {}

/* ----------------------------------------------------------------------
   Apply randomised-BHC to the provided data set.

   Arguments:
      m - the dimension of the randomly selected subset D_m.
      t - when a partition becomes less than max(m+1,t) in size,
          greedy clustering is applied instead of another randomised
	  recursion. If you're not sure what value to give it then
	  t=0 is perfectly reasonable.
---------------------------------------------------------------------- */

vector<Node> DirichletProcessMixture::
RandomisedClustering(DataSet& dataSet,
		     const int m,
		     const int t,
		     const bool verbose)
{
  // Declarations
  int i,depth;
  int n=dataSet.nDataItems;
  vector<Node> D, D_clone;

  // Build our starting tree of nodes
  D=vector<Node>(n);
  D_clone=vector<Node>(n);

  for(i=0; i<n; i++)
    {
      D[i]=Node::CreateDataNode(dataSet, i);
      D_clone[i]=D[i]; // copy
    }

  // Resize D to accommodate the entire tree
  D.resize(2*n-1);
  
  // Seed the pseudo-random number generator
  srand(time(NULL));

  // Begin the recursion
  depth=2*n-2;
  RandomisedClustering(D_clone, D, dataSet, m, t, depth, verbose);

  // Dump the final clusterings
  if(verbose)
    {
      // w is the max width (number of digits) for the node indices,
      // this is used to keep the output below tidy and aligned
      int mult=10;
      int w=1;
      while(2*n-1 >= mult)
	{
	  ++w;
	  mult*=10;
	}
      cout << "\nClus ID\tMerge IDs\tLogEv\n"
	   << "-------\t---------\t-----"
	   << endl;
      cout.precision(4);
      for(i=n; i<int(D.size()); i++)
	{
	  cout << i+1 << "\t(" << std::setw(w) << D[i].GetLeftChildNodeID()+1 << ", "
	       << std::setw(w) << D[i].GetRightChildNodeID()+1 << ")\t"
	       << fixed << D[i].GetClusterLogEvidence() << endl;
	}
      cout<<"\n(Lower bound) model Log Ev: "<<D.back().GetLowerBoundLogEvidence()<<endl;
    }
  
  return D;
}

/* ----------------------------------------------------------------------
   Applies RBHC to D_sub. Creates a node of index 'depth' in D which is
   the top-level node of the tree of D_sub. All deeper node merges are
   created at depth-1, depth-2, etc. The value of depth is updated to
   the index that the next recursion should be passed.
---------------------------------------------------------------------- */

void DirichletProcessMixture::
RandomisedClustering(vector<Node>& D_sub,
		     vector<Node>& D,
		     DataSet& dataSet,
		     const int m,
		     const int t,
		     int& depth,
		     const bool verbose)
{
  int n=D_sub.size();
  int depth0=depth; // root depth
  int i,j,id,leftChildNodeID,rightChildNodeID;
  vector<Node> D_m, D_L, D_R;
  vector<int> random_id;
  Node toplevel_L,toplevel_R;
  vector<int> toplevel_L_childDataIDs,toplevel_R_childDataIDs;
  bool flag;
  double prob_L, prob_R;
  Node tempMerge;
  const static int LEFT=0;
  const static int RIGHT=1;

  if(n==2)
    {
      // Just merge a single pair of data nodes
      D[depth]=Node::CreateMergerNode(dataSet,
				      D[D_sub[0].GetDataID()],
				      D[D_sub[1].GetDataID()],
				      depth);
      --depth;
    }else if(n<=max(m, t-1))
    {
      // Adjust labels of D_sub
      for(i=0; i<n; i++) D_sub[i].SetNodeID(i);
      
      // Perform the greedy clustering (BHC) algorithm on D_sub
      GreedyClustering(D_sub, dataSet, false);

      // Map the tree in D_sub onto the tree in D
      for(i=2*n-2; i>=n; i--)
	{
	  D[depth]=D_sub[i];
	  D[depth].SetNodeID(depth);

	  // Update the child indices
	  leftChildNodeID=D[depth].GetLeftChildNodeID();
	  rightChildNodeID=D[depth].GetRightChildNodeID();
	  if(leftChildNodeID<n)
	    D[depth].SetLeftChildNodeID(D_sub[leftChildNodeID].GetDataID());
	  else
	    D[depth].SetLeftChildNodeID(depth0 + leftChildNodeID - 2*n + 2);
	  if(rightChildNodeID<n)
	    D[depth].SetRightChildNodeID(D_sub[rightChildNodeID].GetDataID());
	  else
	    D[depth].SetRightChildNodeID(depth0 + rightChildNodeID - 2*n + 2);

	  --depth;
	}
    }else
    {
      /// Construct the randomised m-dimensional subset D_m
      D_m=vector<Node>(m);
      random_id=vector<int>(n);
      for(i=0; i<n; i++) random_id[i]=i;
      std::random_shuffle(random_id.begin(), random_id.end());
      for(i=0; i<m; i++)
	{
	  id=random_id[i];
	  D_m[i]=D_sub[id];
	  D_m[i].SetNodeID(i);
	}

      /// Perform the greedy clustering on D_m
      GreedyClustering(D_m, dataSet, false);

      /// Filter D_sub through tree(D_m) to get partitions D_sub=D_L union D_R
      toplevel_L=D_m[D_m.back().GetLeftChildNodeID()];
      toplevel_R=D_m[D_m.back().GetRightChildNodeID()];
      // Explicitly partition the nodes of D_m first ...
      toplevel_L_childDataIDs=toplevel_L.GetChildDataIDs();
      toplevel_R_childDataIDs=toplevel_R.GetChildDataIDs();
      for(i=0; i<int(toplevel_L_childDataIDs.size()); i++)
	D_L.push_back(D[toplevel_L_childDataIDs[i]]);
      for(i=0; i<int(toplevel_R_childDataIDs.size()); i++)
	D_R.push_back(D[toplevel_R_childDataIDs[i]]);
      // ... and now the rest, D_sub\D_m
      vector<int> apartheid=vector<int>(n);
#pragma omp parallel for default(shared) private(i) schedule(dynamic,1)
      for(i=0; i<n; i++)
	{
	  // Book-keeping
	  apartheid[i]=-1;
	  
	  // Exclude elements from D_m since we have already dealt with them
	  flag=false;
	  for(j=0; j<m; j++) if(i==random_id[j]) {flag=true; break;}
	  if(flag) continue;

	  // Does D_sub[i] want to go in D_L or D_R?
	  tempMerge=Node::CreateMergerNode(dataSet,toplevel_L,D_sub[i]);
	  prob_L=tempMerge.GetClusterLogEvidence();
	  tempMerge=Node::CreateMergerNode(dataSet,toplevel_R,D_sub[i]);
	  prob_R=tempMerge.GetClusterLogEvidence();
	  if(prob_L > prob_R)
	    apartheid[i]=LEFT;
	  else
	    apartheid[i]=RIGHT;
	}
      // Carry out the segregation (in serial)
      for(i=0; i<n; i++)
	{
	  if(apartheid[i]==LEFT)
	    D_L.push_back(D_sub[i]);
	  else if(apartheid[i]==RIGHT)
	    D_R.push_back(D_sub[i]);
	}
      
      // Adjust node ID's
      for(i=0; i<int(D_L.size()); i++) D_L[i].SetNodeID(i);
      for(i=0; i<int(D_R.size()); i++) D_R[i].SetNodeID(i);

      if(verbose)
	{
	  cout << "Partitioned (" << n << ") into L=(" << D_L.size()
	       << ") and R=(" << D_R.size() << ")"
	       << endl;
	}
      
      /// Recursion: RBHC(D_L) and RBHC(D_R)
      if(D_L.size()==1)
	{
	  --depth;
	  RandomisedClustering(D_R, D, dataSet, m, t, depth, verbose);
	  D[depth0]=Node::CreateMergerNode(dataSet,
					   D[depth0-1],
					   D[D_L[0].GetDataID()],
					   depth0);
	}else if(D_R.size()==1)
	{
	  --depth;
	  RandomisedClustering(D_L, D, dataSet, m, t, depth, verbose);
	  D[depth0]=Node::CreateMergerNode(dataSet,
					   D[depth0-1],
					   D[D_R[0].GetDataID()],
					   depth0);
	}else
	{
	  --depth;
	  RandomisedClustering(D_L, D, dataSet, m, t, depth, verbose);
	  int temp=depth;
	  RandomisedClustering(D_R, D, dataSet, m, t, depth, verbose);
	  D[depth0]=Node::CreateMergerNode(dataSet,
					   D[depth0-1],
					   D[temp],
					   depth0);
	}
    }
}

/* ----------------------------------------------------------------------
   Run the greedy clustering algorithm on the provided data set.
---------------------------------------------------------------------- */

vector<Node> DirichletProcessMixture::
GreedyClustering(DataSet& dataSet,
		 const bool verbose)
{
  int n=dataSet.nDataItems;
  int i;
  vector<Node> tree=vector<Node>(n);

  // Build our starting tree of nodes
#pragma omp parallel for default(shared) private(i) schedule(dynamic,1)
  for(i=0; i<n; i++) tree[i]=Node::CreateDataNode(dataSet, i);

  // Cluster the tree
  GreedyClustering(tree, dataSet, verbose);
  
  return tree;
}

/* ----------------------------------------------------------------------
   This is where the greedy clustering work is actually done. It acts
   on a vector of nodes; this allows code reuse by the randomised method.
---------------------------------------------------------------------- */

void DirichletProcessMixture::
GreedyClustering(vector<Node>& tree,
		 DataSet& dataSet,
		 const bool verbose)
{
  int n=tree.size();//=dataSet.nDataItems
  vector<Node> trialTree=vector<Node>(n*(n-1)/2);
  int i,j;
  double maxClusterLogEvidence, tempClusterLogEvidence;
  int maxClusterID=-1;
  int w=1;

  // Perform all possible pairwise merges; O(n^2)
#pragma omp parallel for default(shared) private(i,j) schedule(dynamic,1)
  for(i=1; i<n; i++)
    {
      int offset=(i-1)*i/2;
      for(j=0; j<i; j++)
	{
	  // Pair nodes i and j
	  trialTree[offset+j]=Node::CreateMergerNode(dataSet,tree[i],tree[j]);
	}
    }

  // Resize the tree to accommodate all final merges
  tree.resize(2*n-1);

  if(verbose)
    {
      // w is the max width (number of digits) for the node indices,
      // this is used to keep the output below tidy and aligned
      int mult=10;
      while(2*n-1 >= mult)
	{
	  ++w;
	  mult*=10;
	}
      cout << "Clus ID\tMerge IDs\tLog Ev" << endl;
      cout << "-------\t---------\t------" << endl;
      cout.precision(4);
    }
  
  // Now we construct our tree; there are n-1 merges to make
  for(i=n; i<2*n-1; i++)
    {
      // Find the most likely merger
      maxClusterLogEvidence=-numeric_limits<double>::infinity();
      maxClusterID=-1;
      for(j=0; j<int(trialTree.size()); j++)
	{
	  if(tree[trialTree[j].GetLeftChildNodeID()].AllowedToMerge() &&
	     tree[trialTree[j].GetRightChildNodeID()].AllowedToMerge())
	    {
	      tempClusterLogEvidence=trialTree[j].GetClusterLogEvidence();
	      if(tempClusterLogEvidence>=maxClusterLogEvidence)
		{
		  maxClusterLogEvidence=tempClusterLogEvidence;
		  maxClusterID=j;
		}
	    }
	}

      if(verbose)
	{
	  cout << i+1 << "\t("
	       << std::setw(w)
	       << trialTree[maxClusterID].GetLeftChildNodeID()+1
	       << ", "
	       << std::setw(w)
	       << trialTree[maxClusterID].GetRightChildNodeID()+1
	       << ")\t"
	       << fixed
	       << maxClusterLogEvidence
	       << endl;
	}
      
      // Add this merged node to our tree
      tree[i]=trialTree[maxClusterID]; // this creates a copy
      tree[i].SetNodeID(i);

      // Make sure these two nodes are not merged again later
      tree[tree[i].GetLeftChildNodeID()].FlagAsMerged();
      tree[tree[i].GetRightChildNodeID()].FlagAsMerged();

      // Now we extend our trialTree to incorporate this new node
      vector<int> ids;
      for(j=0; j<i; j++) if(tree[j].AllowedToMerge()) ids.push_back(j);
      int offset=trialTree.size();
      trialTree.resize(offset+ids.size());
#pragma omp parallel for default(shared) private(j) schedule(dynamic,1)
      for(j=0; j<int(ids.size()); j++)
	{
	  trialTree[offset+j]=Node::CreateMergerNode(dataSet,tree[ids[j]],tree[i]);
	}
    }

  // Quick cleanup
  for(i=0; i<int(tree.size()); i++) tree[i].UnflagAsMerged();

  if(verbose)
    {
      cout << "\n(Lower bound) model Log Ev: "
	   << fixed << tree.back().GetLowerBoundLogEvidence()
	   << endl;
    }
}

