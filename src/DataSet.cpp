#include "DataSet.h"

//----------------------------------------------------------------------
// PLACE-HOLDER: READ DATA IN FROM A FILE ------------------------------
//----------------------------------------------------------------------
void DataSet::ReadInData(string dataFile)
{
  cout<<"RID: This method should be replaced by one for the subclass being used!"<<endl;
}

void DataSet::ReadInTimePoints(vector<double> timePoints)
{
  cout<<"RITP: This method should be replaced by one for the subclass being used!"<<endl;
}
//----------------------------------------------------------------------
// PLACE-HOLDER: READ NOISE DATA IN FROM A FILE ------------------------------
//----------------------------------------------------------------------
void DataSet::ReadInNoise(string dataFile)
{
  cout<<"RIN: This method should be replaced by one for the subclass being used!"<<endl;
}
//----------------------------------------------------------------------
// PLACE-HOLDER: READ NOISE DATA IN FROM R INPUTS ------------------------------
//----------------------------------------------------------------------
void DataSet::ReadInNoise(vector<double> noise)
{
  cout<<"RIN: This method should be replaced by one for the subclass being used!"<<endl;
}
//----------------------------------------------------------------------
// PLACE-HOLDER: GET CLUSTER NOISE ------------------------------
//----------------------------------------------------------------------
double DataSet::GetClusterNoise(int NodeID)
{
  cout<<"GCN: This method should be replaced by one for the subclass being used!"<<endl;
  return 0;
}
//----------------------------------------------------------------------
// PLACE-HOLDER: GET MERGED CLUSTER NOISE ------------------------------
//----------------------------------------------------------------------
double DataSet::GetMergedClusterNoise(vector<int>)
{
  cout<<"GMCN: This method should be replaced by one for the subclass being used!"<<endl;
  return 0;
}
//----------------------------------------------------------------------
// DATA INPUT: ACCEPT DATA VECTOR AS INPUT -----------------------------
//----------------------------------------------------------------------
void DataSet::ReadInDataVector(vector<int> inputData, int nDataItems, int nFeatures)
{

  //NEed to fill this in!
  cout<<"RIDV: Need to fill this method in"<<endl;

}

double DataSet::SingleClusterLogEvidence(const vector<int>& itemIndex,
				       double& lengthScale,
				       double& noiseFreeScale,
				       double& noiseSigma,
					 double& mixtureComponent){return 0;}
  double DataSet::GetClusterSE(vector<int> genes) {return 0;}
double DataSet::GetClusterSEMS(vector<int> genes) {return 0;}
void DataSet::SetNoiseMode(int noisemode) {}
void DataSet::SetRobustMode(int robustmode) {}
void DataSet::SetReps(int num_reps) {}
int DataSet::GetRobustMode() {return 0;}
void DataSet::SetDataType(string type) {}
double DataSet::GetMLIINoise(vector<int> itemIndex) {return 0;}
int DataSet::GetNoiseMode() {return 0;}
BlockCovarianceMatrix DataSet::SquareExponentialCovarianceFunctionMissingSingleObservation(vector< vector<double> > KnFC, vector <double> KnC, int KblockSize, int KnRank, int timePoint) {BlockCovarianceMatrix noWarning; return noWarning;}
vector<double> DataSet::GetDataForCluster(vector<int> itemIndex) {vector<double> noWarning; return noWarning;}


//----------------------------------------------------------------------
// DATA INPUT - FIND THE SIZE OF DATA IN A FILE ------------------------
//----------------------------------------------------------------------
void DataSet::FindDataSize(string dataFile)
{
  //DEFINITIONS
  int     featureFlag=1;
  double  inputValue;
  string  line;
  fstream file;
  //INITIALISE THE VARIABLES
  nDataItems = 0;
  nFeatures  = 0;
  //OPEN THE FILE
  file.open (dataFile.c_str(), fstream::in);//dataFile is a string object, so need to convert to a char array here!
  if(file.fail())
  {
    cout<<"Failed to open file "<<dataFile<<"."<<endl;
    if(!system("PAUSE")) /* could not pause */ ;
  }
  //COUNT THE NUMBER OF LINES OF DATA
  while ( getline(file, line) )
  {
    nDataItems++;
    if (featureFlag)
    {
      //count the number of features (only do this once!)
      istringstream iss(line, istringstream::in);//treat the string as an I/O stream
      while (iss >> inputValue)
        nFeatures++;
      featureFlag=0;
    }
  }
  file.close();
}
//----------------------------------------------------------------------
// GET METHODS ---------------------------------------------------------
//----------------------------------------------------------------------
int DataSet::Get_nDataItems()
{
  return (nDataItems);
}
int DataSet::Get_nFeatures()
{
  return (nFeatures);
}
//*****************************************************************************
//*****************************************************************************
//----------------------------------------------------------------------
// ----------------------------------------
//----------------------------------------------------------------------



