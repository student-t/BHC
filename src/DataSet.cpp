/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include "DataSet.h"

/* ----------------------------------------------------------------------
   Find the size of the data in a file.
------------------------------------------------------------------------- */

void DataSet::FindDataSize(string dataFile)
{
  // Definitions
  int featureFlag=1;
  double inputValue;
  string line;
  fstream file;

  // Initialise the variables
  nDataItems = 0;
  nFeatures  = 0;

  // Open the file
  file.open(dataFile.c_str(), fstream::in);
  if(file.fail())
  {
    cout<<"Failed to open file "<<dataFile<<"."<<endl;
    if(!system("PAUSE")) /* could not pause */ ;
  }

  // Count the number of lines of data
  while(getline(file, line))
  {
    nDataItems++;
    if (featureFlag)
    {
      // count the number of features (only do this once!)
      istringstream iss(line, istringstream::in);//treat the string as an I/O stream
      while (iss >> inputValue)
        nFeatures++;
      featureFlag=0;
    }
  }
  file.close();
}

/* ---------------------------------------------------------------------- */

int DataSet::Get_nDataItems()
{
  return nDataItems;
}

/* ---------------------------------------------------------------------- */

int DataSet::Get_nFeatures()
{
  return nFeatures;
}

/* ---------------------------------------------------------------------- */
