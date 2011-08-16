#include "BlockCovarianceMatrix.h"
#include <numeric>

#ifdef ENABLE_LAPACK
#include "LapackWrapper.h"
#endif // ENABLE_LAPACK


// Default constructor
BlockCovarianceMatrix::BlockCovarianceMatrix()
{}



// Factory method. Builds "E", the rank N-1 sub-matrix.
BlockCovarianceMatrix BlockCovarianceMatrix::Build_E_SubMatrix(const BlockCovarianceMatrix& inputMatrix)
{
  //THE OBJET WE ARE GOING TO CONSTRUCT
  BlockCovarianceMatrix bcm;

  //DECLARATIONS
  double noisyScaling, newNoiseFree, temp2_newNoiseFree;
  BlockCovarianceMatrix subMatrix_A;

  //THIS WON'T WORK WITH A RANK-1 INPUT MATRIX
  assert(inputMatrix.nRank>1);

  //ASSIGN SOME OBVIOUS VALUES
  bcm.nRank = inputMatrix.nRank - 1;
  bcm.blockSize = inputMatrix.blockSize;

  //FIND the SUB-MATRIX "A" AND INVERT IT
  subMatrix_A = BlockCovarianceMatrix(inputMatrix.noisyCoeff[0], inputMatrix.noiseFreeCoeff[0][0], bcm.blockSize);
  subMatrix_A.InvertRankOneMatrix();

  //INITIALISE THE ARRAYS
  bcm.noisyCoeff = vector<double>(bcm.nRank);
  bcm.noiseFreeCoeff=vector<vector<double> >(bcm.nRank, vector<double>(bcm.nRank));

  //ASSIGN THE SUB-MATRIX VALUES
  const double temp_newNoiseFree = (subMatrix_A.noisyCoeff[0] + bcm.blockSize)*subMatrix_A.noiseFreeCoeff[0][0]*bcm.blockSize;
  for (int i=0; i<bcm.nRank; i++)
  {
    temp2_newNoiseFree=temp_newNoiseFree*inputMatrix.noiseFreeCoeff[i+1][0]; // B
    for (int j=0; j<bcm.nRank; j++)
    {
      //find adjusting noiseFree term
      newNoiseFree = temp2_newNoiseFree * inputMatrix.noiseFreeCoeff[0][j+1]; // C

      //make adjustments
      bcm.noiseFreeCoeff[i][j] = inputMatrix.noiseFreeCoeff[i+1][j+1] - newNoiseFree;
    }
    //find and adjust the noisey term for this row
    noisyScaling  = inputMatrix.noiseFreeCoeff[i+1][i+1] / bcm.noiseFreeCoeff[i][i];
    bcm.noisyCoeff[i] = inputMatrix.noisyCoeff[i+1] * noisyScaling;
  }

  return bcm;
}



// Constructor. "A", THE RANK 1 SUB-MATRIX.
BlockCovarianceMatrix::BlockCovarianceMatrix(const double noisyValue, const double noiseFreeValue, const int inputBlockSize)
{
  nRank = 1;
  blockSize = inputBlockSize;
  noisyCoeff.push_back(noisyValue);
  noiseFreeCoeff.push_back(vector<double>(1, noiseFreeValue));
}



// Factory method. Builds "E", the rank N-1 sub-matrix, when [A] has a
// single missing observation.
BlockCovarianceMatrix BlockCovarianceMatrix::Build_E_SubMatrixMissingSingleObservation(const BlockCovarianceMatrix& inputMatrix)
{
  //THE OBJECT WE ARE GOING TO CONSTRUCT
  BlockCovarianceMatrix bcm;

  //DECLARATIONS
  double noisyScaling, newNoiseFree, temp2_newNoiseFree;
  BlockCovarianceMatrix subMatrix_A;

  //THIS WON'T WORK WITH A RANK-1 INPUT MATRIX
  assert(inputMatrix.nRank>1);

  //ASSIGN SOME OBVIOUS VALUES
  bcm.nRank     = inputMatrix.nRank - 1;
  bcm.blockSize = inputMatrix.blockSize;

  //FIND THE SUB-MATRIX "A" AND INVERT IT
  subMatrix_A = BlockCovarianceMatrix(inputMatrix.noisyCoeff[0], inputMatrix.noiseFreeCoeff[0][0], bcm.blockSize - 1);
  subMatrix_A.InvertRankOneMatrix();

  //INITIALISE THE ARRAYS
  bcm.noisyCoeff = vector<double>(bcm.nRank);
  bcm.noiseFreeCoeff = vector<vector<double> >(bcm.nRank, vector<double>(bcm.nRank));

  //ASSIGN THE SUB-MATRIX VALUES
  //a pre-computation
  const double temp_newNoiseFree = (subMatrix_A.noisyCoeff[0]+(bcm.blockSize - 1))*subMatrix_A.noiseFreeCoeff[0][0]*(bcm.blockSize - 1);
  for (int i=0; i<bcm.nRank; i++)
  {
    temp2_newNoiseFree = temp_newNoiseFree * inputMatrix.noiseFreeCoeff[i+1][0]; // B

    for (int j=0; j<bcm.nRank; j++)
    {
      //find adjusting noiseFree term
      newNoiseFree = temp2_newNoiseFree * inputMatrix.noiseFreeCoeff[0][j+1]; // C

      //make the adjustments
      bcm.noiseFreeCoeff[i][j] = inputMatrix.noiseFreeCoeff[i+1][j+1] - newNoiseFree;
    }
    //find and adjust the noisey term for this row
    noisyScaling  = inputMatrix.noiseFreeCoeff[i+1][i+1] / bcm.noiseFreeCoeff[i][i];
    bcm.noisyCoeff[i] = inputMatrix.noisyCoeff[i+1] * noisyScaling;
  }

  return bcm;
}



// Perform a recursive matrix inversion.
void BlockCovarianceMatrix::InvertMatrix()
{
  if(nRank==1)
  {
    InvertRankOneMatrix();
  }
  else
  {
#ifdef ENABLE_LAPACK
    // This case currently does not work, so ENABLE_LAPACK is not defined
    if(blockSize==1)
    {
      // Use LAPACK to perform general matrix inversion
      InvertGeneralMatrix();
    }
    else
#endif // ENABLE_LAPACK
    {
      InvertBlockMatrix();
    }
  }
}



// Use LAPACK to perform a general matrix inversion.
//
// N.B. This function seems quite inefficient; the use of a
// 'workingMatrix' seems superfluous. However, since the code
// is currently not using this function - and I'm lazy - I shall leave it.
void BlockCovarianceMatrix::InvertGeneralMatrix()
{
#ifdef ENABLE_LAPACK
  //DECLARATIONS
  int    i, j;
  double currentElement;
  double *workingMatrix;

  //ASSERT THAT THIS ISN'T A BLOCK MATRIX
  assert(blockSize==1);
  assert(nRank>0);

  //MEMORY ALLOCATION
  workingMatrix=new double[nRank*nRank];

  //CONSTRUCT A VERSION OF THE COVARIANCE MATRIX (FOR LAPACK)
  for (i=0; i<nRank; i++)
  {
    for (j=1; j<nRank; j++)
    {
      currentElement = noiseFreeCoeff[i][j];
      if (i==j) currentElement *= noisyCoeff[i] + 1;
      workingMatrix[i*nRank+j] = currentElement;
    }
  }

  //PERFORM THE MATRIX INVERSION USING LAPACK
  LAPACK_MatrixInverse(workingMatrix, nRank);

  //COPY THE INVERTED MATRIX BACK INTO THIS OBJECT
  for (i=0; i<nRank; i++)
  {
    for (j=1; j<nRank; j++)
    {
      //copy all the noise-free values
      noiseFreeCoeff[i][j] = workingMatrix[i*nRank+j];
    }

    //make adjustments for the diagonal, noisy terms
    noisyCoeff[i] = workingMatrix[i*nRank+j] - 1;
    noiseFreeCoeff[i][i] = 1;
  }

  // Cleanup
  delete [] workingMatrix;
#endif // ENABLE_LAPACK
}



// Use the identity for sub-dividing a block matrix:
//   K = [A B]  and  E = D - C*A_inverse*B
//       [C D]
void BlockCovarianceMatrix::InvertBlockMatrix()
{
  //DECLARATIONS
  int i, j;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;
  vector<double>        noiseFree_B, offDiagonal;
  double                diagonal_noisy, diagonal_noiseFree, factor_AB;

  //IF WE REACH HERE, USE THE BLOCK-MATRIX STRUCTURE TO OUR ADVANTAGE
  //FIND the SUB-MATRICES "A" AND "E"
  subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
  subMatrix_E = Build_E_SubMatrix(*this);

  //INVERT "A" AND "E"
  subMatrix_A.InvertRankOneMatrix();
  subMatrix_E.InvertMatrix();

  //FIND THE SUB-MATRICES B AND C
  //symmetry => only need one of these
  noiseFree_B = noiseFreeCoeff[0];
  noiseFree_B.erase(noiseFree_B.begin(), noiseFree_B.begin()+1);//remove the first element as this is part of A, not B

  //FIND THE FACTOR COMING FROM A_inv * B
  factor_AB  = subMatrix_A.noiseFreeCoeff[0][0];
  factor_AB *= blockSize + subMatrix_A.noisyCoeff[0];

  //HENCE CONSTRUCT THE OVERALL INVERSE MATRIX
  //copy E_inverse into the D slots
  for (i=1; i<nRank; i++)
  {
    noisyCoeff[i] = subMatrix_E.noisyCoeff[i-1];
    for (j=1; j<nRank; j++)
    {
      noiseFreeCoeff[i][j] = subMatrix_E.noiseFreeCoeff[i-1][j-1];
    }
  }
  //construct the remaining off-diagonal elements, using B * E_inv
  offDiagonal = subMatrix_E.BlockMultiply(noiseFree_B);
  for (i=1; i<nRank; i++)
  {
    //using the fact our matrix is symmetric
    noiseFreeCoeff[0][i] = noiseFreeCoeff[i][0] = -factor_AB * offDiagonal[i-1];
  }
  //construct the final diagonal element
  diagonal_noiseFree = 0;
  for (i=0; i<(nRank-1); i++) //need to find B * E_inv * C here...
  {
    //using the fact that C_transpose = B
    diagonal_noiseFree += blockSize * offDiagonal[i] * noiseFree_B[i];
  }
  diagonal_noiseFree *= factor_AB * factor_AB;
  //also need to add on A_inverse
  diagonal_noiseFree += subMatrix_A.noiseFreeCoeff[0][0];
  //hence compute the new noisy term
  diagonal_noisy = subMatrix_A.noisyCoeff[0] * subMatrix_A.noiseFreeCoeff[0][0];
  diagonal_noisy /= diagonal_noiseFree;
  //and store the values in this object
  noiseFreeCoeff[0][0] = diagonal_noiseFree;
  noisyCoeff[0] = diagonal_noisy;
}



// Invert analytically a rank-one matrix.
void BlockCovarianceMatrix::InvertRankOneMatrix()
{
  //DECLARATIONS
  double newNoisy, newNoiseFree;

  //COMPUTE THE NEW NOISY VALUE (FOR THE DIAGONAL ELEMENTS)
  newNoisy  = -noisyCoeff[0] - blockSize;

  //COMPUTE THE NEW NOISE-FREE VALUE
  newNoiseFree  = -1.0 /
    (noiseFreeCoeff[0][0] * (noisyCoeff[0] * (noisyCoeff[0] + blockSize)));

  //UPDATE THE MATRIX VALUES, SO THIS MATRIX OBJECT IS NOW INVERTED
  noisyCoeff[0] = newNoisy;
  noiseFreeCoeff[0][0] = newNoiseFree;
}



// Compute the matrix log-determinant for the matrix missing a single observation.
double BlockCovarianceMatrix::ComputeMatrixMissingSingleObservationLogDeterminant() const
{
  //DECLARATIONS
  double logDeterminant;
  BlockCovarianceMatrix subMatrix_A, subMatrix_Efirst; //use the maths notation here

  //extract small A and first sub E matrix
  subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize-1);
  subMatrix_Efirst = Build_E_SubMatrixMissingSingleObservation(*this);

  //compute the contributions from the sub-matrices
  logDeterminant  = subMatrix_A.ComputeRankOneMatrixDeterminant();
  logDeterminant += subMatrix_Efirst.ComputeMatrixDeterminant();

  return logDeterminant;
}



// Compute the matrix inverse for the matrix missing a single observation.
// (The matrix loses its first row and column.)
void BlockCovarianceMatrix::InvertMatrixMissingSingleObservation()
{
  //DECLARATIONS
  int i, j;
  BlockCovarianceMatrix subMatrix_A1, subMatrix_E1;
  vector<double>        noiseFree_B, offDiagonal;
  double                diagonal_noisy, diagonal_noiseFree, factor_AB;

  //M = [A,B][C,D], A is a smaller rank 1 matrix, since missing observation

  //RECURSION TO FIND THE INVERSE MATRIX
  //FIND the SUB-MATRICES "A" AND "E"
  subMatrix_A1 = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize-1); //A is missing an observation, so has a smaller blockSize
  subMatrix_E1 = Build_E_SubMatrixMissingSingleObservation(*this);

  //INVERT "A" AND "E"
  subMatrix_A1.InvertRankOneMatrix();
  subMatrix_E1.InvertMatrix();

  //FIND THE SUB-MATRICES B AND C
  //symmetry => only need one of these
  noiseFree_B = noiseFreeCoeff[0]; //the first row of the M.noiseFreeCoeff
  noiseFree_B.erase(noiseFree_B.begin(), noiseFree_B.begin()+1);//remove the first element as this is part of A, not B

  //FIND THE FACTOR COMING FROM A_inv * B
  factor_AB  = subMatrix_A1.noiseFreeCoeff[0][0];
  factor_AB *= (blockSize - 1) + subMatrix_A1.noisyCoeff[0]; //not included any values from B matrix yet

  //HENCE CONSTRUCT THE OVERALL INVERSE MATRIX
  //copy E_inverse into the D slots
  for (i=1; i<nRank; i++)
  {
    noisyCoeff[i] = subMatrix_E1.noisyCoeff[i-1]; //copying the diagonal elements
    for (j=1; j<nRank; j++)
    {
      noiseFreeCoeff[i][j] = subMatrix_E1.noiseFreeCoeff[i-1][j-1];
    }
  }
  //construct the remaining off-diagonal elements, using B * E_inv
  offDiagonal = subMatrix_E1.BlockMultiply(noiseFree_B);
  for (i=1; i<nRank; i++)
  {
    //using the fact our matrix is symmetric
    noiseFreeCoeff[0][i] = noiseFreeCoeff[i][0] = -factor_AB * offDiagonal[i-1];
  }
  //construct the final diagonal element
  diagonal_noiseFree = 0;
  for (i=0; i<(nRank-1); i++) //need to find B * E_inv * C here...
  {
    //using the fact that C_transpose = B, dont subtract 1 from blockSize
    //size going across rows
    diagonal_noiseFree += blockSize * offDiagonal[i] * noiseFree_B[i];
  }
  diagonal_noiseFree *= factor_AB * factor_AB; // as is
  //also need to add on A_inverse
  diagonal_noiseFree += subMatrix_A1.noiseFreeCoeff[0][0]; // as is
  //hence compute the new noisy term
  diagonal_noisy = subMatrix_A1.noisyCoeff[0] * subMatrix_A1.noiseFreeCoeff[0][0];
  diagonal_noisy /= diagonal_noiseFree;
  //and store the values in this object
  noiseFreeCoeff[0][0] = diagonal_noiseFree;
  noisyCoeff[0] = diagonal_noisy;
}



// Compute the matrix log determinant
// Using det[A B][C D] = det(A)det(E) where E = D-C*inv(A)*B
double BlockCovarianceMatrix::ComputeMatrixDeterminant() const
{
  //DECLARATIONS
  double logDeterminant;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;//use the maths notation here

  //RECURSION TO FIND THE LOG-DET
  if (nRank==1)
  {
    logDeterminant = ComputeRankOneMatrixDeterminant();
  }
  else
  {
    //extract rank n-1 and rank-1 sub-matrices
    //remember that we want "E"
    subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
    subMatrix_E = Build_E_SubMatrix(*this);

    //compute the contributions from the sub-matrices
    logDeterminant  = subMatrix_A.ComputeRankOneMatrixDeterminant();
    logDeterminant += subMatrix_E.ComputeMatrixDeterminant();
  }

  return logDeterminant;
}



// Compute analytically the rank-one matrix determinant.
double BlockCovarianceMatrix::ComputeRankOneMatrixDeterminant() const
{
  //DECLARATIONS
  double logDeterminant;

  //CHECK HERE THAT THIS DOUBLELY IS A RANK-ONE MATRIX
  assert(nRank==1);

  if (noiseFreeCoeff[0][0] < 0.0 || noisyCoeff[0] < 0.0)
  {
    cout << "problem with log Determinant: nan" << endl;
  }

  //COMPUTE THE LOG-DETERMINANT
  logDeterminant  = log(noiseFreeCoeff[0][0]) * blockSize;
  logDeterminant += log(noisyCoeff[0]) * (blockSize - 1);
  logDeterminant += log(noisyCoeff[0] + blockSize);

  return logDeterminant;
}



// Method to block-multiply a block-matrix and a vector
vector<double> BlockCovarianceMatrix::BlockMultiply(const vector<double>& inputVector) const
{
  //DECLARATIONS
  int i;
  double currentElement;
  const int input_size=inputVector.size();
  vector<double> outputVector=vector<double>(input_size);

  //COMPUTE THE OUTPUT VECTOR ELEMENTS
  for (i=0; i<input_size; i++)
  {
    // Inner product: <noiseFreeCoeff row, inputVector>
    currentElement = inner_product(inputVector.begin(),
                                   inputVector.end(),
                                   noiseFreeCoeff[i].begin(),
                                   0.0);

    // Since it's a block matrix, we use each matrix term blockSize-many times:
    currentElement *= blockSize;

    // Add the (diagonal) noise term:
    currentElement += noisyCoeff[i] * noiseFreeCoeff[i][i] * inputVector[i];

    // Store the result
    outputVector[i] = currentElement;
  }

  return outputVector;
}



// Method to find the (i,j)-th element of the block-matrix
double BlockCovarianceMatrix::GetElement(const int i, const int j) const
{
  //DECLARATIONS
  int block_i, block_j;
  double matrixElement;

  //FIND WHICH BLOCK THE ELEMENT BELONGS TO
  block_i = i / blockSize;
  block_j = j / blockSize;

  //FIND THE MATRIX ELEMENT
  matrixElement = noiseFreeCoeff[block_i][block_j];

  //IF IT'S A DIAGONAL ELEMENT, MAKE AN ADJUSTMENT
  if (i==j)
    matrixElement *= 1.0 + noisyCoeff[block_i];

  return matrixElement;
}



// Method to find the i-th row of the block-matrix when a single observation is missing
vector<double> BlockCovarianceMatrix::GetRowForMatrixMissingSingleObservation(const int index) const
{
  //DECLARATIONS
  int newindex, block_i;
  vector<double> rowVector = vector<double>(nRank*blockSize);

  //Make alteration to the index
  newindex = index + 1;

  //FIND WHICH BLOCK THE ROW BELONGS TO
  block_i = newindex / blockSize;

  //FIND THE MATRIX ELEMENTS (noise-free component)
  vector<double>::const_iterator inIt = noiseFreeCoeff[block_i].begin();
  vector<double>::iterator rowVectorIt = rowVector.begin();
  while (inIt != noiseFreeCoeff[block_i].end())
  {
    vector<double>::iterator outItEnd = rowVectorIt + blockSize;
    fill(rowVectorIt, outItEnd, *inIt++);
    rowVectorIt = outItEnd;
  }
  // ^ this is a much more efficient version of the following:
  /*
  for (i=0; i<nRank; i++)
  {
    for (j=0; j<blockSize; j++)
    {
      rowVector.push_back(noiseFreeCoeff[block_i][i]);
    }
  }
  */

  //ADJUST THE DIAGONAL ELEMENT
  rowVector[newindex] *= 1 + noisyCoeff[block_i];
  rowVector.erase(rowVector.begin(), rowVector.begin()+1);//Delete the first element of rowVector, since this corresponds to the missing observation

  return rowVector;
}



// Method to build the full i-th row of the block-matrix
vector<double> BlockCovarianceMatrix::GetRow(const int index) const
{
  //DECLARATIONS
  vector<double> rowVector=vector<double>(nRank*blockSize);
  
  //FIND WHICH BLOCK THE ROW BELONGS TO
  int block_i = index / blockSize; //index goes down the rows in the matrix

  //FIND THE MATRIX ELEMENTS (noise-free component)
  vector<double>::const_iterator inIt = noiseFreeCoeff[block_i].begin();
  vector<double>::iterator rowVectorIt = rowVector.begin();
  while (inIt != noiseFreeCoeff[block_i].end())
    {
      vector<double>::iterator outItEnd = rowVectorIt + blockSize;
      fill(rowVectorIt, outItEnd, *inIt++);
      rowVectorIt = outItEnd;
    }
  // ^ this is a more efficient version of the following:
  /*
    for (i=0; i<nRank; i++)
    {
    for (j=0; j<blockSize; j++)
    {
    rowVector.push_back(noiseFreeCoeff[block_i][i]);
    }
    }
  */

  //ADJUST THE DIAGONAL ELEMENT
  rowVector[index] *= 1.0 + noisyCoeff[block_i];

  return rowVector;
}



// Method to form the product y_transpose * K * y for K missing a single observation
// N.B. This assumes that we have already inverted K.
double BlockCovarianceMatrix::ComputeLogLikelihoodProductMissingSingleObservation(const vector<double>& data) const
{
  //DECLARATIONS
  int k;
  double logLike, currentValue;
  vector<double> currentRow;
  vector<double>::const_iterator dataIt=data.begin();

  //ASSERT THAT THE data VECTOR IS THE CORRECT SIZE
  assert(data_size==(nRank*blockSize-1)); //Should be missing an observation

  //FORM FIRST PRODUCT, K*y
  logLike=0;
  //going across the rows in the K matrix
  for(k=0; dataIt != data.end(); ++dataIt, ++k)
  {
    //get the current row from the covariance matrix
    currentRow = GetRowForMatrixMissingSingleObservation(k);
    //loop to sum the terms in a vector*matrix product
    currentValue = std::inner_product(data.begin(), data.end(), currentRow.begin(), 0.0);
    //now form the contribution to the overall product, y_t * K * y
    logLike += (*dataIt) * currentValue;
  }

  return logLike;
}



// Method to compute the bilinear form, y_transpose * K * y, where y=data.
// In practice, K should have been inverted before calling this.
double BlockCovarianceMatrix::ComputeLogLikelihoodProduct(const vector<double>& data) const
{
  // Here is the simple and concise way of doing it:
  // Compute K * y
  vector<double> temp=VectorMultiply(data);
  // Compute and return y^t * (K * y)
  return std::inner_product(temp.begin(), temp.end(), data.begin(), 0.0);

  //////////////////////////////////////////////////////////////////////

  // Here is the long and convoluted way which, in theory, should be faster,
  // but isn't!
  /*
  vector<double>::const_iterator dataIt_e;
  vector<double> blockData = vector<double>(nRank);
  vector<double>::iterator blockDataIt;
  vector<double>::const_iterator dataIt;
  double ans=0;
 

  // Build a block-vector version of the data
  blockDataIt = blockData.begin();
  dataIt = data.begin();
  while(blockDataIt != blockData.end())
    {
      dataIt_e = dataIt + blockSize;
      *blockDataIt++ = std::accumulate(dataIt, dataIt_e, 0.0);
      dataIt = dataIt_e;
    }
  
  // Compute blockData^t * noiseFreeCoeff * blockData
  vector<vector<double> >::const_iterator noiseFreeCoeffIt = noiseFreeCoeff.begin();
  int id=0;
  while(noiseFreeCoeffIt != noiseFreeCoeff.end())
    {
      vector<double> row = *noiseFreeCoeffIt++;
      ans += std::inner_product(row.begin(), row.end(), blockData.begin(), 0.0)*blockData[id];
      id++;
    }

  // Noise correction
  double noise=0;
  for(int i=0; i<nRank; i++)
    {
      double ip=0;
      for(int j=i*blockSize; j<(i+1)*blockSize; j++)
	{
	  ip += data[j]*data[j];
	}
      noise += noiseFreeCoeff[i][i] * noisyCoeff[i] * ip;
    }
  
  return ans+noise;*/
}



// Compute the product K * y, where y=inputVector.
vector<double> BlockCovarianceMatrix::VectorMultiply(const vector<double>& inputVector) const
{
  //DECLARATIONS
  vector<double> outputVector=vector<double>(inputVector.size());
  vector<double> row;
  vector<double>::const_iterator inIt, inIt_e;
  int i,j, out_id,out_id_e;
  double dot, noisyFactor, temp;

  out_id=0;
  for(i=0; i<nRank; i++)
    {
      row = noiseFreeCoeff[i];
      noisyFactor = noisyCoeff[i]*row[i];
      inIt = inputVector.begin();
      dot=0;
      
      // Iterate over each block in this row
      for(j=0; j<nRank; j++)
	{
	  // Apply the coeff to blockSize-many elements from the input vector
	  inIt_e = inIt+blockSize;
	  temp = std::accumulate(inIt, inIt_e, 0.0);
	  dot += temp*row[j];

	  // Move our iterator to the next block of our vector
	  inIt = inIt_e;
	}
      
      // Store the result
      out_id_e = out_id + blockSize;
      for(; out_id<out_id_e; out_id++)
	{
	  outputVector[out_id] = dot + inputVector[out_id]*noisyFactor;
	}
    }
  
  return outputVector;
}



/* I have made no changes to this function - Robert
//METHOD TO FORM THE PRODUCT y_transpose * K * y
//(note that this assumes we've already inverted K)
double BlockCovarianceMatrix::ComputeLogLikelihoodProduct(vector<double> data){
  //DECLARATIONS
  int            i, j, k, counter;
  double         logLike=0, currentValue, currentElement;
  vector<double> workingVector;
  //ASSERT THAT THE data VECTOR IS THE CORRECT SIZE
  assert(data.size()==nRanks*blockSize);
  //FORM FIRST PRODUCT, K*y
  for (k=0; k<data.size(); k++){
    currentValue = 0; //initialise for the new element
    counter      = 0;
    //loop to sum the terms in a vector*matrix product
    for (i=0; i<nRank; i++){
      for (j=0; j<blockSize; j++){
	currentValue += data[counter] * GetElement(k, counter);
	counter++;
      }
    }
    //store the completed value
    workingVector.push_back(currentValue);
  }
  //FORM THE OVERALL PRODUCT, y_t * K * y
  for (i=0; i<data.size(); i++)
    logLike += data[i] * workingVector[i];

  return(logLike);
}
*/
