/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <numeric>

#include "BlockCovarianceMatrix.h"

/* ----------------------------------------------------------------------
   Default constructor
---------------------------------------------------------------------- */

BlockCovarianceMatrix::BlockCovarianceMatrix() {}

/* ----------------------------------------------------------------------
   Constructor: builds "A", the rank-1 sub-matrix
---------------------------------------------------------------------- */

BlockCovarianceMatrix::
BlockCovarianceMatrix(const double noisyValue,
		      const double noiseFreeValue,
		      const int inputBlockSize)
{
  nRank = 1;
  blockSize = inputBlockSize;
  noisyCoeff.push_back(noisyValue);
  noiseFreeCoeff.push_back(vector<double>(1, noiseFreeValue));
}

/* ----------------------------------------------------------------------
   Factory method: builds "E", the rank N-1 sub-matrix
---------------------------------------------------------------------- */

BlockCovarianceMatrix BlockCovarianceMatrix::
Build_E_SubMatrix(const BlockCovarianceMatrix& inputMatrix)
{
  // The object we are going to construct
  BlockCovarianceMatrix bcm;

  // Declarations
  double noisyScaling, newNoiseFree, temp2_newNoiseFree;
  BlockCovarianceMatrix subMatrix_A;

  // This won't work with a rank-1 input matrix
  assert(inputMatrix.nRank > 1);

  // Assign some obvious values
  bcm.nRank = inputMatrix.nRank - 1;
  bcm.blockSize = inputMatrix.blockSize;

  // Find the sub-matrix "A" and invert it
  subMatrix_A = BlockCovarianceMatrix(inputMatrix.noisyCoeff[0],
				      inputMatrix.noiseFreeCoeff[0][0],
				      bcm.blockSize);
  subMatrix_A.InvertRankOneMatrix();

  // Initialise the arrays
  bcm.noisyCoeff = vector<double>(bcm.nRank);
  bcm.noiseFreeCoeff=vector<vector<double> >(bcm.nRank, vector<double>(bcm.nRank));

  // Assign the sub-matrix values
  const double temp_newNoiseFree = (subMatrix_A.noisyCoeff[0] + bcm.blockSize)*subMatrix_A.noiseFreeCoeff[0][0]*bcm.blockSize;
  for (int i=0; i<bcm.nRank; i++)
  {
    temp2_newNoiseFree=temp_newNoiseFree*inputMatrix.noiseFreeCoeff[i+1][0]; // B
    for (int j=0; j<bcm.nRank; j++)
    {
      // find adjusting noiseFree term
      newNoiseFree = temp2_newNoiseFree * inputMatrix.noiseFreeCoeff[0][j+1]; // C

      // make adjustments
      bcm.noiseFreeCoeff[i][j] = inputMatrix.noiseFreeCoeff[i+1][j+1] - newNoiseFree;
    }
    // find and adjust the noisy term for this row
    noisyScaling  = inputMatrix.noiseFreeCoeff[i+1][i+1] / bcm.noiseFreeCoeff[i][i];
    bcm.noisyCoeff[i] = inputMatrix.noisyCoeff[i+1] * noisyScaling;
  }

  return bcm;
}

/* ----------------------------------------------------------------------
   Factory method: builds "E", the rank N-1 sub-matrix, when [A] has a
   single missing observation.
---------------------------------------------------------------------- */

BlockCovarianceMatrix BlockCovarianceMatrix::
Build_E_SubMatrixMissingSingleObservation(const BlockCovarianceMatrix& inputMatrix)
{
  // The object we are going to construct
  BlockCovarianceMatrix bcm;

  // Declarations
  double noisyScaling, newNoiseFree, temp2_newNoiseFree;
  BlockCovarianceMatrix subMatrix_A;

  // This won't work with a rank-1 matrix
  assert(inputMatrix.nRank > 1);

  // Assign some obvious values
  bcm.nRank     = inputMatrix.nRank - 1;
  bcm.blockSize = inputMatrix.blockSize;

  // Find the sub-matrix "A" and invert it
  subMatrix_A = BlockCovarianceMatrix(inputMatrix.noisyCoeff[0],
				      inputMatrix.noiseFreeCoeff[0][0],
				      bcm.blockSize - 1);
  subMatrix_A.InvertRankOneMatrix();

  // Initialise the arrays
  bcm.noisyCoeff = vector<double>(bcm.nRank);
  bcm.noiseFreeCoeff = vector<vector<double> >(bcm.nRank, vector<double>(bcm.nRank));

  // Assign the sub-matrix values (pre-computation)
  const double temp_newNoiseFree = (subMatrix_A.noisyCoeff[0]+(bcm.blockSize - 1))*subMatrix_A.noiseFreeCoeff[0][0]*(bcm.blockSize - 1);
  for (int i=0; i<bcm.nRank; i++)
  {
    temp2_newNoiseFree = temp_newNoiseFree * inputMatrix.noiseFreeCoeff[i+1][0]; // B
    for (int j=0; j<bcm.nRank; j++)
    {
      // find adjusting noiseFree term
      newNoiseFree = temp2_newNoiseFree * inputMatrix.noiseFreeCoeff[0][j+1]; // C

      // make the adjustments
      bcm.noiseFreeCoeff[i][j] = inputMatrix.noiseFreeCoeff[i+1][j+1] - newNoiseFree;
    }
    // find and adjust the noisy term for this row
    noisyScaling  = inputMatrix.noiseFreeCoeff[i+1][i+1] / bcm.noiseFreeCoeff[i][i];
    bcm.noisyCoeff[i] = inputMatrix.noisyCoeff[i+1] * noisyScaling;
  }

  return bcm;
}

/* ----------------------------------------------------------------------
   Perform a recursive matrix inversion on this matrix.
---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertMatrix()
{
  if(nRank==1)
  {
    InvertRankOneMatrix();
  }
  else
  {
    InvertBlockMatrix();
  }
}

/* ----------------------------------------------------------------------
   Consider the matrix
      K = [A B]
          [C D]
   and define E = D-C*(A^-1)*B, then we use the identity
      K^-1 = [(A-B*(D^-1)*C)^-1  -(A^-1)*B*(E^-1)]
             [-(D^-1)*C*(A-B*(D^-1)*C)^-1  (E^-1)]
   to perform a block matrix pseudoinversion on this object.
---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertBlockMatrix()
{
  // Declarations
  int i, j;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;
  vector<double>        noiseFree_B, offDiagonal;
  double                diagonal_noisy, diagonal_noiseFree, factor_AB;

  // Use the block-matrix structure to our advantage;
  // find the sub-matrices A and E
  subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
  subMatrix_E = Build_E_SubMatrix(*this);

  // Invert A and E
  subMatrix_A.InvertRankOneMatrix();
  subMatrix_E.InvertMatrix();

  // Find the sub-matrices B and C
  // (symmetry => only need one of these)
  noiseFree_B = noiseFreeCoeff[0];
  // remove the first element as this is part of A, not B
  noiseFree_B.erase(noiseFree_B.begin(), noiseFree_B.begin()+1);

  // Find the factor coming from (A^-1)*B
  factor_AB  = subMatrix_A.noiseFreeCoeff[0][0];
  factor_AB *= blockSize + subMatrix_A.noisyCoeff[0];

  // Construct the overall inverse matrix;
  // copy E^-1 into the D slots
  for (i=1; i<nRank; i++)
  {
    noisyCoeff[i] = subMatrix_E.noisyCoeff[i-1];
    for (j=1; j<nRank; j++)
    {
      noiseFreeCoeff[i][j] = subMatrix_E.noiseFreeCoeff[i-1][j-1];
    }
  }
  // construct the remaining off-diagonal elements, using B * E^-1
  offDiagonal = subMatrix_E.BlockMultiply(noiseFree_B);
  for (i=1; i<nRank; i++)
  {
    // using the fact that our matrix is symmetric
    noiseFreeCoeff[0][i] = noiseFreeCoeff[i][0] = -factor_AB * offDiagonal[i-1];
  }
  // construct the final diagonal element
  diagonal_noiseFree = 0;
  for (i=0; i<(nRank-1); i++) // need to find B * E^-1 * C here...
  {
    // using the fact that C_transpose = B
    diagonal_noiseFree += blockSize * offDiagonal[i] * noiseFree_B[i];
  }
  diagonal_noiseFree *= factor_AB * factor_AB;
  // also need to add on A^-1
  diagonal_noiseFree += subMatrix_A.noiseFreeCoeff[0][0];
  // hence compute the new noisy term
  diagonal_noisy = subMatrix_A.noisyCoeff[0] * subMatrix_A.noiseFreeCoeff[0][0];
  diagonal_noisy /= diagonal_noiseFree;
  // and store the values in this object
  noiseFreeCoeff[0][0] = diagonal_noiseFree;
  noisyCoeff[0] = diagonal_noisy;
}

/* ----------------------------------------------------------------------
   Invert analytically a rank-1 matrix.
---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertRankOneMatrix()
{
  // Declarations
  double newNoisy, newNoiseFree;

  // Compute the new noisy value (for the diagonal elements)
  newNoisy  = -noisyCoeff[0] - blockSize;

  // Compute the new noise-free value
  newNoiseFree  = -1.0 /
    (noiseFreeCoeff[0][0] * (noisyCoeff[0] * (noisyCoeff[0] + blockSize)));

  // Update the matrix values, so this matrix object is now inverted
  noisyCoeff[0] = newNoisy;
  noiseFreeCoeff[0][0] = newNoiseFree;
}

/* ----------------------------------------------------------------------
   Compute the matrix log-determinant for the matrix missing a single
   observation.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::
ComputeMatrixMissingSingleObservationLogDeterminant() const
{
  // Declarations
  double logDeterminant;
  BlockCovarianceMatrix subMatrix_A, subMatrix_Efirst;

  // extract small A and first sub E matrix
  subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize-1);
  subMatrix_Efirst = Build_E_SubMatrixMissingSingleObservation(*this);

  // compute the contributions from the sub-matrices
  logDeterminant  = subMatrix_A.ComputeRankOneMatrixDeterminant();
  logDeterminant += subMatrix_Efirst.ComputeMatrixDeterminant();

  return logDeterminant;
}

/* ----------------------------------------------------------------------
   Compute the matrix inverse for the matrix missing a single observation
   (the matrix loses its first row and column).
---------------------------------------------------------------------- */

void BlockCovarianceMatrix::InvertMatrixMissingSingleObservation()
{
  // Declarations
  int i, j;
  BlockCovarianceMatrix subMatrix_A1, subMatrix_E1;
  vector<double>        noiseFree_B, offDiagonal;
  double                diagonal_noisy, diagonal_noiseFree, factor_AB;

  //M = [A,B][C,D], A is a smaller rank-1 matrix, since missing observation

  // Recursion to find the inverse matrix;
  // find the sub-matrices A and E
  subMatrix_A1 = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize-1);
  subMatrix_E1 = Build_E_SubMatrixMissingSingleObservation(*this);

  // Invert A and E
  subMatrix_A1.InvertRankOneMatrix();
  subMatrix_E1.InvertMatrix();

  // Find the sub-matrices B and C
  // (symmetry => only need one of these)
  //the first row of the M.noiseFreeCoeff
  noiseFree_B = noiseFreeCoeff[0];
  // remove the first element as this is part of A, not B
  noiseFree_B.erase(noiseFree_B.begin(), noiseFree_B.begin()+1);

  // Find the factor coming from (A^-1)*B
  factor_AB  = subMatrix_A1.noiseFreeCoeff[0][0];
  factor_AB *= (blockSize - 1) + subMatrix_A1.noisyCoeff[0];

  // Construct the overall inverse matrix;
  // copy E^-1 into the D slots
  for (i=1; i<nRank; i++)
  {
    noisyCoeff[i] = subMatrix_E1.noisyCoeff[i-1]; // copying the diagonal elements
    for (j=1; j<nRank; j++)
    {
      noiseFreeCoeff[i][j] = subMatrix_E1.noiseFreeCoeff[i-1][j-1];
    }
  }
  // construct the remaining off-diagonal elements, using B * E^-1
  offDiagonal = subMatrix_E1.BlockMultiply(noiseFree_B);
  for (i=1; i<nRank; i++)
  {
    // using the fact our matrix is symmetric
    noiseFreeCoeff[0][i] = noiseFreeCoeff[i][0] = -factor_AB * offDiagonal[i-1];
  }
  // construct the final diagonal element
  diagonal_noiseFree = 0;
  for (i=0; i<(nRank-1); i++) // need to find B * E_inv * C here...
  {
    // using the fact that C_transpose = B, dont subtract 1 from blockSize
    // size going across rows
    diagonal_noiseFree += blockSize * offDiagonal[i] * noiseFree_B[i];
  }
  diagonal_noiseFree *= factor_AB * factor_AB; // as is
  // also need to add on A_inverse
  diagonal_noiseFree += subMatrix_A1.noiseFreeCoeff[0][0]; // as is
  // hence compute the new noisy term
  diagonal_noisy = subMatrix_A1.noisyCoeff[0] * subMatrix_A1.noiseFreeCoeff[0][0];
  diagonal_noisy /= diagonal_noiseFree;
  // and store the values in this object
  noiseFreeCoeff[0][0] = diagonal_noiseFree;
  noisyCoeff[0] = diagonal_noisy;
}

/* ----------------------------------------------------------------------
   Compute the matrix log-determinant using the identity
   
      det[A B][C D] = det(A)det(E)
      
   where E=D-C*(A^-1)*B.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::ComputeMatrixDeterminant() const
{
  // Declarations
  double logDeterminant;
  BlockCovarianceMatrix subMatrix_A, subMatrix_E;//use the maths notation here

  // Recursion to find the log-det
  if (nRank==1)
  {
    logDeterminant = ComputeRankOneMatrixDeterminant();
  }
  else
  {
    // extract rank n-1 and rank-1 sub-matrices
    // remember that we want "E"
    subMatrix_A = BlockCovarianceMatrix(noisyCoeff[0], noiseFreeCoeff[0][0], blockSize);
    subMatrix_E = Build_E_SubMatrix(*this);

    // compute the contributions from the sub-matrices
    logDeterminant  = subMatrix_A.ComputeRankOneMatrixDeterminant();
    logDeterminant += subMatrix_E.ComputeMatrixDeterminant();
  }

  return logDeterminant;
}

/* ----------------------------------------------------------------------
   Compute analytically the rank-1 matrix determinant.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::ComputeRankOneMatrixDeterminant() const
{
  // Declarations
  double logDeterminant;

  assert(nRank==1);

  if (noiseFreeCoeff[0][0] < 0.0 || noisyCoeff[0] < 0.0)
  {
    cout << "problem with log Determinant: nan" << endl;
  }

  //Compute the log-det
  logDeterminant  = log(noiseFreeCoeff[0][0]) * blockSize;
  logDeterminant += log(noisyCoeff[0]) * (blockSize - 1);
  logDeterminant += log(noisyCoeff[0] + blockSize);

  return logDeterminant;
}

/* ----------------------------------------------------------------------
   Method to block-multiply this block-matrix and a vector.
---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::
BlockMultiply(const vector<double>& inputVector) const
{
  // Declarations
  int i;
  double currentElement;
  const int input_size=inputVector.size();
  vector<double> outputVector=vector<double>(input_size);

  // Compute the output vector elements
  for (i=0; i<input_size; i++)
  {
    // Inner product: <noiseFreeCoeff row, inputVector>
    currentElement = std::inner_product(inputVector.begin(),
					inputVector.end(),
					noiseFreeCoeff[i].begin(),
					0.0);

    // Since it's a block matrix, we use each matrix term blockSize-many times
    currentElement *= blockSize;

    // Add the (diagonal) noise term
    currentElement += noisyCoeff[i] * noiseFreeCoeff[i][i] * inputVector[i];

    // Store the result
    outputVector[i] = currentElement;
  }

  return outputVector;
}

/* ----------------------------------------------------------------------
   Method to find the (i,j)-th element of the block-matrix.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::GetElement(const int i, const int j) const
{
  // Declarations
  int block_i, block_j;
  double matrixElement;

  // Find which block the element belongs to
  block_i = i / blockSize;
  block_j = j / blockSize;

  // Find the matrix element
  matrixElement = noiseFreeCoeff[block_i][block_j];

  // If it's a diagonal element we make an adjustment
  if (i==j) matrixElement *= 1.0 + noisyCoeff[block_i];

  return matrixElement;
}

/* ----------------------------------------------------------------------
   Method to find the i-th row of the block matrix when a single
   observation is missing.
---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::
GetRowForMatrixMissingSingleObservation(const int index) const
{
  // Declarations
  int newindex, block_i;
  vector<double> rowVector = vector<double>(nRank*blockSize);

  // Make alteration to the index
  newindex = index + 1;

  // Find which block the row belongs to
  block_i = newindex / blockSize;

  // Find the matrix elements (noise-free component
  vector<double>::const_iterator inIt = noiseFreeCoeff[block_i].begin();
  vector<double>::iterator rowVectorIt = rowVector.begin();
  while (inIt != noiseFreeCoeff[block_i].end())
  {
    vector<double>::iterator outItEnd = rowVectorIt + blockSize;
    fill(rowVectorIt, outItEnd, *inIt++);
    rowVectorIt = outItEnd;
  }
  // the above is a much more efficient version of the following:
  /*
  for (i=0; i<nRank; i++)
  {
    for (j=0; j<blockSize; j++)
    {
      rowVector.push_back(noiseFreeCoeff[block_i][i]);
    }
  }
  */

  // Adjust the diagonal element
  rowVector[newindex] *= 1 + noisyCoeff[block_i];
  // the first element of rowVector corresponds to the missing observation
  rowVector.erase(rowVector.begin(), rowVector.begin()+1);

  return rowVector;
}

/* ----------------------------------------------------------------------
   Method to build the full i-th row of the block-matrix.
---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::GetRow(const int index) const
{
  // Declarations
  vector<double> rowVector=vector<double>(nRank*blockSize);

  // Find the block this row belongs to
  int block_i = index / blockSize;

  // Find the matrix elements (noise-free component)
  vector<double>::const_iterator inIt = noiseFreeCoeff[block_i].begin();
  vector<double>::iterator rowVectorIt = rowVector.begin();
  while (inIt != noiseFreeCoeff[block_i].end())
    {
      vector<double>::iterator outItEnd = rowVectorIt + blockSize;
      fill(rowVectorIt, outItEnd, *inIt++);
      rowVectorIt = outItEnd;
    }
  // the above is a much more efficient version of the following:
  /*
    for (i=0; i<nRank; i++)
    {
    for (j=0; j<blockSize; j++)
    {
    rowVector.push_back(noiseFreeCoeff[block_i][i]);
    }
    }
  */

  // Adjust the diagonal element
  rowVector[index] *= 1.0 + noisyCoeff[block_i];

  return rowVector;
}

/* ----------------------------------------------------------------------
   Method to compute the bilinear form, y_transpose * K * y, where y=data,
   for K missing a single observation.
   In practice, K should have been inverted before calling this.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::
ComputeLogLikelihoodProductMissingSingleObservation(const vector<double>& data) const
{
  // Declarations
  int k;
  double logLike, currentValue;
  vector<double> currentRow;
  vector<double>::const_iterator dataIt=data.begin();

  assert(data_size==(nRank*blockSize-1)); //Should be missing an observation

  // Form the product K*y
  logLike=0;
  // going across the rows in the K matrix
  for(k=0; dataIt != data.end(); ++dataIt, ++k)
  {
    // get the current row from the covariance matrix
    currentRow = GetRowForMatrixMissingSingleObservation(k);
    // loop to sum the terms in a vector*matrix product
    currentValue = std::inner_product(data.begin(), data.end(), currentRow.begin(), 0.0);
    // now form the contribution to the overall product, y_t * K * y
    logLike += (*dataIt) * currentValue;
  }

  return logLike;
}

/* ----------------------------------------------------------------------
   Method to compute the bilinear form, y_transpose * K * y, where y=data.
   In practice, K should have been inverted before calling this.
---------------------------------------------------------------------- */

double BlockCovarianceMatrix::
ComputeLogLikelihoodProduct(const vector<double>& data) const
{
  // Compute K * y
  vector<double> temp=VectorMultiply(data);
  
  // Compute and return y^t * (K * y)
  return std::inner_product(temp.begin(), temp.end(), data.begin(), 0.0);
}

/* ----------------------------------------------------------------------
   Compute the product K * y, where y=inputVector.
---------------------------------------------------------------------- */

vector<double> BlockCovarianceMatrix::
VectorMultiply(const vector<double>& inputVector) const
{
  // Declarations
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

/* ---------------------------------------------------------------------- */

