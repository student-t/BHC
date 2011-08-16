#ifndef BLOCKCOVARIANCEMATRIX_H
#define BLOCKCOVARIANCEMATRIX_H

#include "header.h"

// This class handles the covariance function matrices for the time-series
// version of BHC. In this case, the data all lie at a set of time-points that
// are common to all time-series. This meams the covaraince function matrix
// has a characteristic block-matrix form that we can use to massively speed up
// the calculation of matrix inverse and determinants.
//
// We assume the following for this matrix:
//  - it is formed from (nRank*nRank) blocks
//  - each block is (blockSize*blockSize)
//  - the diagonal blocks take the form a(bI + J)
//    (I is the identity matrix; J is the unit matrix)
//  - the off-diagonal blocks take the form aJ
class BlockCovarianceMatrix
{
public:

  //CONSTRUCTORS
  BlockCovarianceMatrix();
  BlockCovarianceMatrix(const double noisyValue,
			const double noiseFreeValue,
			const int inputBlockSize);

  //FACTORY METHODS
  static BlockCovarianceMatrix
    Build_E_SubMatrix(const BlockCovarianceMatrix& inputMatrix);
  static BlockCovarianceMatrix
    Build_E_SubMatrixMissingSingleObservation(const BlockCovarianceMatrix& inputMatrix);

  //OTHER METHODS
  void InvertMatrix();
  void InvertGeneralMatrix();
  void InvertBlockMatrix();
  void InvertRankOneMatrix();
  void InvertMatrixMissingSingleObservation();
  double ComputeMatrixDeterminant() const;
  double ComputeMatrixMissingSingleObservationLogDeterminant() const;
  double ComputeRankOneMatrixDeterminant() const;
  vector<double> BlockMultiply(const vector<double>& inputVector) const;
  double ComputeLogLikelihoodProduct(const vector<double>& data) const;
  double 
    ComputeLogLikelihoodProductMissingSingleObservation(const vector<double>& data) const;
  vector<double> VectorMultiply(const vector<double>& inputVector) const;

  //GETS AND SETS
  double GetElement(const int i, const int j) const;
  vector<double> GetRow(const int i) const;
  vector<double> GetRowForMatrixMissingSingleObservation(const int i) const;

  //TAGS
  //the number of blocks in a single dimension
  int nRank;
  //the number of elements in a single dimension of a block
  int blockSize;
  //'a' values that define the noiseless component of the blocks
  vector<vector<double> > noiseFreeCoeff;
  //'b' values that define the noisy, diagonal component
  vector<double> noisyCoeff;
};

#endif // BLOCKCOVARIANCEMATRIX_H
