/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef BLOCKCOVARIANCEMATRIX_H
#define BLOCKCOVARIANCEMATRIX_H

#include "header.h"

/* ----------------------------------------------------------------------
   This class handles the covariance function matrices for the time-series
   version of BHC. In this case, the data all lie at a set of time-points that
   are common to all time-series. This means the covaraince function matrix
   has a characteristic block-matrix form that we can use to massively speed up
   the calculation of matrix inverses and determinants.
   
   We assume the following for this matrix:
      - it is formed from (nRank*nRank) blocks
      - each block is (blockSize*blockSize)
      - the diagonal blocks take the form a(bI + J)
        where I is the identity matrix and J the unit matrix
      - the off-diagonal blocks take the form aJ
---------------------------------------------------------------------- */

class BlockCovarianceMatrix
{
 public:
  BlockCovarianceMatrix();
  BlockCovarianceMatrix(const double noisyValue,
			const double noiseFreeValue,
			const int inputBlockSize);
  static BlockCovarianceMatrix
    Build_E_SubMatrix(const BlockCovarianceMatrix& inputMatrix);
  static BlockCovarianceMatrix
    Build_E_SubMatrixMissingSingleObservation(const BlockCovarianceMatrix& inputMatrix);
  void InvertMatrix();
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
  double GetElement(const int i, const int j) const;
  vector<double> GetRow(const int i) const;
  vector<double> GetRowForMatrixMissingSingleObservation(const int i) const;

  
  int nRank; // the number of blocks in a single dimension
  int blockSize; // the number of elements in a single dimension of a block
  vector<vector<double> > noiseFreeCoeff; // 'a' values -  noiseless component
  vector<double> noisyCoeff; // 'b' values -  noisy, diagonal component
};

#endif // BLOCKCOVARIANCEMATRIX_H
