/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef HEADER_H
#define HEADER_H

#define NDEBUG

// Handy macros
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

extern bool fast_switch; // declared in header.cpp

// Standard includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <climits>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <numeric>
#include <ctime>

// From this package
#include "gammaln.h"

// TO ENABLE OPENMP
// ----------------
// If building the R package, then uncomment this line:
#include <Rconfig.h>
// else, uncomment this line:
//#define SUPPORT_OPENMP 1

#ifdef SUPPORT_OPENMP // R-support
#ifdef _OPENMP // Compiler-support
#include <omp.h>
#else
#undef SUPPORT_OPENMP
#endif
#endif

using namespace std;

#endif // HEADER_H

