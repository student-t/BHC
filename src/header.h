#ifndef HEADER_H
#define HEADER_H

#define NDEBUG // no debug - disable things like asserts

// Handy macros
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

extern bool fast_switch; // Currently not used

// Includes from the STL
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
// From this package
#include "gammaln.h"
// R includes
#include "Rconfig.h"

#ifdef SUPPORT_OPENMP
#ifdef _OPENMP
#include <omp.h>
#else
#undef SUPPORT_OPENMP
#endif
#endif

using namespace std;

#endif // HEADER_H

