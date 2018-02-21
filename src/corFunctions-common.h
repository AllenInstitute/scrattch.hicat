/* 
 * Common functions for fast calculations of correlations
 *
 */


#ifndef __corFunctions_common_h__
#define __corFunctions_common_h__

#define LDOUBLE 	long double

#include "pivot.h"

enum { noWarning, warnZeroMAD };

double median(double * x, int n, int copy, int * err);

double quantile(double * x, int n, double q, int copy, int * err);

void testMedian(double *x, int * n, double * res);

void prepareColBicor(double * col, int nr, double maxPOutliers, int fallback, int cosine, double * res, 
                     int * nNAentries, int * NAmed, volatile int *zeroMAD, double * aux, double * aux2);

void prepareColCor(double * x, int nr, int cosine, double * res, int * nNAentries, int * NAmean);

#endif
