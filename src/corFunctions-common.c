/* 
 * Common functions for fast calculations of correlations
 *
 */


/*
 * General notes about handling missing data, zero MAD etc:
 * The idea is that bicor should switch to cor whenever it is feasible, it helps, and it is requested:
 * (1) if median is NA, the mean would be NA as well, so there's no point in switching to Pearson
 * (2) In the results, columns and rows corresponding to input with NA means/medians are NA'd out.
 * (3) The convention is that if zeroMAD is set to non-zero, it is the index of the column in which MAD is
 *     zero (plus one for C indexing)
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>
#define LDOUBLE 	long double

#include "pivot.h"
#include "corFunctions-common.h"

#define RefUX	0.5

/*===================================================================================
 *
 * median
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the median of the remaining
// (finite) entries.

double median(double * x, int n, int copy, int * err)
{
  double * xx, res;
  if (copy)
  {
    if ( (xx=malloc(n * sizeof(double)))==NULL ) 
    {
      Rprintf("Memory allocation error in median(). Could not allocate %d kB.\n", 
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;

    
  *err = 0;
  // Put all NA's at the end.
  int bound = n-1;
  for (int i=n-1; i>=0; i--) if (ISNAN(xx[i]))
  {
     xx[i] = xx[bound];
     xx[bound] = NA_REAL;
     bound--;
  }

  // Rprintf("Median: n: %d, bound: %d\n", n, bound);
  // Any non-NA's left?

  if (bound<0)
    res = NA_REAL;
  else 
  // yes, return the appropriate pivot. 
    res = pivot(xx, bound+1, ( 1.0 * bound)/2);

  if (copy) free(xx);

  return res;

}


/*===================================================================================
 *
 * quantile
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the appropriate
// quantile of the remaining (finite) entries.

// q is the quantile: 1/2 will give exactly the median above.

double quantile(double * x, int n, double q, int copy, int * err)
{
  double * xx;
  double res;

  if (copy)
  {
    if ( (xx=malloc(n * sizeof(double)))==NULL ) 
    {
      Rprintf("Memory allocation error in quantile(). Could not allocate %d kB.\n", 
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;

    
  *err = 0;
  // Put all NA's at the end.
  int bound = n-1;
  for (int i=n-1; i>=0; i--) if (ISNAN(xx[i]))
  {
     xx[i] = xx[bound];
     xx[bound] = NA_REAL;
     bound--;
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound<0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot. 
    res = pivot(xx, bound+1, ( 1.0 * bound)*q);

  if (copy) free(xx);

  return res;

}



/*==========================================================================================
 *
 * testMedian
 *
 * =========================================================================================*/


void testMedian(double *x, int * n, double * res)
{
  int err;
  *res = median(x, *n, 0, &err);
} 

/*==========================================================================================
 *
 * testQuantile
 *
 * =========================================================================================*/


void testQuantile(double *x, int *n, double *q, double *res)
{
  int err;
  *res = quantile(x, *n, *q, 0, &err);
} 


/*==========================================================================================
 *
 * prepareColBicor
 *
 * =========================================================================================*/

// prepareColBicor: calculate median and mad of x and put 
// (1-u^2)^2 * (x - median(x))/(9.0 * qnorm75 * mad(x))/ appropriate normalization 
// into res. 
// res must have enough space allocated to hold the result; 
// aux and aux2 each must also have enough space to hold a copy of x.

// maxQoutliers is the maximum allowed proportion of outliers on either side of the median.

// fallback: 1: none, 2: individual, 3: all, 4: force Pearson calculation. 4 is necessary for remedial
// calculations.

// In this case: Pearson pre-calculation entails normalizing columns by mean and variance. 

void prepareColBicor(double * col, int nr, double maxPOutliers, int fallback,
                     int cosine,
                     double * res, int * nNAentries, 
                     int * NAmed, volatile int * zeroMAD,
                     double * aux, double * aux2)
{
  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;
  // Note to self: asymptCorr * qnorm75 is very close to 1 and should equal 1 theoretically. Should
  // probably leave them out completely. 

  if (fallback==4)
  {
    prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
    return;
  }

  int err = 0;

  // Calculate the median of col

  memcpy((void *)res, (void *)col, nr * sizeof(double));
  double med = median(res, nr, 0, &err);

  // Create a conditional copy of the median
  double medX;
  if (cosine) medX = 0; else medX = med;

  *zeroMAD = 0;
  // calculate absolute deviations from the median

  if (ISNAN(med))
  {
    *NAmed = 1;
    for (int k=0; k<nr; k++) *(res + k) = 0;
  } else {
    *NAmed = 0;
    *nNAentries = 0;
    for (int k=0; k<nr; k++)
      if (ISNAN(col[k]))
      {
        (*nNAentries)++;
        res[k] = NA_REAL;
        aux[k] = NA_REAL;
      } else {
        res[k] = col[k] - medX;
        aux[k] = fabs(col[k] - med);
      }

    // calculate mad, i.e. median absolute deviation
    double mad = median(aux, nr, 0, &err);

    // If mad is zero, value of fallback decides what is it we will do.
    if (mad==0)
    {
       *zeroMAD = 1;
       switch (fallback)
       {
          case 1: 
          {
             // Return after zeoring out results and setting the NAmed flag
             for (int k=0; k<nr; k++) res[k] = 0;
             *NAmed = 1;
             return;
          }
          case 2: 
             // Switch to Pearson correlation and return
             // Rprintf("mad is zero in a column. Switching to Pearson for this column.\n");
             prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
          case 3: 
             // Do nothing: the setting of *zeroMAD above is enough.
             return;
       }
    } 

    // We now re-use aux to store a copy of the weights ux. To calculate them, first get (x-med)/(9*mad).

    // Rprintf("median: %6.4f, mad: %6.4f, cosine: %d\n", med, mad, cosine);

    double denom = 9.0 * mad;
    for (int k=0; k<nr; k++)
      if (!ISNAN(col[k]))
        aux[k] = (col[k] - med) / denom;
      else
        aux[k] = NA_REAL;

    // Get the low and high quantiles  of ux
    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double lowQ = quantile(aux2, nr, maxPOutliers, 0, &err);

    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double hiQ = quantile(aux2, nr, 1-maxPOutliers, 0, &err);

    // Rprintf("prepareColBicor: lowQ=%f, hiQ = %f\n", lowQ, hiQ);
    // If the low quantile is below -1, rescale the aux (that serve as ux below)
    // such that the low quantile will fall at -1; similarly for the high quantile

    if (lowQ > -RefUX) lowQ = -RefUX;
    if (hiQ < RefUX) hiQ = RefUX;
    lowQ = fabs(lowQ);

    for (int k=0; k<nr; k++) if (!ISNAN(aux[k]))
    {
      if (aux[k] < 0)
        aux[k] = aux[k] * RefUX / lowQ;
      else
        aux[k] = aux[k] * RefUX / hiQ;
    }

    // Calculate the (1-ux^2)^2 * (x-median(x)) 

    LDOUBLE sum = 0;
    for (int k=0; k<nr; k++)
      if (!ISNAN(res[k]))
      {
        double ux = aux[k];
        if (fabs(ux) > 1) ux = 1;  // sign of ux doesn't matter.
        ux = 1-ux*ux;
        res[k] *= ux*ux ;
        sum += res[k]*res[k];
      } else
        res[k] = 0;
    sum = sqrtl(sum);
    if (sum==0)
    {
       for (int k=0; k<nr; k++)
          res[k] = 0;
       *NAmed = 1;
    } else {
       for (int k=0; k<nr; k++)
          res[k] = res[k] / sum;
    }
  }
}

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for Pearson correlation fallback
// and when bicor is called with robustsX or robustY = 0
// if cosine is not zero, the cosine correlation will be calculated.

void prepareColCor(double * x, int nr, int cosine, double * res, int * nNAentries, int * NAmean)
{
  *nNAentries = 0;
  int count = 0;
  LDOUBLE mean = 0, sum = 0;
  for (int k = 0; k < nr; k++)
    if (!ISNAN(x[k]))
    {
      mean += x[k];
      sum += ((LDOUBLE) x[k])*( (LDOUBLE) x[k]);
      count ++;
    }
  if (count > 0)
  {
    *NAmean = 0;
    *nNAentries = nr-count;
    if (cosine) mean = 0; else mean = mean/count;
    sum = sqrtl(sum - count * mean*mean);
    if (sum > 0)
    {
      // Rprintf("sum: %Le\n", sum);
       for (int k=0; k<nr; k++)
         if (!ISNAN(x[k]))
            res[k] = (x[k] - mean)/sum;
         else
            res[k] = 0;
    } else {
       // Rprintf("prepareColCor: have zero variance.\n");
       *NAmean = 1;
       for (int k=0; k<nr; k++) res[k] = 0;
    }
  } else {
    *NAmean = 1;
    *nNAentries = nr;
    for (int k=0; k<nr; k++)
       res[k] = 0;
  }
}

