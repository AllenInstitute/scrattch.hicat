/*
 * Compilation:
 *  gcc --std=c99 -fPIC -O3 -o pivot.so -shared pivot.c
 *
 */


#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

#include "pivot.h"


void RprintV(double * v, size_t l)
{
  for (size_t i=0; i<l; i++) Rprintf("%5.3f ", v[i]);
  Rprintf("\n");
}

double vMax(double * v, size_t len)
{
  double mx = v[0];
  for (size_t i=1; i<len; i++)
    if (v[i] > mx) mx = v[i];
  return mx;
}

double vMin(double * v, size_t len)
{
  double mn = v[0];
  for (size_t i=1; i<len; i++)
    if (v[i] < mn) mn = v[i];
  return mn;
}


double pivot(double * v, size_t len, double target)
{
  // Rprintf("Entering pivot with len=%d and target=%f\n   ", len, target);
  // RprintV(v, len);

  if (len > 2)
  {
    // pick the pivot, say as the median of the first, middle and last
    size_t i1 = 0, i2 = len-1, i3 = (len-1)/2, ip;
    if (v[i1] <= v[i2])
    {
      if (v[i2] <= v[i3])
        ip = i2;
      else if (v[i3] >= v[i1])
         ip = i3;
      else 
         ip = i1;
    } else {
      if (v[i1] <= v[i3])
        ip = i1;
      else if (v[i2] <= v[i3])
        ip = i3;
      else
        ip = i2;
    }

    // put ip at the end
    double vp = v[ip];
    v[ip] = v[len-1];
    v[len-1] = vp;

    // Rprintf("   pivot value: %5.3f, index: %d\n", vp, ip);

    // pivot everything else
    size_t bound = 0;
    for (size_t i=0; i<len; i++) if (v[i] < vp)
    {
      double x = v[bound];
      v[bound] = v[i];
      v[i] = x;
      bound++;
    }

    v[len-1] = v[bound];
    v[bound] = vp;

    // Rprintf("   After pivoting: bound:%d and vector: ", bound); // RprintV(v, len);

    // Did we find the target?
    
    double crit = target - bound;
    // Rprintf("   crit: %5.3f\n", crit);
    if (fabs(crit) > 1.0)
    {
      if (crit < 0)
        return pivot(v, bound, target);
      else
        return pivot(v+bound+1, len-bound-1, target-bound-1);
    }
    // Rprintf("vMax(v, bound): %5.3f, vMin(v+bound+1, len-bound-1): %5.3f, vp: %5.3f\n", vMax(v, bound),
                // vMin(v+bound+1, len-bound-1), vp);
    if (crit < 0)
    {
       double v1 = vMax(v, bound);
       return (v1 *(-crit) + vp * (1+crit));
    } // else
    double v2 = vMin(v+bound+1, len-bound-1);
    return (vp * (1-crit) + v2 * crit);
  } 
  else if (len==2)
  {
      // Rprintf("  Short v, returning a direct value.\n"); 
      double v1 = vMin(v, 2);
      double v2 = vMax(v, 2);
      if (target < 0) return v1;
      else if (target > 1) return v2;
      else return (target * v2 + (1-target) * v1);
  }
  else 
  {
     // Rprintf("  length 1 v, returning a direct value.\n"); 
     return v[0];
  }
}


/*
 *
 * This isn't needed for now.
 *
 *
void testPivot(double * v, size_t * len, double * target, double * result)
{
   * result = pivot(v, *len, *target);
}
*/

/*****************************************************************************************************
 *
 * Implement order via qsort.
 *
 *****************************************************************************************************/

int compareOrderStructure(const orderStructure * os1, const orderStructure * os2)
{
  if (ISNAN(os1->val)) return 1;
  if (ISNAN(os2->val)) return -1;
  if (os1->val < os2->val) return -1;
  if (os1->val > os2->val) return 1;
  return 0;
}

void qorder_internal(double * x, size_t n, orderStructure * os)
{
  for (R_xlen_t i = 0; i<n; i++)
  {
    (os+i)->val = *(x+i);
    (os+i)->index = i;
  }

  // Rprintf("qorder: calling qsort..");
  qsort(os, (size_t) n, sizeof(orderStructure), 
           ((int (*) (const void *, const void *)) compareOrderStructure));
}
  

SEXP qorder(SEXP data)
{
  R_xlen_t n = Rf_xlength(data);

  // Rprintf("qorder: length of input data is %ld.\n", n);

  double * x = REAL(data);

  orderStructure * os = Calloc((size_t) n, orderStructure);

  qorder_internal(x, (size_t) n, os);

  SEXP ans;
  if (n<(size_t) 0x80000000)
  {
    // Rprintf("..returning integer order.\n");
    PROTECT (ans = allocVector(INTSXP, n));
    int * ansp = INTEGER(ans);
    for (R_xlen_t i = 0; i<n; i++) ansp[i] = (int) ( (os+i)->index+1);
  } else {
    // Rprintf("..returning floating point (double) order.\n");
    PROTECT (ans = allocVector(REALSXP, n));
    double * ansp = REAL(ans);
    for (R_xlen_t i = 0; i<n; i++) ansp[i] = (double) ((os+i)->index+1);
  }
  Free(os);
  UNPROTECT(1);
  return ans;
}

