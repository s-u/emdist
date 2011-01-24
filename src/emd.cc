#include <cmath>

#include "em_api.h"

typedef struct location {
  double x;
  double y;
} location_t, *location_ptr;

typedef EarthMover<location_ptr> EM;

static double dist(location_ptr a, location_ptr b)
{
  return sqrt((a->x-b->x)*(a->x-b->x) + (a->y - b->y)*(a->y - b->y));
}

#define R_NO_REMAP 1
#define USE_RINTERNALS 1
#include <Rinternals.h>

extern "C" SEXP emd_2d(SEXP sBase, SEXP sCur, SEXP sXDist, SEXP sYDist) {
  SEXP sBaseDim = Rf_getAttrib(sBase, R_DimSymbol);
  if (sBaseDim == R_NilValue) Rf_error("base must be a matrix");
  int *baseDim = INTEGER(sBaseDim);
  int m = baseDim[0], n = baseDim[1], l = LENGTH(sBase);
  if (l != LENGTH(sCur)) Rf_error("base and current length mismatch");
  double *baseVal = REAL(sBase);
  double *curVal = REAL(sCur);
  double *xDist = REAL(sXDist);
  double *yDist = REAL(sYDist);
  int xl = LENGTH(sXDist);
  int yl = LENGTH(sYDist);

  if (xl != 1 && xl != n)
    Rf_error("column distance must be a vector of length one or the number of columns");
  if (yl != 1 && yl != m)
    Rf_error("row distance must be a vector of length one or the number of rows");

  std::vector<location_t> location(l);
  std::vector< std::pair<location_ptr,int> > baseweights(l);
  std::vector< std::pair<location_ptr,int> > currentweights(l);

  l = 0;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      location[l].x = (xl == 1) ? (xDist[0] * (double)j) : xDist[j];
      location[l].y = (yl == 1) ? (yDist[0] * (double)i) : yDist[i];
      baseweights[l].first = &location[l];
      baseweights[l].second = (int) (baseVal[l] * 1000000);
      currentweights[l].first = &location[l];
      currentweights[l].second = (int) (curVal[l] * 1000000);
      l++;
    }
  
  double d = EM(baseweights, currentweights).distance();
  
  return Rf_ScalarReal(d);
}
