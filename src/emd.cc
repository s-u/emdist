#include <cmath>

typedef struct location {
  double x;
  double y;
} location_t, *location_ptr;

static double dist(location_ptr a, location_ptr b)
{
  return sqrt((a->x-b->x)*(a->x-b->x) + (a->y - b->y)*(a->y - b->y));
}

#include "em_api.h"

typedef EarthMover<location_ptr> EM;

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
  int nzb = 0, nzc = 0;

  if (xl != 1 && xl != n)
    Rf_error("column distance must be a vector of length one or the number of columns");
  if (yl != 1 && yl != m)
    Rf_error("row distance must be a vector of length one or the number of rows");

  for (int i = 0; i < l; i++) {
    if (baseVal[i] != 0.0) nzb++;
    if (curVal[i] != 0.0) nzc++;
  }

  std::vector<location_t> baseloc(nzb);
  std::vector<location_t> curloc(nzc);
  std::vector< std::pair<location_ptr,int> > baseweights(nzb);
  std::vector< std::pair<location_ptr,int> > currentweights(nzc);

  l = 0;
  int bi = 0, ci = 0;
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      if (baseVal[l] != 0.0) {
	baseloc[bi].x = (xl == 1) ? (xDist[0] * (double)j) : xDist[j];
	baseloc[bi].y = (yl == 1) ? (yDist[0] * (double)i) : yDist[i];
	baseweights[bi].first = &baseloc[bi];
	baseweights[bi].second = (int) (baseVal[l] * 1000000);
	bi++;
      }
      if (curVal[l] != 0.0) {
	curloc[ci].x = (xl == 1) ? (xDist[0] * (double)j) : xDist[j];
	curloc[ci].y = (yl == 1) ? (yDist[0] * (double)i) : yDist[i];
	currentweights[ci].first = &curloc[ci];
	currentweights[ci].second = (int) (curVal[l] * 1000000);
	ci++;
      }
      l++;
    }
  
  double d = EM(baseweights, currentweights).distance();
  
  return Rf_ScalarReal(d);
}
