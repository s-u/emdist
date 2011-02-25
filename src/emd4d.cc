#include <cmath>

typedef struct location {
  double x, y, z, a;
} location_t, *location_ptr;

static double dist(location_ptr a, location_ptr b)
{
  return sqrt((a->x - b->x)*(a->x - b->x) +
	      (a->y - b->y)*(a->y - b->y) +
	      (a->z - b->z)*(a->z - b->z) +
	      (a->a - b->a)*(a->a - b->a)
	      );
}

#include "em_api.h"

typedef EarthMover<location_ptr> EM;

#define R_NO_REMAP 1
#define USE_RINTERNALS 1
#include <Rinternals.h>

extern "C" SEXP emd_4d(SEXP sBase, SEXP sCur) {
  SEXP sBaseDim = Rf_getAttrib(sBase, R_DimSymbol);
  SEXP sCurDim = Rf_getAttrib(sCur, R_DimSymbol);
  if (sBaseDim == R_NilValue || LENGTH(sBaseDim) != 2) Rf_error("base must be a matrix");
  if (sCurDim  == R_NilValue || LENGTH(sCurDim)  != 2) Rf_error("cur must be a matrix");
  int *baseDim = INTEGER(sBaseDim);
  int *curDim = INTEGER(sCurDim);
  int baseRows = baseDim[0], baseCol = baseDim[1];
  int curRows = curDim[0], curCol = curDim[1];
  sBase = Rf_coerceVector(sBase, REALSXP);
  sCur = Rf_coerceVector(sCur, REALSXP);
  double *baseVal = REAL(sBase);
  double *curVal = REAL(sCur);

  if (baseCol != curCol) Rf_error("base and current sets must have the same dimensionality");
  if (baseCol < 2) Rf_error("at least two columns (weight and location) are required");
  if (baseCol > 5) Rf_warning("more than four dimensions are used, those will be ignored");

  std::vector<location_t> baseloc(baseRows);
  std::vector<location_t> curloc(curRows);
  std::vector< std::pair<location_ptr,int> > baseweights(baseRows);
  std::vector< std::pair<location_ptr,int> > currentweights(curRows);

  for (int i = 0; i < baseRows; i++) {
    baseloc[i].x = baseVal[i + baseRows];
    baseloc[i].y = (baseCol > 2) ? baseVal[i + 2 * baseRows] : 0.0;
    baseloc[i].z = (baseCol > 3) ? baseVal[i + 3 * baseRows] : 0.0;
    baseloc[i].a = (baseCol > 4) ? baseVal[i + 4 * baseRows] : 0.0;
    baseweights[i].first = &baseloc[i];
    baseweights[i].second = (int) (baseVal[i] * 1000000.0);
#ifdef EMD_DEBUG
    Rprintf("A: (%g,%g,%g,%g)->%d\n", baseloc[i].x, baseloc[i].y, baseloc[i].z, baseloc[i].a, baseweights[i].second);
#endif
  }
  for (int i = 0; i < curRows; i++) {
    curloc[i].x = curVal[i + curRows];
    curloc[i].y = (curCol > 2) ? curVal[i + 2 * curRows] : 0.0;
    curloc[i].z = (curCol > 3) ? curVal[i + 3 * curRows] : 0.0;
    curloc[i].a = (curCol > 4) ? curVal[i + 4 * curRows] : 0.0;
    currentweights[i].first = &curloc[i];
    currentweights[i].second = (int) (curVal[i] * 1000000.0);
#ifdef EMD_DEBUG
    Rprintf("B: (%g,%g,%g,%g)->%d\n", curloc[i].x, curloc[i].y, curloc[i].z, curloc[i].a, currentweights[i].second);
#endif
  }
  
  double d = EM(baseweights, currentweights).distance();
#ifdef EMD_DEBUG
  Rprintf("<A, B> = %g\n", d);
#endif
  
  return Rf_ScalarReal(d);
}
