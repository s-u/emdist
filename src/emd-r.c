/* this version of the EMD uses Rubner's code */
#include "emd-rubner.h"

#include <stdlib.h>

#define R_NO_REMAP 1
#define USE_RINTERNALS 1
#include <Rinternals.h>

SEXP emd_r(SEXP sBase, SEXP sCur) {
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
  if (baseCol > FDIM + 1) Rf_warning("more than %d dimensions are used, those will be ignored", FDIM);

  signature_t baseSig, curSig;
  
  baseSig.n = baseRows;
  baseSig.Features = calloc(baseRows, sizeof(signature_t));
  baseSig.Weights  = calloc(baseRows, sizeof(float));
  curSig.n = curRows;
  curSig.Features = calloc(curRows, sizeof(signature_t));
  curSig.Weights  = calloc(curRows, sizeof(float));

  int i, j;
  for (i = 0; i < baseRows; i++) {
    for (j = 0; j < baseCol; j++)
      if (j < FDIM)
	baseSig.Features[i].loc[j] = baseVal[i + (j + 1) * baseRows];
    baseSig.Weights[i] = baseVal[i];
  }
  for (i = 0; i < curRows; i++) {
    for (j = 0; j < curCol; j++)
      if (j < FDIM)
	curSig.Features[i].loc[j] = curVal[i + (j + 1) * curRows];
    curSig.Weights[i] = curVal[i];
  }
  
  double d = emd_rubner(&baseSig, &curSig, NULL,NULL);
  
  return Rf_ScalarReal(d);
}
