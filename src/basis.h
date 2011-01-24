/* RCSINFO $Id: basis.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __BASIS_H
#define __BASIS_H

#include "dstruct.h"

typedef struct var_data{
   int nartif;
   int nslacks;
   int nfree;
   int nbndone;
   int nbounded;
   int nfixed;
   double cmax;
} var_data; 

int
    ILLbasis_build_basisinfo (lpinfo *lp),
    ILLbasis_get_initial (lpinfo *lp, int algorithm),
    ILLbasis_get_cinitial (lpinfo *lp, int algorithm),
    ILLbasis_load (lpinfo *lp, ILLlp_basis *B),
    ILLbasis_factor (lpinfo *lp, int *singular),
    ILLbasis_refactor (lpinfo *lp),
    ILLbasis_update (lpinfo *lp, svector *y, int lindex, int *refactor, int *singular);

void
    ILLbasis_column_solve (lpinfo *lp, svector *rhs, svector *soln),
    ILLbasis_column_solve_update (lpinfo *lp, svector *rhs, svector *upd, svector *soln),
    ILLbasis_row_solve (lpinfo *lp, svector *rhs, svector *soln),
    ILLbasis_free_basisinfo (lpinfo *lp),
    ILLbasis_free_fbasisinfo (lpinfo *lp),
    ILLbasis_init_basisinfo (lpinfo *lp);

#endif /* __BASIS_H */
