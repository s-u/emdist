/* $RCSfile: simplex.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include "lpdata.h"

/* possible values of nextstep */
#define SIMPLEX_CONTINUE   1
#define SIMPLEX_TERMINATE  2
#define SIMPLEX_RESUME     3

/* reason for resuming simplex */
#define SIMPLEX_RESUME_SING     1
#define SIMPLEX_RESUME_UNSHIFT  2
#define SIMPLEX_RESUME_NUMER    3

/* values for newphase */
#define SIMPLEX_PHASE_RECOMP  1
#define SIMPLEX_PHASE_NEW     2

#define SIMPLEX_PIVOTINROW 1
#define SIMPLEX_PIVOTINCOL 2

typedef struct param_info{
   int  origalgo;
   int  pphaseI;
   int  pphaseII;
   int  dphaseI;
   int  dphaseII;
   int  p_strategy;
   int  d_strategy;
} param_info;

typedef struct iter_info{
   int  newphase;
   int  nextphase;
   int  nextstep;
   int  sdisplay;
   int  itercnt;
   int  solstatus;
   int  curtime;
   int  rounds;
   int  chkobj;
   int  nosolve;
   int  noprog;
   int  inner;
   int  algorithm;
   int  resumeid;
   int  pricetype;
   double  prevobj;
   double  objtol;
   param_info  oldinfo;
} iter_info;

void
    ILLsimplex_init_lpinfo (lpinfo *lp),
    ILLsimplex_free_lpinfo (lpinfo *lp),
    ILLsimplex_load_lpinfo (ILLlpdata *qslp, lpinfo *lp),
    ILLsimplex_set_bound (lpinfo *lp, double objbound, int sense);

int
    ILLsimplex_retest_psolution (lpinfo *lp, price_info *p, int phase, feas_info *fs),
    ILLsimplex_retest_dsolution (lpinfo *lp, price_info *p, int phase, feas_info *fs),
    ILLsimplex_solution (lpinfo *lp, double *xz, double *piz, double *dz, double *objval),
    ILLsimplex_infcertificate (lpinfo *lp, double *pi),
    ILLsimplex (lpinfo *lp, int algorithm, ILLlp_basis *B, price_info *pinf,
                int *sol_status, int sdisplay),
    ILLsimplex_pivotin (lpinfo *lp, price_info *pinf, int rcnt, int *rlist,
                        int pivot_opt, int *basis_mod);

#endif /* __SIMPLEX_H */
