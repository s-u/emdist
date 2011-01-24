/*  $RCSfile: ratio.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __RATIO_H
#define __RATIO_H

/* defs for phase I ratio test */
#define BBOUND    1
#define BATOLOWER 2
#define BATOUPPER 3
#define BBTOLOWER 4
#define BBTOUPPER 5
#define BSKIP     6

/* result of ratio test */
#define RATIO_UNBOUNDED 1
#define RATIO_NOBCHANGE 2
#define RATIO_BCHANGE   3
#define RATIO_FAILED    4
#define RATIO_NEGATIVE  5

typedef struct ratio_res {
    double tz;
    int eindex; 
    int lindex;
    int lvstat;
    int ratio_stat;
    int boundch;
    int coeffch; 
    double lbound; 
    double ecoeff; 
    double pivotval; 
} ratio_res; 

void
    ILLratio_pI_test (lpinfo *lp, int eindex, int dir, ratio_res *rs), 
    ILLratio_pII_test (lpinfo *lp, int eindex, int dir, ratio_res *rs), 
    ILLratio_dI_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs), 
    ILLratio_dII_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs),
    ILLratio_longdII_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs),
    ILLratio_pivotin_test (lpinfo *lp, int *rlist, int rcnt, ratio_res *rs);

#endif /* __RATIO_H */
