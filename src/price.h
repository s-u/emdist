/*  $RCSfile: price.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __PRICE_H
#define __PRICE_H

#include "dstruct.h"

/*define as > 0 if heap is to be used */
#define USEHEAP 1

/*result of pricing */
#define PRICE_OPTIMAL 1
#define PRICE_NONOPTIMAL  2

/*type of pricing */
#define ROW_PRICING 1
#define COL_PRICING 2

typedef struct price_res{
   int eindex;
   int dir;
   int lindex;
   int lvstat;
   int price_stat;
   double dinfeas;
   double pinfeas;
} price_res;

int
    ILLprice_test_for_heap (lpinfo *lp, price_info *pinf, int nkeys,
                            double *keylist, int algo, int upd),
    ILLprice_build_heap (price_info *pinf, int nkeys, double *keylist),
    ILLprice_build_pricing_info (lpinfo *lp, price_info *pinf, int phase),
    ILLprice_update_pricing_info (lpinfo *lp, price_info *pinf, int phase,
                                svector *wz, int eindex, int lindex, double y),
    ILLprice_get_price (price_info *p, int phase),
    ILLprice_build_mpartial_info (lpinfo *lp, price_info *pinf, int pricetype),
    ILLprice_build_pdevex_norms (lpinfo *lp, p_devex_info *pdinfo, int reinit),
    ILLprice_update_pdevex_norms (lpinfo *lp, p_devex_info *pdinfo, int eindex,
                                  double yl),
    ILLprice_build_psteep_norms (lpinfo *lp, p_steep_info *psinfo),
    ILLprice_build_ddevex_norms (lpinfo *lp, d_devex_info *ddinfo, int reinit),
    ILLprice_update_ddevex_norms (lpinfo *lp, d_devex_info *ddinfo, int eindex,
                                  double yl),
    ILLprice_build_dsteep_norms (lpinfo *lp, d_steep_info *dsinfo),
    ILLprice_get_dsteep_norms (lpinfo *lp, int count, int *rowind, double *norms),
    ILLprice_get_rownorms (lpinfo *lp, price_info *pinf, double *rnorms),
    ILLprice_get_colnorms (lpinfo *lp, price_info *pinf, double *cnorms),
    ILLprice_get_newnorms (lpinfo *lp, int nelems, double *norms, int *matcnt,
                        int *matbeg, int *matind, double *matval, int option),
    ILLprice_get_new_rownorms (lpinfo *lp, int newrows, double *rnorms,
                    int *rmatcnt, int *rmatbeg, int *rmatind, double *rmatval),
    ILLprice_get_new_colnorms (lpinfo *lp, int newrows, double *rnorms,
                        int *matcnt, int *matbeg, int *matind, double *matval),
    ILLprice_load_rownorms (lpinfo *lp, double *rnorms, price_info *pinf),
    ILLprice_load_colnorms (lpinfo *lp, double *cnorms, price_info *pinf);


void
    ILLprice_free_heap (price_info *pinf),
    ILLprice_init_pricing_info (price_info *pinf),
    ILLprice_free_pricing_info (price_info *pinf),
    ILLprice_free_mpartial_info (mpart_info *p),
    ILLprice_init_mpartial_price (lpinfo *lp, price_info *pinf, int phase, int pricetype),
    ILLprice_update_mpartial_price (lpinfo *lp, price_info *pinf, int phase, int pricetype),
    ILLprice_delete_onempart_price (lpinfo *lp, price_info *pinf, int indx, int pricetype),
    ILLprice_mpartial_group (lpinfo *lp, mpart_info *p, int phase, int g, int pricetype),
    ILLprice_column (lpinfo *lp, int ix, int phase, price_res *pr),
    ILLprice_row (lpinfo *lp, int ix, int phase, price_res *pr),
    ILLprice_update_psteep_norms (lpinfo *lp, p_steep_info *psinfo,
                                   svector *wz, int eindex, double yl),
    ILLprice_update_dsteep_norms (lpinfo *lp, d_steep_info *dsinfo,
                                   svector *wz, int lindex, double yl),
    ILLprice_compute_dual_inf (lpinfo *lp, price_info *p, int *ix, int icnt, int phase),
    ILLprice_primal (lpinfo *lp, price_info *pinf, price_res *pr, int phase),
    ILLprice_compute_primal_inf (lpinfo *lp, price_info *p, int *ix, int icnt, int phase),
    ILLprice_dual (lpinfo *lp, price_info *pinf, int phase, price_res *pr);
 
void
    test_dsteep_norms (lpinfo *lp, price_info *p);

#endif /* __PRICE_H */
