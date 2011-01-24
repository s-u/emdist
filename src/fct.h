/* RCSINFO $Id: fct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#define CNT_YNZ           1 /* nz in entering columns */
#define CNT_ZNZ           2 /* nz in ith row of B^{-1}, ie z_i */
#define CNT_ZANZ          3 /* nz in ith row of B^{-1}, ie z_i */
#define CNT_PINZ          4 /* nz in phase II pi (solve) */
#define CNT_P1PINZ        5 /* nz in phase I pi (solve) */
#define CNT_UPNZ          6 /* nz in ftran_updates */
#define CNT_PPHASE1ITER   7 /* primal phase I iterations */
#define CNT_PPHASE2ITER   8
#define CNT_DPHASE1ITER   9 /* dual phase I iterations */
#define CNT_DPHASE2ITER   10
#define CNT_PIPIV         11
#define CNT_PIIPIV        12
#define CNT_DIPIV         13
#define CNT_DIIPIV        14
#define CNT_YRAVG         15
#define CNT_ZARAVG        16

#define ROW_PIVOT         0
#define COL_PIVOT         1

int
    ILLfct_compute_zA (lpinfo *lp, svector *z, svector *zA),
    ILLfct_compute_wz (lpinfo *lp, double *wz),
    ILLfct_adjust_viol_bounds (lpinfo *lp),
    ILLfct_perturb_bounds (lpinfo *lp),
    ILLfct_bound_shift (lpinfo *lp, int col, int bndtype, double newbnd),
    ILLfct_adjust_viol_coefs (lpinfo *lp),
    ILLfct_perturb_coefs (lpinfo *lp),
    ILLfct_coef_shift (lpinfo *lp, int col, double newcoef),
    ILLfct_test_pivot (lpinfo *lp, int indx, int indxtype, double piv_val);

void
    ILLfct_load_workvector (lpinfo *lp, svector *s),
    ILLfct_zero_workvector (lpinfo *lp),
    ILLfct_set_variable_type (lpinfo *lp),
    ILLfct_compute_pobj (lpinfo *lp),
    ILLfct_compute_dobj (lpinfo *lp),
    ILLfct_compute_xbz (lpinfo *lp),
    ILLfct_compute_piz (lpinfo *lp),
    ILLfct_compute_phaseI_xbz (lpinfo *lp),
    ILLfct_compute_phaseI_piz (lpinfo *lp),
    ILLfct_compute_dz (lpinfo *lp),
    ILLfct_compute_phaseI_dz (lpinfo *lp),
    ILLfct_compute_yz (lpinfo *lp, svector *yz, svector *updz, int ecol),
    ILLfct_compute_zz (lpinfo *lp, svector *zz, int lindex),
    ILLfct_compute_psteep_upv (lpinfo *lp, svector *swz),
    ILLfct_compute_dsteep_upv (lpinfo *lp, svector *swz),
    ILLfct_update_basis_info (lpinfo *lp, int eindex, int lindex, int lvstat),
    ILLfct_update_xz (lpinfo *lp, double tz, int eindex, int lindex),
    ILLfct_update_piz (lpinfo *lp, double alpha),
    ILLfct_update_pIpiz (lpinfo *lp, svector *z, double alpha),
    ILLfct_update_dz (lpinfo *lp, int eindex, double alpha),
    ILLfct_update_pIdz (lpinfo *lp, svector *zA, int eindex, double alpha),
    ILLfct_unroll_bound_change (lpinfo *lp),
    ILLfct_unroll_coef_change (lpinfo *lp),
    ILLfct_check_pfeasible (lpinfo *lp, feas_info *fs, double ftol),
    ILLfct_check_pIpfeasible (lpinfo *lp, feas_info *fs, double ftol),
    ILLfct_check_dfeasible (lpinfo *lp, feas_info *fs, double ftol),
    ILLfct_check_pIdfeasible (lpinfo *lp, feas_info *fs, double ftol),
    ILLfct_dual_adjust (lpinfo *lp, double ftol),
    ILLfct_dphaseI_simple_update (lpinfo *lp, double ftol),
    ILLfct_set_status_values (lpinfo *lp, int pstatus, int dstatus, int ptype,
                              int dtype),
    ILLfct_init_counts (lpinfo *lp),
    ILLfct_update_counts (lpinfo *lp, int f, int upi, double upd),
    ILLfct_print_counts (lpinfo *lp),
    ILLfct_check_pIdfeasible (lpinfo *lp, feas_info *fs, double ftol),
    ILLfct_update_pfeas (lpinfo *lp, int lindex, svector *srhs),
    ILLfct_compute_ppIzz (lpinfo *lp, svector *srhs, svector *ssoln),
    ILLfct_update_ppI_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                         svector *ssoln, int eindex, int lindex, double alpha),
    ILLfct_update_dfeas (lpinfo *lp, int eindex, svector *srhs),
    ILLfct_compute_dpIy (lpinfo *lp, svector *srhs, svector *ssoln),
    ILLfct_update_dpI_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                         svector *ssoln, int lindex, double alpha),
    ILLfct_update_dIIfeas (lpinfo *lp, int eindex, svector *srhs),
    ILLfct_compute_dpIIy (lpinfo *lp, svector *srhs, svector *ssoln),
    ILLfct_update_dpII_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                               svector *ssoln, int eindex, int lindex, double eval, double alpha);

void
    fct_test_workvector (lpinfo *lp),
    fct_test_pfeasible (lpinfo *lp),
    fct_test_dfeasible (lpinfo *lp),
    fct_test_pI_x (lpinfo *lp, price_info *p),
    fct_test_pII_x (lpinfo *lp, price_info *p),
    fct_test_pI_pi_dz (lpinfo *lp, price_info *p),
    fct_test_pII_pi_dz (lpinfo *lp, price_info *p);

#endif /* __FUNCTIONS_H */
