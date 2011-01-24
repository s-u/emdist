/* RCSINFO $Id: lib.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef ILL_LIB_H
#define ILL_LIB_H

#include "lpdefs.h"
#include "lpdata.h"
#include "price.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

#define ILL_LP_OPTIMAL           1
#define ILL_LP_NONOPTIMAL        2
#define ILL_LP_PRIMAL_FEASIBLE   3
#define ILL_LP_PRIMAL_INFEASIBLE 4
#define ILL_LP_PRIMAL_UNBOUNDED  5
#define ILL_LP_DUAL_FEASIBLE     6
#define ILL_LP_DUAL_INFEASIBLE   7
#define ILL_LP_DUAL_UNBOUNDED    8


/****************************************************************************/
/*                                                                          */
/*                               lib.c                                      */
/*                                                                          */
/****************************************************************************/


int
    ILLlib_optimize (lpinfo *lp, ILLlp_basis *B, price_info *pinf, int algo,
        int *status, int simplex_display),
    ILLlib_cache_solution (lpinfo *lp, ILLlp_cache *C),
    ILLlib_solution (lpinfo *lp, ILLlp_cache *C, double *val, double *x,
        double *pi, double *slack, double *rc),
    ILLlib_get_x (lpinfo *lp, ILLlp_cache *C, double *x),
    ILLlib_get_slack (lpinfo *lp, ILLlp_cache *C, double *slack),
    ILLlib_objval (lpinfo *lp, ILLlp_cache *C, double *val),
    ILLlib_newrow (lpinfo *lp, ILLlp_basis *B, double rhs, char sense,
        double range, const char *name),
    ILLlib_newrows (lpinfo *lp, ILLlp_basis *B, int num, double *rhs,
        char *sense, double *range, const char **names),
    ILLlib_addrow (lpinfo *lp, ILLlp_basis *B, int cnt, int *ind, double *val,
        double rhs, char sense, double range, const char *rowname),
    ILLlib_addrows (lpinfo *lp, ILLlp_basis *B, int num, int *rmatcnt,
        int *rmatbeg, int *rmatind, double *rmatval, double *rhs, char *sense,
        double *range, const char **names, int *nofactor),
    ILLlib_delrows (lpinfo *lp, ILLlp_basis *B, ILLlp_cache *C, int num,
        int *dellist, int *basis_ok, int *cache_ok),
    ILLlib_newcol (lpinfo *lp, ILLlp_basis *B, double obj, double lower,
        double upper, const char *name, int factorok),
    ILLlib_newcols (lpinfo *lp, ILLlp_basis *B, int num, double *obj,
        double *lower, double *upper, const char **names, int factorok),
    ILLlib_addcol (lpinfo *lp, ILLlp_basis *B, int cnt, int *ind, double *val,
        double obj, double lower, double upper, const char *name,
        int factorok),
    ILLlib_addcols (lpinfo *lp, ILLlp_basis *B, int num, int *cmatcnt,
       int *cmatbeg, int *cmatind, double *cmatval, double *obj,
       double *lower, double *upper, const char **names, int factorok),
    ILLlib_delcols (lpinfo *lp, ILLlp_basis *B, int num, int *dellist,
       int *basis_ok),
    ILLlib_chgcoef (lpinfo *lp, int rowindex, int colindex, double coef),
    ILLlib_chgsense (lpinfo *lp, int num, int *rowlist, char *sense),
    ILLlib_getrows (lpinfo *lp, int num, int *rowlist, int **rowcnt,
        int **rowbeg, int **rowind, double **rowval, double **rhs,
        char **sense, char ***names),
    ILLlib_getcols (lpinfo *lp, int num, int *collist, int **colcnt,
        int **colbeg, int **colind, double **colval, double **obj,
        double **lower, double **upper, char ***names),
    ILLlib_getobj (lpinfo *lp, double *obj),
    ILLlib_chgobj (lpinfo *lp, int indx, double coef),
    ILLlib_getrhs (lpinfo *lp, double *rhs),
    ILLlib_chgrhs (lpinfo *lp, int indx, double coef),
    ILLlib_getintflags (lpinfo *lp, int *intflags),
    ILLlib_rownames (lpinfo *lp, char **rownames),
    ILLlib_colnames (lpinfo *lp, char **colnames),
    ILLlib_colindex (lpinfo *lp, const char *name, int *colindex),
    ILLlib_rowindex (lpinfo *lp, const char *name, int *rowindex),
    ILLlib_chgbnd  (lpinfo *lp, int indx, char lu, double bnd), 
    ILLlib_chgbnds (lpinfo *lp, int cnt, int *indx, char *lu, double *bnd),
    ILLlib_getbnd (lpinfo *lp, int indx, char lu, double *bnd),
    ILLlib_getbnds (lpinfo *lp, double *lower, double *upper),
    ILLlib_strongbranch (lpinfo *lp, price_info *pinf, int *candidatelist,
        int ncand, double *xlist, double *downpen, double *uppen,
        int iterations, double objbound),
    ILLlib_getbasis (lpinfo *lp, char *cstat, char *rstat),
    ILLlib_loadbasis (ILLlp_basis *B, int nstruct, int nrows, char *cstat,
        char *rstat),
    ILLlib_readbasis (lpinfo *lp, ILLlp_basis *B, const char *fname),
    ILLlib_writebasis (lpinfo *lp, ILLlp_basis *B, const char *fname),
    ILLlib_getrownorms (lpinfo *lp, price_info *pinf, double *rownorms),
    ILLlib_loadrownorms (lpinfo *lp, price_info *pinf, double *rownorms),
    ILLlib_recompute_rownorms (lpinfo *lp, price_info *pinf),
    ILLlib_iter (lpinfo *lp),
    ILLlib_print_x (FILE *fd, lpinfo *lp, ILLlp_cache *C, double *x,
        int nonZerosOnly), 
	ILLwrite_lp_file (ILLlpdata *lp, FILE *eout, qserror_collector *c);


extern int ILLlib_findName(ILLlpdata *qslp, int forRow, const char *name,
                            int id, char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int
    ILLpresolve_add_logicals(ILLlpdata *lp);


/****************************************************************************/
/*                                                                          */
/*                            binary.c                                      */
/*                                                                          */
/****************************************************************************/

int
    ILLmip_binary_dfs (lpinfo *lp);

#endif /* ILL_LIB_H */

