/* RCS_INFO = "$RCSfile: lib.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

/****************************************************************************/
/*                                                                          */
/*               Interface Routines to Core LP Solver                       */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLlib_optimize (lpinfo *lp, ILLlp_basis *B, price_info *pinf,    */
/*            int algo, int *status, int simplex_display)                   */
/*    int ILLlib_cache_solution (lpinfo *lp, ILLlp_cache *C)                */
/*    int ILLlib_solution (lpinfo *lp, ILLlp_cache *C, double *val,         */
/*            double *x, double *pi, double *slack, double *rc)             */
/*    int ILLlib_get_x (lpinfo *lp, ILLlp_cache *C, double *x)              */
/*    int ILLlib_get_slack (lpinfo *lp, ILLlp_cache *C, double *slack)      */
/*    int ILLlib_objval (lpinfo *lp, ILLlp_cache *C, double *val)           */
/*    int ILLlib_newrow (lpinfo *lp, ILLlp_basis *B, double rhs,            */
/*            char sense, double range, const char *name)                   */
/*        -range can specify a rangeval for the row (if sense is not 'R',   */
/*         then range is ignored); it should be 0 if no range is needed;    */
/*         if sense is 'R' but no rangeval array exists for the LP, the     */
/*         array will be allocated and initialized.                         */
/*    int ILLlib_newrows (lpinfo *lp, ILLlp_basis *B, int num, double *rhs, */
/*            char *sense, double *range, const char **names)               */
/*        -range is an array specifying the rangevals for the rows; range   */
/*         should be NULL if no rangevals are needed.                       */
/*    int ILLlib_addrow (lpinfo *lp, ILLlp_basis *B, int cnt, int *ind,     */
/*            double *val, double rhs, char sense, double range,            */
/*            const char *name)                                             */
/*    int ILLlib_addrows (lpinfo *lp, ILLlp_basis *B, int num,              */
/*            int *rmatcnt, int *rmatbeg, int *rmatind, double *rmatval,    */
/*            double *rhs, char *sense, double *range, const char **names,  */
/*            int *factorok)                                                */
/*    int ILLlib_delrows (lpinfo *lp, ILLlp_basis *B,                       */
/*            int num, int *dellist, int *basis_ok)                         */
/*    int ILLlib_newcol (lpinfo *lp, ILLlp_basis *B,                        */
/*            double obj, double lower, double upper, const char *name,     */
/*            int factorok)                                                 */
/*    int ILLlib_newcols (lpinfo *lp, ILLlp_basis *B,                       */
/*            int num, double *obj, double *lower, double *upper,           */
/*            const char **names, int factorok)                             */
/*    int ILLlib_addcol (lpinfo *lp, ILLlp_basis *B,                        */
/*            int cnt, int *ind, double *val, double obj, double lower,     */
/*            double upper, const char *name, int factorok)                 */
/*    int ILLlib_addcols (lpinfo *lp, ILLlp_basis *B,                       */
/*            int num, int *cmatcnt, int *cmatbeg, int *cmatind,            */
/*            double *cmatval, double *obj, double *lower, double *upper,   */
/*            const char **names, int factorok)                             */
/*    int ILLlib_delcols (lpinfo *lp, ILLlp_basis *B, int num, int *dellist */
/*            int *basis_ok)                                                */
/*    int ILLlib_chgcoef (lpinfo *lp, int rowindex, int colindex,           */
/*            double coef)                                                  */
/*    int ILLlib_chgsense (lpinfo *lp, int num, int *rowlist, char *sense)  */
/*    int ILLlib_getrows (lpinfo *lp, int num, int *rowlist, int **rowcnt,  */
/*            int **rowbeg, int **rowind, double **rowval, double **rhs,    */
/*            char **sense, char ***names)                                  */
/*    int ILLlib_getcols (lpinfo *lp, int num, int *collist, int **colcnt,  */
/*            int **colbeg, int **colind, double **colval, double **obj,    */
/*            double **lower, double **upper, char ***names)                */
/*    int ILLlib_getobj (lpinfo *lp, double *obj)                           */
/*    int ILLlib_chgobj (lpinfo *lp, int indx, double coef)                 */
/*    int ILLlib_getrhs (lpinfo *lp, double *rhs)                           */
/*    int ILLlib_chgrhs (lpinfo *lp, int indx, double coef)                 */
/*    int ILLlib_getintflags (lpinfo *lp, int *intflags)                    */
/*    int ILLlib_rownames (lpinfo *lp, char **rownames)                     */
/*    int ILLlib_colnames (lpinfo *lp, char **colnames)                     */
/*    int ILLlib_colindex (lpinfo *lp, char *name, int *colindex)           */
/*    int ILLlib_rowindex (lpinfo *lp, char *name, int *rowindex)           */
/*    int ILLlib_chgbnd  (lpinfo *lp, int indx, char lu, double bnd)        */
/*    int ILLlib_chgbnds (lpinfo *lp, int cnt, int *indx, char *lu,         */
/*            double *bnd)                                                  */
/*    int ILLlib_getbnd (lpinfo *lp, int indx, char lu, double *bnd)        */
/*    int ILLlib_getbnds (lpinfo *lp, double *lower, double *upper)         */
/*    int ILLlib_strongbranch (lpinfo *lp, price_info *pinf,                */
/*      int *candidatelist, int ncand, double *xlist, double *downpen,      */
/*      double *uppen, int iterations, double objbound)                     */
/*    int ILLlib_getbasis (lpinfo *lp, char *cstat, char *rstat)            */
/*    int ILLlib_loadbasis (ILLlp_basis *B, int nstruct, int nrows,         */
/*      char *cstat, char *rstat)                                           */
/*    int ILLlib_readbasis (lpinfo *lp, ILLlp_basis *B, char *fname)        */
/*    int ILLlib_writebasis (lpinfo *lp, const char *fname)                 */
/*    int ILLlib_getrownorms (lpinfo *lp, price_info *pinf,                 */
/*            double *rownorms)                                             */
/*    int ILLlib_loadrownorms (lpinfo *lp, price_info *pinf,                */
/*            double *rownorms)                                             */
/*    int ILLlib_recompute_rownorms (lpinfo *lp, price_info *pinf)          */
/*    int ILLlib_print_x (FILE *fd, lpinfo *lp, ILLlp_cache *C, double *x,  */
/*            int nonZerosOnly)                                             */
/*    int ILLlib_print_x (lpinfo *lp, ILLlp_cache *C)                       */
/*    int ILLlib_iter (lpinfo *lp)                                          */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"
#include "simplex.h"
#include "price.h"
#include "basis.h"
#include "lib.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lp.h"
#include "mps.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static void
    check_pinf (price_info *pinf, int *it_exists);

static int
    matrix_addrow (ILLmatrix *A, int rowcnt, int *rowind, double *rowval),
    matrix_addrow_end (ILLmatrix *A, int row, int rowcnt, int *rowind,
        double *rowval),
    matrix_addcoef (lpinfo *lp, ILLmatrix *A, int row, int col, double val),
    matrix_addcol (ILLmatrix *A, int colcnt, int *colind, double *colval),
    delcols_work (lpinfo *lp, char *colmark),
    reset_colindex (lpinfo *lp),
    reset_rowindex (lpinfo *lp);

int ILLlib_optimize (lpinfo *lp, ILLlp_basis *B, price_info *pinf, int algo,
        int *status, int simplex_display)
{
    int rval = 0;
    int sol_status;

    if (status) *status = QS_LP_UNSOLVED;

    /* ILLprice_free_pricing_info (pinf); */   /* Should be removed later */

    rval = ILLsimplex (lp, algo, B, pinf, &sol_status, simplex_display);
    ILL_CLEANUP_IF (rval);

    if (status) *status = sol_status;

CLEANUP:

    if (rval == E_SIMPLEX_ERROR) {
        FILE *eout = (FILE *) NULL;
        int tval;

        printf ("write bad lp to error.lp\n"); fflush (stdout);
        eout = fopen ("error.lp", "w");
        if (!eout) {
            fprintf (stderr, "could not open file to write bad lp\n");
        } else {
            tval = ILLwrite_lp(lp->O, NULL); 
            if (tval) {
                fprintf (stderr, "error while writing bad lp\n");
            }
            fclose (eout);
        }
        
        printf ("write bad basis to error.bas\n"); fflush (stdout);
        tval = ILLlib_writebasis (lp, (ILLlp_basis *) NULL, "error.bas");
        if (tval) {
            fprintf (stderr, "error while writing bad basis\n");
        }
    }

    ILL_RETURN (rval, "ILLlib_optimize");
}

int ILLlib_cache_solution (lpinfo *lp, ILLlp_cache *C)
{
    int rval = 0;

    if (C) {
        if (C->nstruct != lp->O->nstruct || C->nrows != lp->O->nrows) {
            fprintf (stderr, "lp_cache does not match size of lp\n");
            rval = 1; ILL_CLEANUP;
        }
        rval = ILLlib_solution (lp, (ILLlp_cache *) NULL, &(C->val), C->x,
                                    C->pi, C->slack, C->rc);
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_cache_solution");
}

int ILLlib_solution (lpinfo *lp, ILLlp_cache *C, double *val, double *x,
        double *pi, double *slack, double *rc)
{
    int i, rval = 0;
    double *tempx  = (double *) NULL;
    double *temprc = (double *) NULL;
    int ncols = lp->O->ncols;
    int nrows = lp->O->nrows;
    int nstruct = lp->O->nstruct;
    ILLlpdata *qslp = lp->O;

    if (C) {
        if (C->nrows != nrows || C->nstruct != nstruct) {
            fprintf (stderr, "cache mismatch in ILLlib_solution\n");
            rval = 0; ILL_CLEANUP;
        }
        if (val) {
            *val = C->val;
        }
        if (x) {
            for (i = 0; i < nstruct; i++) {
                x[i] = C->x[i];
            }
        }
        if (pi) {
            for (i = 0; i < nrows; i++) {
                pi[i] = C->pi[i];
            }
        }
        if (slack) {
            for (i = 0; i < nrows; i++) {
                slack[i] = C->slack[i];
            }
        }
        if (rc) {
            for (i = 0; i < nstruct; i++) {
                rc[i] = C->rc[i];
            }
        }
    } else {
        if (x || slack) {
            ILL_SAFE_MALLOC (tempx, ncols, double);
        }

        if (rc) {
            ILL_SAFE_MALLOC (temprc, ncols, double);
        }

        rval = ILLsimplex_solution (lp, tempx, pi, temprc, val);
        ILL_CLEANUP_IF (rval);

        if (x) {
            for (i = 0; i < nstruct; i++) {
                x[i] = tempx[qslp->structmap[i]];
            }
        }
        if (slack) {
            for (i = 0; i < nrows; i++) {
                slack[i] = tempx[qslp->rowmap[i]];
            }
        }

        if (rc) {
            for (i = 0; i < nstruct; i++) {
                rc[i] = temprc[qslp->structmap[i]];
            }
        }

        if (val) {
            if (lp->O->objsense == ILL_MAX) (*val) *= -1.0;
        }
    }

CLEANUP:

    ILL_IFFREE (tempx, double);
    ILL_IFFREE (temprc, double);

    ILL_RETURN (rval, "ILLlib_solution");
}

int ILLlib_get_x (lpinfo *lp, ILLlp_cache *C, double *x)
{
    int rval = 0;

    rval = ILLlib_solution (lp, C, (double *) NULL, x, (double *) NULL,
                                  (double *) NULL, (double *) NULL);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_get_x");
}

int ILLlib_get_slack (lpinfo *lp, ILLlp_cache *C, double *slack)
{
    int rval = 0;

    rval = ILLlib_solution (lp, C, (double *) NULL, (double *) NULL,
                                  (double *) NULL, slack, (double *) NULL);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_get_slack");
}


int ILLlib_objval (lpinfo *lp, ILLlp_cache *C, double *val)
{
    int rval = 0;

    if (lp->basisstat.optimal) {
        rval = ILLlib_solution (lp, C, val, (double *) NULL, (double *) NULL,
                                           (double *) NULL, (double *) NULL);
        ILL_CLEANUP_IF (rval);
    } else {
        *val = lp->dobjval;  /* Ask Sanjeeb */
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_objval");
}

int ILLlib_chgbnd  (lpinfo *lp, int indx, char lu, double bnd)
{
    int rval = 0;
    int col;

    if (!lp) {
        fprintf (stderr, "ILLlib_chgbnd called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (indx < 0 || indx > lp->O->nstruct) {
        fprintf (stderr, "ILLlib_chgbnd called with bad indx: %d\n", indx);
        rval = 1; ILL_CLEANUP;
    }

    if (lp->O->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (lp->O->sinfo);
        ILL_IFFREE (lp->O->sinfo, ILLlp_sinfo);
    }

    col = lp->O->structmap[indx];

    switch (lu) {
    case 'L':
        lp->O->lower[col] = bnd; break;
    case 'U':
        lp->O->upper[col] = bnd; break;
    case 'B':
        lp->O->lower[col] = bnd; lp->O->upper[col] = bnd; break;
    default:
        fprintf (stderr, "ILLlib_chgbnd called with lu: %c\n", lu);
        rval = 1; ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgbnd");
}

int ILLlib_chgbnds (lpinfo *lp, int cnt, int *indx, char *lu, double *bnd)
{
    int rval = 0;
    int i;

    for (i = 0; i < cnt; i++) {
        rval = ILLlib_chgbnd (lp, indx[i], lu[i], bnd[i]);
        if (rval) ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgbnds");
}

int ILLlib_getbnd (lpinfo *lp, int indx, char lu, double *bnd)
{
    int rval = 0;
    int col;

    if (!lp) {
        fprintf (stderr, "ILLlib_getbnd called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (indx < 0 || indx > lp->O->nstruct) {
        fprintf (stderr, "ILLlib_getbnd called with bad indx: %d\n", indx);
        rval = 1; ILL_CLEANUP;
    }

    col = lp->O->structmap[indx];

    switch (lu) {
    case 'L':
        *bnd = lp->O->lower[col]; break;
    case 'U':
        *bnd = lp->O->upper[col]; break;
    default:
        fprintf (stderr, "ILLlib_getbnd called with lu: %c\n", lu);
        rval = 1; ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getbnd");
}

int ILLlib_getbnds (lpinfo *lp, double *lower, double *upper)
{
    int rval = 0;
    ILLlpdata *qslp;
    int nstruct;
    int j, col;

    if (!lp) {
        fprintf (stderr, "ILLlib_getbnd called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nstruct = qslp->nstruct;

    for (j = 0; j < nstruct; j++) {
        col = qslp->structmap[j];
        if (lower) lower[j] = qslp->lower[col];
        if (upper) upper[j] = qslp->upper[col];
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getbnds");
}

int ILLlib_strongbranch (lpinfo *lp, price_info *pinf, int *candidatelist,
        int ncand, double *xlist, double *downpen, double *uppen,
        int iterations, double objbound)
{
    int rval = 0;
    int i, k, status, have_norms;
    int olditer = lp->maxiter;
    int nstruct = lp->O->nstruct;
    int nrows = lp->O->nrows;
    double *rnorms = (double *) NULL;
    double *myx = (double *) NULL;
    double xi, t, oldbnd;
    price_info lpinf;
    ILLlp_basis B, origB;

    ILLlp_basis_init (&B);
    ILLlp_basis_init (&origB);
    ILLprice_init_pricing_info (&lpinf);
    lpinf.dI_price  = QS_PRICE_DSTEEP;
    lpinf.dII_price = QS_PRICE_DSTEEP;

    if (xlist == (double *) NULL) {
        ILL_SAFE_MALLOC (myx, nstruct, double);
        rval = ILLlib_get_x (lp, (ILLlp_cache *) NULL, myx);
        ILL_CLEANUP_IF (rval);
    } 

    rval = ILLlp_basis_alloc (&origB, nstruct, nrows);
    ILL_CLEANUP_IF (rval);

    rval = ILLlib_getbasis (lp, origB.cstat, origB.rstat);
    ILL_CLEANUP_IF (rval);

    check_pinf (pinf, &have_norms);
    if (have_norms == 0) {
        ILL_SAFE_MALLOC (origB.rownorms, nrows, double);
        rval = ILLlib_getrownorms (lp, pinf, origB.rownorms);
        ILL_CLEANUP_IF (rval);
    } else {
        lp->basisid = -1;
        rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, &lpinf, DUAL_SIMPLEX,
                               &status, 0);
        ILL_CLEANUP_IF (rval);
    }

    rval = ILLlp_basis_alloc (&B, nstruct, nrows);  /* Note: B and orgiB may */
                                                    /* differ.               */
    ILL_CLEANUP_IF (rval);

    rval = ILLlib_getbasis (lp, B.cstat, B.rstat);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (B.rownorms, nrows, double);

    if (have_norms == 0) {
        rval = ILLlib_getrownorms (lp, pinf, B.rownorms);
        ILL_CLEANUP_IF (rval);
    } else {
        rval = ILLlib_getrownorms (lp, &lpinf, B.rownorms);
        ILL_CLEANUP_IF (rval);
    }

    lp->maxiter = iterations;

    for (i = 0; i < ncand; i++) {
        k = candidatelist[i];
        rval = ILLlib_getbnd (lp, k, 'U', &oldbnd);
        ILL_CLEANUP_IF (rval);
        if (xlist) xi = xlist[i];
        else       xi = myx[k];
        t = ILLutil_our_floor (xi);
        if (t >= -0.1 && t <= 0.1) {
            t = 0.0;   /* Helps to have a hard 0. */
        }

        rval = ILLlib_chgbnd (lp, k, 'U', t);
        ILL_CLEANUP_IF (rval);

        rval = ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX,
                               &status, 0);
        ILL_CLEANUP_IF (rval);

        downpen[i] = lp->dobjval;
        rval = ILLlib_chgbnd (lp, k, 'U', oldbnd);
        ILL_CLEANUP_IF (rval);

        rval = ILLlib_getbnd (lp, k, 'L', &oldbnd);
        ILL_CLEANUP_IF (rval);
        t = ILLutil_our_ceil (xi);
        if (t >= 0.9 && t <= 1.1) {
            t = 1.0;   /* Helps to have a hard 1. */
        }
        rval = ILLlib_chgbnd (lp, k, 'L', t);
        ILL_CLEANUP_IF (rval);

        rval = ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX,
                               &status, 0);
        ILL_CLEANUP_IF (rval);

        uppen[i] = lp->dobjval;
        rval = ILLlib_chgbnd (lp, k, 'L', oldbnd);
        ILL_CLEANUP_IF (rval);
    }

    if (lp->O->objsense == ILL_MAX) {

    } else {
        for (i = 0; i < ncand; i++) {
            if (downpen[i] > objbound) downpen[i] = objbound;
            if (  uppen[i] > objbound) uppen[i]   = objbound;
        }
    }

    /* Restore the old optimal solution */

    lp->maxiter = olditer;
    rval = ILLlib_optimize (lp, &origB, pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    lp->maxiter = olditer;
    ILLprice_free_pricing_info (&lpinf);
    ILLlp_basis_free (&B);
    ILLlp_basis_free (&origB);
    ILL_IFFREE (rnorms, double);
    if (xlist == (double *) NULL) {
        ILL_IFFREE (myx, double);
    }
    ILL_RETURN (rval, "ILLlib_strongbranch");
}

#define EXTRA_ROWS (100)
#define EXTRA_COLS (100)
#define EXTRA_MAT  (1000)

int ILLlib_newrow (lpinfo *lp, ILLlp_basis *B, double rhs, char sense,
        double range, const char *name)
{
    int rval = 0;

    rval = ILLlib_addrow (lp, B, 0, (int *) NULL, (double *) NULL, rhs,
                         sense, range, name);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_newrow");
}

int ILLlib_newrows (lpinfo *lp, ILLlp_basis *B, int num, double *rhs,
        char *sense, double *range, const char **names)
{
    int rval = 0;
    int *rmatcnt = (int *) NULL;
    int *rmatbeg = (int *) NULL;
    int i;

    if (!num) ILL_CLEANUP;

    ILL_SAFE_MALLOC (rmatcnt, num, int);
 
    ILL_SAFE_MALLOC (rmatbeg, num, int);
    
    for (i = 0; i < num; i++) {
        rmatcnt[i] = 0;
        rmatbeg[i] = 0;
    }

    rval = ILLlib_addrows (lp, B, num, rmatcnt, rmatbeg, (int *) NULL,
                          (double *) NULL, rhs, sense, range, names,
                          (int *) NULL);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_IFFREE (rmatcnt, int);
    ILL_IFFREE (rmatbeg, int);

    ILL_RETURN (rval, "ILLlib_newrows");
}

int ILLlib_addrows (lpinfo *lp, ILLlp_basis *B, int num, int *rmatcnt,
        int *rmatbeg, int *rmatind, double *rmatval, double *rhs, char *sense,
        double *range, const char **names, int *factorok)
{
    int rval = 0;
    int i, j, total, bsing;
    int *imap = (int *) NULL;
    int *bbeg = (int *) NULL;
    int *bcnt = (int *) NULL;
    int *bindi = (int *) NULL;
    int *rindi = (int *) NULL;
    int *jstat = (int *) NULL;
    double *bval = (double *) NULL;
    double rng;
    int badfactor = 0;

    if (B == (ILLlp_basis *) NULL || B->rownorms == (double *) NULL) {
        if (factorok) *factorok = 0;
    }

    if (B != (ILLlp_basis *) NULL) {
        ILL_IFFREE (B->colnorms, double);
    }

    if (B != (ILLlp_basis *) NULL && B->rownorms != (double *) NULL &&
        factorok != (int *) NULL && *factorok == 1) {
        int *structmap = lp->O->structmap;
        
        lp->matbeg = lp->O->A.matbeg;
        lp->matcnt = lp->O->A.matcnt;
        lp->matind = lp->O->A.matind;
        lp->matval = lp->O->A.matval;
  
        lp->nrows = lp->O->nrows;
        lp->ncols = lp->O->ncols;

        rval = ILLutil_reallocrus_count ((void **) &(B->rownorms),
                                 lp->O->nrows + num, sizeof (double));
        ILL_CLEANUP_IF (rval);

        ILL_SAFE_MALLOC (bcnt, num, int);
        ILL_SAFE_MALLOC (bbeg, num, int);
        ILL_SAFE_MALLOC (imap, lp->O->nstruct, int);

        ILL_SAFE_MALLOC (jstat, lp->ncols, int);

        for (i = 0; i < lp->ncols; i++) {
            jstat[i] = -1;
        }
        for (i = 0; i < lp->O->nstruct; i++) {
            jstat[structmap[i]] = i;
        }

        for (i = 0; i < lp->O->nstruct; i++) {
            imap[i] = -1;
        }
        for (i = 0; i < lp->nrows; i++) {
            if (jstat[lp->baz[i]] != -1) {
                imap[jstat[lp->baz[i]]] = i;
            }
        }
             
        for (i = 0, total = 0; i < num; i++) {
            bcnt[i] = 0;
            bbeg[i] = total;
            for (j = 0; j < rmatcnt[i]; j++) {
                if (imap[rmatind[rmatbeg[i]+j]] != -1) {
                    bcnt[i]++;
                    total++;
                }
            }
        }
        if (total) {
            ILL_SAFE_MALLOC (bindi, total, int);
            ILL_SAFE_MALLOC (bval, total, double);
        }
        for (i = 0, total = 0; i < num; i++) {
            for (j = 0; j < rmatcnt[i]; j++) {
                if (imap[rmatind[rmatbeg[i]+j]] != -1) {
                    bval[total] = rmatval[rmatbeg[i]+j];
                    bindi[total] = imap[rmatind[rmatbeg[i]+j]];
                    total++;
                }
            }
        }

        rval = ILLprice_get_new_rownorms(lp, num, B->rownorms + lp->O->nrows,
                                bcnt, bbeg, bindi, bval);
        ILL_CLEANUP_IF (rval);

        ILL_IFFREE (bcnt, int);
        ILL_IFFREE (bbeg, int);
        ILL_IFFREE (bindi, int);
        ILL_IFFREE (bval, double);
        ILL_IFFREE (imap, int);

        badfactor = 1;
    }

    for (i = 0; i < num; i++) {
        rng = (range != (double *) NULL ? range[i] : 0.0);
        if (names) { 
            rval = ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i], 
                  rmatval + rmatbeg[i], rhs[i], sense[i], rng, names[i]);
        } else {
            rval = ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i], 
                  rmatval + rmatbeg[i], rhs[i], sense[i], rng, (char *) NULL);
        }
        ILL_CLEANUP_IF (rval);
    }


    if (B != (ILLlp_basis *) NULL && B->rownorms != (double *) NULL &&
       (factorok != (int *) NULL && *factorok == 0)) {
        lp->matbeg  = lp->O->A.matbeg;
        lp->matcnt  = lp->O->A.matcnt;
        lp->matind  = lp->O->A.matind;
        lp->matval  = lp->O->A.matval;
        lp->nrows   = lp->O->nrows;
        lp->ncols   = lp->O->ncols;
        lp->bz      = lp->O->rhs;
        lp->nnbasic = lp->ncols - lp->nrows;

        rval = ILLbasis_load (lp, B);
        ILL_CLEANUP_IF (rval);

        if (lp->f) ILLfactor_free_factor_work (lp->f);

        rval = ILLbasis_factor (lp, &bsing);
        ILL_CLEANUP_IF (rval);
        *factorok = 1;

        rval = ILLutil_reallocrus_count ((void **) &(B->rownorms),
                                 lp->O->nrows , sizeof (double));
        ILL_CLEANUP_IF (rval);

        ILL_SAFE_MALLOC (rindi, lp->O->nrows /* num */, int);

        for (i = 0; i < num; i++) {
            rindi[i] = lp->O->nrows - num + i;
        }

        rval = ILLprice_get_dsteep_norms (lp, num, rindi,
                                  B->rownorms + lp->O->nrows - num);
        ILL_CLEANUP_IF (rval);
    }

    if (factorok != (int *) NULL && badfactor == 1) {
        *factorok = 0;
    }


CLEANUP:

    ILL_IFFREE (bcnt, int);
    ILL_IFFREE (bbeg, int);
    ILL_IFFREE (bindi, int);
    ILL_IFFREE (bval, double);
    ILL_IFFREE (imap, int);
    ILL_IFFREE (jstat, int);
    ILL_IFFREE (rindi, int);
    ILL_RETURN (rval, "ILLlib_addrows");
}

int ILLlib_addrow (lpinfo *lp, ILLlp_basis *B, int cnt, int *ind, double *val,
        double rhs, char sense, double range, const char *name)
{
    int rval = 0;
    ILLlpdata *qslp;
    ILLmatrix *A;
    int i, nrows, ncols; 
    char buf[ILL_namebufsize];
    int tind[1];
    double tval[1];
    int *tempind = (int *) NULL;
    int pind, hit;

    if (!lp) {
        fprintf (stderr, "ILLlib_addrow called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    A = &qslp->A;

    if (qslp->rA) {  /* After an addrow call, needs to be updated */
        ILLlp_rows_clear (qslp->rA);
        ILL_IFFREE (qslp->rA, ILLlp_rows);
    }

    if (qslp->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (qslp->sinfo);
        ILL_IFFREE (qslp->sinfo, ILLlp_sinfo);
    }

    nrows = qslp->nrows;
    ncols = qslp->ncols;

    /* If the row has a range, create the rangeval array if needed  */

    if (sense == 'R' && qslp->rangeval == (double *) NULL &&
                        qslp->rowsize > 0) {
        ILL_SAFE_MALLOC(qslp->rangeval, qslp->rowsize, double);
        for (i = 0; i < qslp->nrows; i++) {
            qslp->rangeval[i] = 0.0;
        }
    }

    /* Add the row to the row structures */

    if (qslp->rowsize < nrows + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(qslp->rhs),
                   qslp->rowsize + EXTRA_ROWS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->sense),
                   qslp->rowsize + EXTRA_ROWS, sizeof (char));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->rowmap),
                   qslp->rowsize + EXTRA_ROWS, sizeof (int));
        ILL_CLEANUP_IF (rval);

        if (qslp->rangeval || sense == 'R') {
            rval = ILLutil_reallocrus_count ((void **) &(qslp->rangeval),
                   qslp->rowsize + EXTRA_ROWS, sizeof (double));
            ILL_CLEANUP_IF (rval);
        }
         
	rval = ILLutil_reallocrus_count ((void **) &(qslp->rownames),
                   qslp->rowsize + EXTRA_ROWS, sizeof (char *));
	ILL_CLEANUP_IF (rval);
        qslp->rowsize += EXTRA_ROWS;
    }

    qslp->rhs[nrows]    = rhs;
    qslp->sense[nrows]  = sense;
    qslp->rowmap[nrows] = ncols;  /* this will be the new logical */
    if (qslp->rangeval) {
        if (sense == 'R') qslp->rangeval[nrows] = range;
        else              qslp->rangeval[nrows] = 0.0;
    }
    ILL_FAILtrue(qslp->rownames == NULL, "must always be non NULL"); 
    ILLlib_findName(qslp, 1 /*row*/, name, nrows, buf); 
    ILL_UTIL_STR(qslp->rownames[nrows], buf);
    ILLsymboltab_register(&qslp->rowtab, buf, qslp->nrows, &pind, &hit); 
    ILL_FAILfalse(hit == 0, "must be new"); 


    /* Add the logical variable to the column structures */

    if (qslp->colsize < ncols + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(qslp->lower),
                   qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->upper),
                   qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->obj),
                   qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        qslp->colsize += EXTRA_COLS;
    }

    qslp->obj[ncols] = 0.0;
    qslp->lower[ncols] = 0.0;
    if (sense == 'E') {
        qslp->upper[ncols] = 0.0;            /* Artificial */
    } else if (sense == 'R') {             
        qslp->upper[ncols] = range;          /* Range      */
    } else {
        qslp->upper[ncols] = ILL_MAXDOUBLE;  /* Slack      */
    }

    /* Add new row and new logical col to matrix */

    /* Need to map the structural indices to their proper place */
    
    if (cnt) {
        ILL_SAFE_MALLOC (tempind, cnt, int);
        for (i = 0; i < cnt; i++) {
            tempind[i] = qslp->structmap[ind[i]];
        }
    }

    rval = matrix_addrow (A, cnt, tempind, val);
    ILL_CLEANUP_IF (rval);

    tind[0] = nrows;
    if (sense == 'G' || sense == 'R') tval[0] = -1;
    else                              tval[0] = 1;

    rval = matrix_addcol (A, 1, tind, tval);
    ILL_CLEANUP_IF (rval);

    if (B != (ILLlp_basis *) NULL) {
        rval = ILLutil_reallocrus_count ((void **) &(B->rstat), nrows + 1,
					 sizeof (char));
        ILL_CLEANUP_IF (rval);
        B->rstat[nrows] = QS_ROW_BSTAT_BASIC;
    }

#if 0
    lp->basisid = -1;  /* To get optimizer to reload the basis */
#endif

    qslp->ncols++;
    qslp->nrows++;
    qslp->nzcount += (cnt+1);

    if (B != (ILLlp_basis *) NULL) {
        B->nrows++;
    }

CLEANUP:
    ILL_IFFREE (tempind, int);

    ILL_RETURN (rval, "ILLlib_addrow");
}

int ILLlib_delrows (lpinfo *lp, ILLlp_basis *B, ILLlp_cache *C, int num,
        int *dellist, int *basis_ok, int *cache_ok)
{
    int rval = 0;
    int i, j, k, nrows, ncols, nstruct, spot, dk, bok = 0, cok = 0;
    ILLlpdata *qslp;
    ILLmatrix *A;
    char *rowmark = (char *) NULL;
    char *colmark = (char *) NULL;
    int *newrowindex = (int *) NULL;	
    int *newcolindex = (int *) NULL;	
    int *ind, *beg, *cnt;
    double *val;

    if (!lp) {
        fprintf (stderr, "ILLlib_delrows called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (num <= 0) {
        if (basis_ok) *basis_ok = 1;
        if (cache_ok) *cache_ok = 1;
        ILL_CLEANUP;
    }

    if (basis_ok) *basis_ok = 0;
    if (cache_ok) *cache_ok = 0;

    qslp = lp->O;
    A = &qslp->A;

    if (qslp->rA) {  /* After a delrow call, needs to be updated */
        ILLlp_rows_clear (qslp->rA);
        ILL_IFFREE (qslp->rA, ILLlp_rows);
    }

    nrows = A->matrows;
    ncols = A->matcols;
    ind   = A->matind;
    beg   = A->matbeg;
    cnt   = A->matcnt;
    val   = A->matval;
    nstruct = qslp->nstruct;

    ILL_SAFE_MALLOC (rowmark, nrows, char);

    for (i = 0; i < nrows; i++) {
        rowmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
        rowmark[dellist[i]] = 1;
    }


    /* Try to update the basis */

    if (B) {
        bok = 1;
        cok = 1;
        for (i = 0; i < num; i++) {
            j = dellist[i];
            if (B->rstat[j] == QS_ROW_BSTAT_LOWER ||
                B->rstat[j] == QS_ROW_BSTAT_UPPER) {
                bok = 0;
                break;
            }
            if (C && C->pi[j] > 1e-6) {
/*
                printf ("XXXX: Postive pi (%f) at basic row\n", C->pi[j]);
                fflush (stdout);
*/
                cok = 0;
            }
        }
        if (bok == 1) {
            ILL_IFFREE (B->colnorms, double);
            if (B->rownorms) {
                for (i = 0, k = 0; i < nstruct; i++) {
                    if (B->cstat[i] == QS_COL_BSTAT_BASIC) k++;
                }
                for (i = 0, j = k; i < nrows; i++) {
                    if (B->rstat[i] == QS_ROW_BSTAT_BASIC) {
                        if (rowmark[i] == 0) {
                            B->rownorms[k++] = B->rownorms[j];
                        }
                        j++;
                    }
                }
                if (k != nrows - num) {
                    fprintf (stderr, "error in  ILLlib_delrows\n");
                    rval = 1; ILL_CLEANUP;
                }
            }

            for (i = 0, j = 0; i < nrows; i++) {
                if (rowmark[i] == 0) {
                    B->rstat[j++] = B->rstat[i];
                }
            }
            B->nrows = j;

            if (C && cok == 1) {
                for (i = 0, j = 0; i < nrows; i++) {
                    if (rowmark[i] == 0) {
                        C->pi[j] = C->pi[i];
                        C->slack[j++] = C->slack[i];
                    }
                }
                C->nrows = j;
                if (cache_ok) *cache_ok = 1;
            }
            if (basis_ok) *basis_ok = 1;
        }
    } 

    ILL_SAFE_MALLOC (newrowindex, nrows, int);


    /* Delete the marked rows */

    ILL_FAILtrue(qslp->rownames == NULL, "must always be non NULL"); 
    for (i = 0, j = 0; i < nrows; i++) {
        if (rowmark[i] == 0) {
            if (i != j) {
                qslp->rhs[j] = qslp->rhs[i];
                qslp->sense[j] = qslp->sense[i];
                if (qslp->rangeval) qslp->rangeval[j] = qslp->rangeval[i];
                if (qslp->rownames) qslp->rownames[j] = qslp->rownames[i];
            }
            newrowindex[i] = j++;
        } else {
            if (qslp->rownames) {
		rval = ILLsymboltab_delete(&qslp->rowtab, qslp->rownames[i]); 
                ILL_CLEANUP_IF(rval);
                ILL_IFFREE (qslp->rownames[i], char);
            }
        }
    }


    /* Delete the logicals */

    ILL_SAFE_MALLOC (colmark, ncols, char);

    for (i = 0; i < ncols; i++) {
        colmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
        colmark[qslp->rowmap[dellist[i]]] = 1;
    }

    rval = delcols_work (lp, colmark);
    ILL_CLEANUP_IF (rval);

    A->matcols -= num;
    qslp->ncols -= num;


    /* Pack the rowmap  */

    for (i = 0, j = 0; i < nrows; i++) {
        if (rowmark[i] == 0) {
            qslp->rowmap[j++] = qslp->rowmap[i];
        }
    }

    /* Remove the entries to deleted rows, and update the indices */

    for (i = 0; i < ncols - num; i++) {
        dk = 0;
        spot = beg[i];
        for (j = 0; j < cnt[i]; j++) {
            if (rowmark[ind[beg[i] + j]] == 1) {
                dk++;
            } else {
                val[spot] = val[beg[i] + j];
                ind[spot] = newrowindex[ind[beg[i] + j]];
                spot++;
            }
        }
        for (; spot < beg[i] + cnt[i]; spot++) {
            ind[spot] = -1;
        }

        cnt[i] -= dk;
        if (cnt[i] == 0) {
            ind[beg[i]] = 1;   /* we always mark the empty cols */
        }
    }

    A->matrows -= num;
    qslp->nrows -= num;

#if 0
    lp->basisid = -1;  /* To get optimizer to reload the basis */
#endif

CLEANUP:

    ILL_IFFREE (rowmark, char);
    ILL_IFFREE (colmark, char);
    ILL_IFFREE (newcolindex, int);
    ILL_IFFREE (newrowindex, int);
    ILL_RETURN (rval, "ILLlib_delrows");
}

int ILLlib_delcols (lpinfo *lp, ILLlp_basis *B, int num, int *dellist,
        int *basis_ok)
{
    int rval = 0;
    int i, j, bok, ncols;
    char *colmark = (char *) NULL;
    ILLlpdata *qslp;

    if (!lp) {
        fprintf (stderr, "ILLlib_delcols called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (basis_ok) *basis_ok = 0;

    if (num <= 0) {
        *basis_ok = 1;
        ILL_CLEANUP;
    }

    qslp = lp->O;
    ncols = qslp->A.matcols;

    if (qslp->rA) {  /* After a delcol call, needs to be updated */
        ILLlp_rows_clear (qslp->rA);
        ILL_IFFREE (qslp->rA, ILLlp_rows);
    }

    ILL_SAFE_MALLOC (colmark, ncols, char);

    for (i = 0; i < ncols; i++) {
        colmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
        colmark[qslp->structmap[dellist[i]]] = 1;
    }

    if (B) {
        B->nstruct -= num;
        bok = 1;
        for (i = 0; i < num; i++) {
            j = dellist[i];
            if (B->cstat[j] == QS_COL_BSTAT_BASIC) {
                bok = 0;
                printf ("BONG\n"); fflush (stdout);
                break;
            }
        }
        if (bok == 1) {
            ILL_IFFREE (B->colnorms, double);
            for (i = 0, j = 0; i < qslp->nstruct; i++) {
                if (colmark[qslp->structmap[i]] == 0) {
                    B->cstat[j++] = B->cstat[i];
                }
            }
            if (basis_ok) *basis_ok = 1;
        }
    }

    rval = delcols_work (lp, colmark);
    ILL_CLEANUP_IF (rval);


    qslp->A.matcols -= num;
    qslp->ncols -= num;
    qslp->nstruct -= num;

#if 0
    lp->basisid = -1;  /* To get optimizer to reload the basis */
#endif

CLEANUP:

    ILL_IFFREE (colmark, char);
    ILL_RETURN (rval, "ILLlib_delcols");
}

static int delcols_work (lpinfo *lp, char *colmark)
{
    int rval = 0;
    int i, j, k, nrows, ncols;
    ILLlpdata *qslp;
    ILLmatrix *A;
    int *newcolindex = (int *) NULL;
    int *ind, *beg, *cnt;

    /* Allows logicals to be deleted, to handle call from delcols. */

    qslp = lp->O;
    A = &qslp->A;
    nrows = A->matrows;
    ncols = A->matcols;
    ind   = A->matind;
    beg   = A->matbeg;
    cnt   = A->matcnt;

    ILL_SAFE_MALLOC (newcolindex, ncols, int);

    /* Delete the columns */

    for (i = 0, j = 0; i < ncols; i++) {
        if (colmark[i] == 0) {
            if (i != j) {
                beg[j] = beg[i];
                cnt[j] = cnt[i];
                qslp->obj[j] = qslp->obj[i];
                qslp->lower[j] = qslp->lower[i];
                qslp->upper[j] = qslp->upper[i];
            }
            newcolindex[i] = j++;
        } else {
            for (k = 0; k < cnt[i]; k++) {
                ind[beg[i] + k] = -1;
            }
            newcolindex[i] = -1;   
        }
    }

    /* Update the struct arrays */

    for (i = 0, j = 0; i < qslp->nstruct; i++) {
        k = qslp->structmap[i];
        if (colmark[k] == 0) {
            qslp->structmap[j] = newcolindex[k];
            qslp->colnames[j] = qslp->colnames[i];
            if (qslp->intmarker) qslp->intmarker[j] = qslp->intmarker[i];
            j++;
        } else {
	    rval = ILLsymboltab_delete(&qslp->coltab, qslp->colnames[i]); 
	    ILL_CLEANUP_IF(rval);
            ILL_IFFREE (qslp->colnames[i], char);
        }
    }

    /* Update the rowmap: note if logicals deleted, map will be -1 */

    for (i = 0; i < nrows; i++) {
        qslp->rowmap[i] = newcolindex[qslp->rowmap[i]];
    }

CLEANUP:

    ILL_IFFREE (newcolindex, int);
    ILL_RETURN (rval, "delcols_work");
}

int ILLlib_chgcoef (lpinfo *lp, int rowindex, int colindex, double coef)
{
    int rval = 0;
    ILLlpdata *qslp;
    ILLmatrix *A;
    int nrows, nstruct, j; 

    if (!lp) {
        fprintf (stderr, "ILLlib_chgcoef called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    A = &qslp->A;

    nrows = qslp->nrows;
    nstruct = qslp->nstruct;

    if (rowindex < 0 || rowindex >= nrows ||
        colindex < 0 || colindex >= nstruct) {
        fprintf (stderr, "ILLlib_chgcoef called with out-of-range index\n");
        rval = 1; ILL_CLEANUP;
    }

    if (qslp->rA) {  /* After a chgcoef call, needs to be updated */
        ILLlp_rows_clear (qslp->rA);
        ILL_IFFREE (qslp->rA, ILLlp_rows);
    }

    if (qslp->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (qslp->sinfo);
        ILL_IFFREE (qslp->sinfo, ILLlp_sinfo);
    }

    j = qslp->structmap[colindex];

    rval = matrix_addcoef (lp, A, rowindex, j, coef);
    ILL_CLEANUP_IF(rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgcoef");
}

int ILLlib_chgsense (lpinfo *lp, int num, int *rowlist, char *sense)
{
    int rval = 0;
    int i, j, k;
    ILLlpdata *qslp = lp->O;
    ILLmatrix *A = &(lp->O->A);

    for (i = 0; i < num; i++) {
        j = qslp->rowmap[rowlist[i]];
        qslp->sense[j] = sense[i];
        if (A->matcnt[j] != 1) {
            fprintf (stderr, "logical variable is not a singleton\n");
            rval = 1; ILL_CLEANUP;
        }
        k = A->matbeg[j];
        switch (sense[i]) {
        case 'E':                         /* Artificial */
            qslp->lower[j] = 0.0;
            qslp->upper[j] = 0.0;
            A->matval[k] = 1.0;
            break;
        case 'G':                         /* Surplus   */
            qslp->lower[j] = 0.0;
            qslp->upper[j] = ILL_MAXDOUBLE;
            A->matval[k] = -1.0;
            break;
        case 'L':                         /* Slack     */
            qslp->lower[j] = 0.0;
            qslp->upper[j] = ILL_MAXDOUBLE;
            A->matval[k] = 1.0;
            break;
        default:
            fprintf (stderr, "illegal sense %c in ILLlib_chgsense\n", sense[i]);
            rval = 1; ILL_CLEANUP;
        }
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgsense");
}

int ILLlib_newcol (lpinfo *lp, ILLlp_basis *B, double obj, double lower,
        double upper, const char *name, int factorok)
{
    int rval = 0;

    rval = ILLlib_addcol (lp, B, 0, (int *) NULL, (double *) NULL, obj,
                         lower, upper, name, factorok);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_newcol");
}

int ILLlib_newcols (lpinfo *lp, ILLlp_basis *B, int num, double *obj,
        double *lower, double *upper, const char **names, int factorok)
{
    int rval = 0;
    int *cmatcnt = (int *) NULL;
    int *cmatbeg = (int *) NULL;
    int i;

    ILL_SAFE_MALLOC (cmatcnt, num, int);

    ILL_SAFE_MALLOC (cmatbeg, num, int);
   
    for (i = 0; i < num; i++) {
        cmatcnt[i] = 0;
        cmatbeg[i] = 0;
    }

    rval = ILLlib_addcols (lp, B, num, cmatcnt, cmatbeg, (int *) NULL,
                         (double *) NULL, obj, lower, upper, names, factorok);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_IFFREE (cmatcnt, int);
    ILL_IFFREE (cmatbeg, int);

    ILL_RETURN (rval, "ILLlib_newcols");
}

int ILLlib_addcols (lpinfo *lp, ILLlp_basis *B, int num, int *cmatcnt,
       int *cmatbeg, int *cmatind, double *cmatval, double *obj, double *lower,
       double *upper, const char **names, int factorok)
{
    int rval = 0;
    int i;

    for (i = 0; i < num; i++) {
        if (names) {
            rval = ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
                        cmatval + cmatbeg[i], obj[i], lower[i], upper[i],
                        names[i], factorok);
        } else {
            rval = ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
                        cmatval + cmatbeg[i], obj[i], lower[i], upper[i],
                        (const char *) NULL, factorok);
        }
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_addcols");
}

int ILLlib_addcol (lpinfo *lp, ILLlp_basis *B, int cnt, int *ind, double *val,
        double obj, double lower, double upper, const char *name, int factorok)
{
    int rval = 0;
    ILLlpdata *qslp;
    ILLmatrix *A;
    int ncols;
    char buf[ILL_namebufsize];
    int pind, hit;
    double l, u;

    if (!lp) {
        fprintf (stderr, "ILLlib_addcol called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    A = &qslp->A;
    ncols = qslp->ncols;

    if (qslp->rA) {  /* After an addcol call, needs to be updated */
        ILLlp_rows_clear (qslp->rA);
        ILL_IFFREE (qslp->rA, ILLlp_rows);
    }

    if (qslp->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (qslp->sinfo);
        ILL_IFFREE (qslp->sinfo, ILLlp_sinfo);
    }


    /* Add the new variable to the column structures */

    if (qslp->colsize < ncols + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(qslp->lower),
                       qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->upper),
                       qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->obj),
                       qslp->colsize + EXTRA_COLS, sizeof (double));
        ILL_CLEANUP_IF (rval);

        qslp->colsize += EXTRA_COLS;
    }

    qslp->obj[ncols] = obj;
    qslp->lower[ncols] = lower;
    qslp->upper[ncols] = upper; 

    /*  Add the variable to the structural arrays */

    if (qslp->structsize < qslp->nstruct + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(qslp->structmap),
                       qslp->structsize + EXTRA_COLS, sizeof (int));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(qslp->colnames),
                       qslp->structsize + EXTRA_COLS, sizeof (char *));
        ILL_CLEANUP_IF (rval);

        if (qslp->intmarker) {
            rval = ILLutil_reallocrus_count ((void **) &(qslp->intmarker),
                           qslp->structsize + EXTRA_COLS, sizeof (char));
            ILL_CLEANUP_IF (rval);
        }
        qslp->structsize += EXTRA_COLS;
    }
    
    qslp->structmap[qslp->nstruct] = ncols;
    if (qslp->intmarker != (char *) NULL) {
        /* NOTE: If we want to add integer variables, this is the place. */
        qslp->intmarker[qslp->nstruct] = (char) 0;
    }

    ILL_FAILtrue(qslp->colnames == NULL, "must always be non NULL"); 
    ILLlib_findName(qslp, 0 /*isRow*/, name, qslp->nstruct, buf); 
    ILLsymboltab_register(&qslp->coltab, buf, qslp->nstruct, &pind, &hit); 
    ILL_FAILfalse((pind == qslp->nstruct) && (hit == 0), "must be new"); 
    ILL_UTIL_STR(qslp->colnames[qslp->nstruct], buf);


    /*  Add col to the matrix */

    rval = matrix_addcol (A, cnt, ind, val);
    ILL_CLEANUP_IF (rval);


    if (B != (ILLlp_basis *) NULL) {
        rval = ILLutil_reallocrus_count ((void **) &(B->cstat),
                                        qslp->nstruct + 1, sizeof (char));
        ILL_CLEANUP_IF (rval);
        if (lower == -ILL_MAXDOUBLE && upper == ILL_MAXDOUBLE) {
            B->cstat[qslp->nstruct] = QS_COL_BSTAT_FREE;
        } else if (upper == ILL_MAXDOUBLE) {
            B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
        } else if (lower == ILL_MAXDOUBLE) {
            B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
        } else {
            l = fabs(lower); 
            u = fabs(upper);
            if (l < u) {
                B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
            } else {
                B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
            }
        }

        /* UPDATE THE PINF PRIMAL NORMS */

        ILL_IFFREE (B->colnorms, double);
    }

    if (factorok == 0) {
#if 0
        lp->basisid = -1;  /* To get optimizer to reload the basis */
#endif
    } else {
        if (!lp->nbaz || !lp->vindex || !lp->vstat) {
            fprintf (stderr, "ERROR: factorok set without a current basis\n");
            rval = 1; ILL_CLEANUP;
        }

        rval = ILLutil_reallocrus_count ((void **) &(lp->nbaz),
                       qslp->nstruct + 1, sizeof (int)); 
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(lp->vindex), 
                       qslp->ncols + 1, sizeof (int));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(lp->vstat),
                       qslp->ncols + 1, sizeof (int));


        lp->nbaz[qslp->nstruct] = qslp->ncols;
        lp->vindex[qslp->ncols] = qslp->nstruct;

        if (lower == -ILL_MAXDOUBLE && upper == ILL_MAXDOUBLE) {
            lp->vstat[qslp->ncols] = STAT_ZERO;
        } else if (upper == ILL_MAXDOUBLE) {
            lp->vstat[qslp->ncols] = STAT_LOWER;
        } else if (lower == ILL_MAXDOUBLE) {
            lp->vstat[qslp->ncols] = STAT_UPPER;
        } else {
            l = fabs(lower); 
            u = fabs(upper);
            if (l < u) {
                lp->vstat[qslp->ncols] = STAT_LOWER;
            } else {
                lp->vstat[qslp->ncols] = STAT_UPPER;
            }
        }
    }


    qslp->ncols++;
    qslp->nstruct++;
    (qslp->nzcount) += cnt;

    if (B != (ILLlp_basis *) NULL) {
        B->nstruct++;
    }

CLEANUP:
    ILL_RETURN (rval, "ILLlib_addcol");
}

static int matrix_addrow (ILLmatrix *A, int rowcnt, int *rowind, double *rowval)
{
    int rval = 0;
    int i, j, k, ind, memo, stop, delta = 0;

    /* matsize will be the length of the array.                   */
    /* matfree will keep track of the free space at end of array. */

    for (i = 0; i < rowcnt; i++) {
        if (rowind[i] >= A->matcols || rowind[i] < 0) {
            fprintf (stderr, "illegal col index in matrix_addrow\n");
            rval= 1; ILL_CLEANUP;
        }
    }

    for (i = 0; i < rowcnt; i++) {
        j = rowind[i];
        if (A->matcnt[j] > 0 && 
           (A->matbeg[j] + A->matcnt[j] + 1 > A->matsize ||
            A->matind[A->matbeg[j] + A->matcnt[j]] != -1)) {
            delta += (A->matcnt[j] + 2);  /* 1 for the new coef and 1 for */
                                          /* an extra space               */
        }
    }

    if (delta < A->matfree) {  
        for (i = 0; i < rowcnt; i++) {
            j = rowind[i];
            if (A->matcnt[j] == 0) {
                A->matind[A->matbeg[j]] = A->matrows;
                A->matval[A->matbeg[j]] = rowval[i];
                A->matcnt[j] = 1;
            } else if (A->matind[A->matbeg[j] + A->matcnt[j]] == -1) {
                /* Since A->matfree is positive, we know that we are not */
                /* sitting at the end of the array.                      */
                A->matind[A->matbeg[j] + A->matcnt[j]] = A->matrows;
                A->matval[A->matbeg[j] + A->matcnt[j]] = rowval[i];
                if ((A->matbeg[j]+A->matcnt[j]) == (A->matsize-A->matfree)) {
                    A->matfree--;   /* at end of used space */
                }
                (A->matcnt[j])++;
            } else {
                ind = A->matsize - A->matfree + 1;  /* leave space for -1 */
                memo = ind;
                stop = A->matbeg[j] + A->matcnt[j];
                for (k = A->matbeg[j]; k < stop; k++) {
                    if (ind >= A->matsize) {
                        printf ("WHAT: %d, %d\n", A->matsize, ind);
                        fflush (stdout);
                        exit (1);
                    }
                    A->matind[ind] = A->matind[k];
                    A->matval[ind] = A->matval[k];
                    A->matind[k] = -1;
                    ind++;
                }
                A->matind[ind] = A->matrows;
                A->matval[ind] = rowval[i];
                A->matbeg[j] = memo;
                (A->matcnt[j])++;
                (A->matfree) -= (A->matcnt[j] + 1);
            }
        }
    } else {
        rval = matrix_addrow_end (A, A->matrows, rowcnt, rowind, rowval);
        ILL_CLEANUP_IF (rval);
    }
    A->matrows++;

CLEANUP:

    ILL_RETURN (rval, "matrix_addrow");
}

static int matrix_addrow_end (ILLmatrix *A, int row, int rowcnt, int *rowind,
        double *rowval)
{
    int rval = 0;
    int i, j, k, start, stop, total;
    int *newbeg = (int *) NULL;
    int *newind = (int *) NULL;
    double *newval = (double *) NULL;
    int ncols = A->matcols;
    
    if (A->matcolsize > 0) {
        ILL_SAFE_MALLOC (newbeg, A->matcolsize, int);
    }
    ILL_SAFE_MALLOC (newind, A->matsize + rowcnt + EXTRA_MAT, int);
    ILL_SAFE_MALLOC (newval, A->matsize + rowcnt + EXTRA_MAT, double);

    A->matsize += (rowcnt + EXTRA_MAT);

    for (i = 0; i < rowcnt; i++) {
        A->matcnt[rowind[i]]++;
    }
    for (total = 0, j = 0; j < ncols; j++) {
        newbeg[j] = total;
        if (A->matcnt[j] > 0) total += A->matcnt[j];
        else                  total += 1;
    }
    for (i = 0; i < rowcnt; i++) {
        A->matcnt[rowind[i]]--;
    }
    for (j = total; j < A->matsize; j++) {
        newind[j] = -1;
    }
    A->matfree = A->matsize - total;

    for (j = 0; j < ncols; j++) {
        if (A->matcnt[j] > 0) {
            stop = A->matbeg[j] + A->matcnt[j];
            start = newbeg[j];
            for (k = A->matbeg[j]; k < stop; k++) {
                newind[start] = A->matind[k];
                newval[start] = A->matval[k];
                start++;
            }
        } else {
            newind[newbeg[j]] = 1;
        }
    }
    for (i = 0; i < rowcnt; i++) {
        j = rowind[i];
        newind[newbeg[j] + A->matcnt[j]] = row;
        newval[newbeg[j] + A->matcnt[j]] = rowval[i];
        (A->matcnt[j])++;
    }

    ILL_IFFREE (A->matbeg, int);
    ILL_IFFREE (A->matind, int);
    ILL_IFFREE (A->matval, double);

    A->matbeg = newbeg;
    A->matind = newind;
    A->matval = newval;

CLEANUP:

    if (rval) {
        ILL_IFFREE (newbeg, int);
        ILL_IFFREE (newind, int);
        ILL_IFFREE (newval, double);
    }

    ILL_RETURN (rval, "matrix_addrow_end");
}

static int matrix_addcoef (lpinfo *lp, ILLmatrix *A, int row, int col,
        double val)
{
    int i, k, delta, ind, stop, memo;
    int tind[1];
    double tval[1];
    int rval = 0;

    if (row >= A->matrows || row < 0) {
        fprintf (stderr, "illegal row index in matrix_addcoef\n");
        rval= 1; ILL_CLEANUP;
    }

    if (col >= A->matcols || col < 0) {
        fprintf (stderr, "illegal col index in matrix_addcoef\n");
        rval= 1; ILL_CLEANUP;
    }

    for (i = A->matbeg[col]; i < A->matbeg[col] + A->matcnt[col]; i++) {
        if (A->matind[i] == row) {
            A->matval[i] = val;
            ILL_CLEANUP;
        }
    }

    /* The coef is new, we need to add it to A */

    lp->O->nzcount++;
    delta = A->matcnt[col] + 2;

    if (A->matcnt[col] == 0) {
        /* First entry, always a free space */
        A->matind[A->matbeg[col]] = row;
        A->matval[A->matbeg[col]] = val;
        A->matcnt[col] = 1;
    } else if (A->matbeg[col] + A->matcnt[col] < A->matsize &&
               A->matind[A->matbeg[col] + A->matcnt[col]] == -1) {
        /* Free space in the column */
        A->matind[A->matbeg[col] + A->matcnt[col]] = row;
        A->matval[A->matbeg[col] + A->matcnt[col]] = val;
        if ((A->matbeg[col] + A->matcnt[col]) == (A->matsize - A->matfree)) {
            A->matfree--;
        }
        (A->matcnt[col])++;
    } else if (A->matfree > delta) {
        /* Enough space to move column to end of array */
        ind = A->matsize - A->matfree + 1;
        memo = ind;
        stop = A->matbeg[col] + A->matcnt[col];
        for (k = A->matbeg[col]; k < stop; k++) {
            A->matind[ind] = A->matind[k];
            A->matval[ind] = A->matval[k];
            A->matind[k] = -1;
            ind++;
        }
        A->matind[ind] = row;
        A->matval[ind] = val;

        A->matbeg[col] = memo;
        (A->matcnt[col])++;
        (A->matfree) -= (A->matcnt[col] + 1);
    } else {
        /* Need to malloc space to move column to end of array */

        tind[0] = col;
        tval[0] = val;

        rval = matrix_addrow_end (A, row, 1, tind, tval);
        ILL_CLEANUP_IF(rval); 
    }

CLEANUP:

    ILL_RETURN (rval, "matrix_addcoef");
}

static int matrix_addcol (ILLmatrix *A, int colcnt, int *colind, double *colval)
{
    int rval = 0;
    int i, ind;

    for (i = 0; i < colcnt; i++) {
        if (colind[i] >= A->matrows || colind[i] < 0) {
            fprintf (stderr, "illegal row index in matrix_addcol\n");
            rval= 1; ILL_CLEANUP;
        }
    }

    if (A->matcolsize < A->matcols + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(A->matbeg),
                       A->matcolsize + EXTRA_COLS, sizeof (int));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(A->matcnt),
                       A->matcolsize + EXTRA_COLS, sizeof (int));
        ILL_CLEANUP_IF (rval);

        (A->matcolsize) += EXTRA_COLS;
    }

    if (A->matfree < colcnt + 1) {
        rval = ILLutil_reallocrus_count ((void **) &(A->matind),
                       A->matsize + colcnt + EXTRA_MAT + 1, sizeof (int));
        ILL_CLEANUP_IF (rval);

        rval = ILLutil_reallocrus_count ((void **) &(A->matval),
                       A->matsize + colcnt + EXTRA_MAT + 1, sizeof (double));
        ILL_CLEANUP_IF (rval);

        for (i = 0; i < colcnt + EXTRA_MAT + 1; i++) {
            A->matind[A->matsize + i] = -1;
        }
        A->matsize += (colcnt + EXTRA_MAT + 1);
        A->matfree += (colcnt + EXTRA_MAT + 1);
    }

    ind = A->matsize - A->matfree;
    A->matbeg[A->matcols] = ind;
    A->matcnt[A->matcols] = colcnt;
    if (colcnt == 0) {
       A->matind[ind] = 1;   /* Dummy value to stop columns from stealing */
                             /* this space in addrows.                    */ 
       A->matfree -= 1;
    } else {
        for (i = 0; i < colcnt; i++) {
            A->matval[ind] = colval[i];
            A->matind[ind] = colind[i];
            ind++;
        }
        A->matfree -= colcnt;
    }
    A->matcols++;

CLEANUP:

    ILL_RETURN (rval, "matrix_addcol");
}

int ILLlib_getrows (lpinfo *lp, int num, int *rowlist, int **rowcnt, 
        int **rowbeg, int **rowind, double **rowval, double **rhs,
        char **sense, char ***names)
{
    int rval = 0;
    int *allbeg = (int *) NULL;
    int *allcnt = (int *) NULL;
    int *allind = (int *) NULL;
    double *allval = (double *) NULL;
    int i, row, k, start, stop, len, tcnt, cnt = 0;
    ILLlpdata *qslp;
    ILLlp_rows lprows; 

    if (rowcnt) *rowcnt = (int *) NULL;
    if (rowbeg) *rowbeg = (int *) NULL;
    if (rowind) *rowind = (int *) NULL;
    if (rowval) *rowval = (double *) NULL;
    if (rhs)    *rhs    = (double *) NULL;
    if (sense)  *sense  = (char *) NULL;
    if (names)  *names  = (char **) NULL;

    if (!lp) {
        fprintf (stderr, "ILLlib_getrows called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    if (!num) ILL_CLEANUP;

    qslp = lp->O;

    rval = ILLlp_rows_init(&lprows, qslp, 0); 
    ILL_CLEANUP_IF(rval);
    allbeg = lprows.rowbeg;
    allcnt = lprows.rowcnt;
    allind = lprows.rowind;
    allval = lprows.rowval;

    for (i = 0; i < num; i++) {
        cnt += allcnt[rowlist[i]];
    }

    if (rowcnt) {
        ILL_SAFE_MALLOC (*rowcnt, num, int);
        for (i = 0; i < num; i++) {
            (*rowcnt)[i] = allcnt[rowlist[i]];  
        }
    }

    if (rowbeg) {
        ILL_SAFE_MALLOC (*rowbeg, num, int);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            (*rowbeg)[i] = tcnt;
            tcnt += allcnt[rowlist[i]];
        }
    }

    if (cnt && rowind) {
        ILL_SAFE_MALLOC (*rowind, cnt, int);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            row = rowlist[i];
            start = allbeg[row];
            stop = start + allcnt[row];
            for (k = start; k < stop; k++) {
                (*rowind)[tcnt++] = allind[k];
            }
        }
    }

    if (cnt && rowval) {
        ILL_SAFE_MALLOC (*rowval, cnt, double);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            row = rowlist[i];
            start = allbeg[row];
            stop = start + allcnt[row];
            for (k = start; k < stop; k++) {
                (*rowval)[tcnt++] = allval[k];
            }
        }
    }

    if (rhs) {
        ILL_SAFE_MALLOC (*rhs, num, double);
        for (i = 0; i < num; i++) {
            (*rhs)[i] = qslp->rhs[rowlist[i]];
        }
    }

    if (sense) {
        ILL_SAFE_MALLOC (*sense, num, char);
        for (i = 0; i < num; i++) {
            (*sense)[i] = qslp->sense[rowlist[i]];
        }
    }

    if (names) {
        if (qslp->rownames == (char **) NULL) {
            fprintf (stderr, "LP does not have row names\n");
            rval = 1; ILL_CLEANUP;
        }
        ILL_SAFE_MALLOC (*names, num, char *);
        for (i = 0; i < num; i++) {
            (*names)[i] = (char *) NULL;
        }
        for (i = 0; i < num; i++) {
            len = strlen (qslp->rownames[rowlist[i]]) + 1;
            ILL_SAFE_MALLOC ((*names)[i], len, char);
            strcpy ((*names)[i], qslp->rownames[rowlist[i]]);
        }
    }

CLEANUP:

    ILL_IFFREE (allbeg, int);
    ILL_IFFREE (allcnt, int);
    ILL_IFFREE (allind, int);
    ILL_IFFREE (allval, double);

    if (rval) {
        if (rowcnt) ILL_IFFREE (*rowcnt, int);
        if (rowbeg) ILL_IFFREE (*rowbeg, int);
        if (rowind) ILL_IFFREE (*rowind, int);
        if (rowval) ILL_IFFREE (*rowval, double);
        if (rhs)    ILL_IFFREE (*rhs, double);
        if (sense)  ILL_IFFREE (*sense, char);
        if (names && (*names)) {
            for (i = 0; i < num; i++) {
                ILL_IFFREE ((*names)[i], char);
            } 
            ILL_IFFREE (*names, char *);
        }
    }

    ILL_RETURN (rval, "ILLlib_getrows");
}

int ILLlib_getcols (lpinfo *lp, int num, int *collist, int **colcnt, 
        int **colbeg, int **colind, double **colval, double **obj,
        double **lower, double **upper, char ***names)
{
    int rval = 0;
    int i, col, k, start, stop, len, tcnt, cnt = 0;
    ILLlpdata *qslp;
    ILLmatrix *A;
    int *tlist = (int *) NULL;

    if (colcnt) *colcnt = (int *) NULL;
    if (colbeg) *colbeg = (int *) NULL;
    if (colind) *colind = (int *) NULL;
    if (colval) *colval = (double *) NULL;
    if (lower)  *lower  = (double *) NULL;
    if (upper)  *upper  = (double *) NULL;
    if (obj)    *obj    = (double *) NULL;
    if (names)  *names  = (char **) NULL;

    if (!lp) {
        fprintf (stderr, "ILLlib_getcols called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    if (!num) ILL_CLEANUP;

    qslp = lp->O;
    A = &(qslp->A);

    ILL_SAFE_MALLOC (tlist, num, int);
    for (i = 0; i < num; i++) {
        tlist[i] = qslp->structmap[collist[i]];
    }

    for (i = 0; i < num; i++) {
        cnt += A->matcnt[tlist[i]];
    }

    if (colcnt) {
        ILL_SAFE_MALLOC (*colcnt, num, int);
        for (i = 0; i < num; i++) {
            (*colcnt)[i] = A->matcnt[tlist[i]];  
        }
    }

    if (colbeg) {
        ILL_SAFE_MALLOC (*colbeg, num, int);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            (*colbeg)[i] = tcnt;
            tcnt += A->matcnt[tlist[i]];
        }
    }

    if (cnt && colind) {
        ILL_SAFE_MALLOC (*colind, cnt, int);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            col = tlist[i];
            start = A->matbeg[col];
            stop = start + A->matcnt[col];
            for (k = start; k < stop; k++) {
                (*colind)[tcnt++] = A->matind[k];
            }
        }
    }

    if (cnt && colval) {
        ILL_SAFE_MALLOC (*colval, cnt, double);
        tcnt = 0;
        for (i = 0; i < num; i++) {
            col = tlist[i];
            start = A->matbeg[col];
            stop = start + A->matcnt[col];
            for (k = start; k < stop; k++) {
                (*colval)[tcnt++] = A->matval[k];
            }
        }
    }

    if (obj) {
        ILL_SAFE_MALLOC (*obj, num, double);
        for (i = 0; i < num; i++) {
            (*obj)[i] = qslp->obj[tlist[i]];
        }
    }

    if (lower) {
        ILL_SAFE_MALLOC (*lower, num, double);
        for (i = 0; i < num; i++) {
            (*lower)[i] = qslp->lower[tlist[i]];
        }
    }

    if (upper) {
        ILL_SAFE_MALLOC (*upper, num, double);
        for (i = 0; i < num; i++) {
            (*upper)[i] = qslp->upper[tlist[i]];
        }
    }

    if (names) {
        if (qslp->colnames == (char **) NULL) {
            fprintf (stderr, "LP does not have col names\n");
            rval = 1; ILL_CLEANUP;
        }
        ILL_SAFE_MALLOC (*names, num, char *);
        for (i = 0; i < num; i++) {
            (*names)[i] = (char *) NULL;
        }
        for (i = 0; i < num; i++) {
            len = strlen (qslp->colnames[collist[i]]) + 1;
            ILL_SAFE_MALLOC ((*names)[i], len, char);
            strcpy ((*names)[i], qslp->colnames[collist[i]]);
        }
    }

CLEANUP:

    if (rval) {
        if (colcnt) ILL_IFFREE (*colcnt, int);
        if (colbeg) ILL_IFFREE (*colbeg, int);
        if (colind) ILL_IFFREE (*colind, int);
        if (colval) ILL_IFFREE (*colval, double);
        if (obj)    ILL_IFFREE (*obj, double);
        if (lower)  ILL_IFFREE (*lower, double);
        if (upper)  ILL_IFFREE (*upper, double);
        if (names && (*names)) {
            for (i = 0; i < num; i++) {
                ILL_IFFREE ((*names)[i], char);
            } 
            ILL_IFFREE (*names, char *);
        }
    }
    ILL_IFFREE (tlist, int);

    ILL_RETURN (rval, "ILLlib_getcols");
}

int ILLlib_getobj (lpinfo *lp, double *obj)
{
    ILLlpdata *qslp;
    int nstruct, j;
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "ILLlib_getobj called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nstruct = qslp->nstruct;

    for (j = 0; j < nstruct; j++) {
        obj[j] = qslp->obj[j];
    } 

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getobj");
}

int ILLlib_chgobj (lpinfo *lp, int indx, double coef)
{
    int rval = 0;
    int col;

    if (!lp) {
        fprintf (stderr, "ILLlib_chgobj called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (indx < 0 || indx >= lp->O->nstruct) {
        fprintf (stderr, "ILLlib_chgrhs called with bad indx: %d\n", indx);
        rval = 1; ILL_CLEANUP;
    }

    if (lp->O->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (lp->O->sinfo);
        ILL_IFFREE (lp->O->sinfo, ILLlp_sinfo);
    }

    col = lp->O->structmap[indx];
    lp->O->obj[col] = coef;

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgobj");
}

int ILLlib_getrhs (lpinfo *lp, double *rhs)
{
    ILLlpdata *qslp;
    int nrows, i;
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "ILLlib_getrhs called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nrows = qslp->nrows;

    for (i = 0; i < nrows; i++) {
        rhs[i] = qslp->rhs[i];
    } 

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getrhs");
}

int ILLlib_chgrhs (lpinfo *lp, int indx, double coef)
{
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "ILLlib_chgrhs called without an lp\n");
        rval = 1; ILL_CLEANUP;
    }

    if (indx < 0 || indx >= lp->O->nrows) {
        fprintf (stderr, "ILLlib_chgrhs called with bad indx: %d\n", indx);
        rval = 1; ILL_CLEANUP;
    }

    if (lp->O->sinfo) {   /* Presolve LP is no longer valid, free the data */
        ILLlp_sinfo_free (lp->O->sinfo);
        ILL_IFFREE (lp->O->sinfo, ILLlp_sinfo);
    }

    lp->O->rhs[indx] = coef;

CLEANUP:

    ILL_RETURN (rval, "ILLlib_chgrhs");
}

int ILLlib_rownames (lpinfo *lp, char **rownames)
{
    ILLlpdata *qslp;
    int nrows, len, i, rcount = 0;
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "ILLlib_rownames called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }
    if (!rownames) {
        fprintf (stderr, "ILLlib_rownames called with NULL rownames\n"); 
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nrows = qslp->nrows;

    if (qslp->rownames == (char **) NULL) {
        fprintf (stderr, "LP does not have rownames assigned\n");
        rval = 1; ILL_CLEANUP;
    }

    for (i = 0; i < nrows; i++) {
        len = strlen (qslp->rownames[i]) + 1;
        ILL_SAFE_MALLOC (rownames[i], len, char);
        strcpy (rownames[i], qslp->rownames[i]);
        rcount++;
    }

CLEANUP:

    if (rval) {
        for (i = 0; i < rcount; i++) {
            ILL_IFFREE (rownames[i], char);
        }
    }
    ILL_RETURN (rval, "ILLlib_rownames");
}

int ILLlib_getintflags (lpinfo *lp, int *intflags)
{
    int j, nstruct, rval = 0;
    ILLlpdata *qslp;

    if (!lp) {
        fprintf (stderr, "ILLlib_getintflags called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nstruct = qslp->nstruct;

    if (qslp->intmarker == (char *) NULL) {
        for (j = 0; j < nstruct; j++) {
            intflags[j] = 0;
        }
    } else {
        for (j = 0; j < nstruct; j++) {
            if (qslp->intmarker[j]) {
                intflags[j] = 1;
            } else {
                intflags[j] = 0;
            }
        }
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getintflags");
}

int ILLlib_colnames (lpinfo *lp, char **colnames)
{
    ILLlpdata *qslp;
    int nstruct, len, i, ccount = 0;
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "ILLlib_colnames called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }
    if (!colnames) {
        fprintf (stderr, "ILLlib_colnames called with NULL colnames\n"); 
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nstruct = qslp->nstruct;

    if (qslp->colnames == (char **) NULL) {
        fprintf (stderr, "LP does not have colnames assigned\n");
        rval = 1; ILL_CLEANUP;
    }

    for (i = 0; i < nstruct; i++) {
        len = strlen (qslp->colnames[i]) + 1;
        ILL_SAFE_MALLOC (colnames[i], len, char);
        strcpy (colnames[i], qslp->colnames[i]);
        ccount++;
    }


CLEANUP:

    if (rval) {
        for (i = 0; i < ccount; i++) {
            ILL_IFFREE (colnames[i], char);
        }
    }

    ILL_RETURN (rval, "ILLlib_colnames");
}

static int reset_colindex (lpinfo *lp)
{
    int rval = 0;
    int test;
    ILLlpdata *qslp = lp->O;
   
    test = ILLsymboltab_index_ok (&qslp->coltab);
    if (!test) {
        rval = ILLsymboltab_index_reset (&qslp->coltab, qslp->nstruct,
                                         qslp->colnames);
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_RETURN (rval, "reset_colindex");
}

static int reset_rowindex (lpinfo *lp)
{
    int rval = 0;
    int test;
    ILLlpdata *qslp = lp->O;
   
    test = ILLsymboltab_index_ok (&qslp->rowtab);
    if (!test) {
        rval = ILLsymboltab_index_reset (&qslp->rowtab, qslp->nrows,
                                         qslp->rownames);
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_RETURN (rval, "reset_rowindex");
}

int ILLlib_colindex (lpinfo *lp, const char *name, int *colindex)
{
    int rval = 0;
    ILLlpdata *qslp;
   
    *colindex = -1;

    if (!lp) {
        fprintf (stderr, "ILLlib_colindex called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;

    rval = reset_colindex (lp);
    ILL_CLEANUP_IF (rval);

    rval = ILLsymboltab_getindex (&qslp->coltab, name, colindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_colindex");
}

int ILLlib_rowindex (lpinfo *lp, const char *name, int *rowindex)
{
    int rval = 0;
    ILLlpdata *qslp;

    *rowindex = -1;

    if (!lp) {
        fprintf (stderr, "ILLlib_rowindex called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;

    rval = reset_rowindex (lp);
    ILL_CLEANUP_IF (rval);

    rval = ILLsymboltab_getindex (&qslp->rowtab, name, rowindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_rowindex");
}

int ILLlib_getbasis (lpinfo *lp, char *cstat, char *rstat)
{
    int rval = 0;
    int i, j, nrows;
    ILLlpdata *qslp;

    if (!lp) {
        fprintf (stderr, "ILLlib_getbasis called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }

    if (lp->basisid == -1) {
        fprintf (stderr, "ILLlib_getbasis called with modifed LP\n");
        rval = 1; ILL_CLEANUP;
    }

    nrows = lp->nrows;
    qslp  = lp->O;

    for (i = 0; i < qslp->nstruct; i++) {
        j = qslp->structmap[i];
        switch (lp->vstat[j]) {
        case STAT_BASIC:
            cstat[i] = QS_COL_BSTAT_BASIC;
            break;
        case STAT_LOWER:
            cstat[i] = QS_COL_BSTAT_LOWER;
            break;
        case STAT_UPPER:
            cstat[i] = QS_COL_BSTAT_UPPER;
            break;
        case STAT_ZERO:
            cstat[i] = QS_COL_BSTAT_FREE;
            break;
        default:
            fprintf (stderr, "unknown vstat in ILLlib_getbasis: %d\n",
                     lp->vstat[j]);
            rval = 1; ILL_CLEANUP;
        }
    }

    for (i = 0; i < nrows; i++) {
        j = qslp->rowmap[i];
        if (qslp->rangeval && qslp->rangeval[i] != 0.0) {
            switch (lp->vstat[j]) {
            case STAT_BASIC:
                rstat[i] = QS_ROW_BSTAT_BASIC;
                break;
            case STAT_LOWER:
                rstat[i] = QS_ROW_BSTAT_LOWER;
                break;
            case STAT_UPPER:
                rstat[i] = QS_ROW_BSTAT_UPPER;
                break;
            default:
                fprintf (stderr, "unknown vstat in ILLlib_getbasis 2\n");
                rval = 1; ILL_CLEANUP;
            }
        } else {
            switch (lp->vstat[j]) {
            case STAT_BASIC:
                rstat[i] = QS_ROW_BSTAT_BASIC;
                break;
            case STAT_UPPER:
            case STAT_LOWER:
                rstat[i] = QS_ROW_BSTAT_LOWER;
                break;
            default:
                fprintf (stderr, "unknown vstat in ILLlib_getbasis 3: %d, %d\n",
                                  i, lp->vstat[j]);
                rval = 1; ILL_CLEANUP;
            }
        }
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_getbasis");
}

int ILLlib_loadbasis (ILLlp_basis *B, int nstruct, int nrows, char *cstat,
        char *rstat)
{
    int i;
    int rval = 0;

    ILLlp_basis_init (B);

    if (!cstat || !rstat) {
        rval = 1;
        ILL_CLEANUP_IF (rval);
    }

    rval = ILLlp_basis_alloc (B, nstruct, nrows);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < nstruct; i++) {
        B->cstat[i] = cstat[i];
    }
    for (i = 0; i < nrows; i++) {
        B->rstat[i] = rstat[i];
    }

CLEANUP:

    ILL_RETURN (rval, "ILLlib_loadbasis");
}

#define READ_BASIS_XL 0
#define READ_BASIS_XU 1
#define READ_BASIS_LL 2
#define READ_BASIS_UL 3

int ILLlib_readbasis (lpinfo *lp, ILLlp_basis *B, const char *fname)
{
    int rval = 0;
    ILLlpdata *qslp = lp->O;
    int nstruct = qslp->nstruct;
    int nrows = qslp->nrows;
    int i, j, end = 0, sec, havename = 0;
    int rowtype, row, col;
    char *bname = (char *) NULL;
    FILE *file_in = (FILE *) NULL;
    ILLread_mps_state state;
	qsline_reader *in = NULL; 

    ILLlp_basis_init (B);

    ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
    ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);
    B->nstruct = nstruct;
    B->nrows   = nrows;

    for (j = 0; j < nstruct; j++) {
        B->cstat[j] = QS_COL_BSTAT_LOWER;
    }
    for (i = 0; i < nrows; i++) {
        B->rstat[i] = QS_ROW_BSTAT_BASIC;
    }

    file_in = fopen (fname, "r");
    if (file_in == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for reading\n", fname);
        rval = 1; ILL_CLEANUP;
    }

    in = ILLline_reader_new((qsread_line_fct) fgets, file_in); 
	rval = ILLmps_state_init (&state, in, fname);
    ILL_CLEANUP_IF (rval);

    while (ILLmps_next_line (&state) == 0) {
        if (ILLmps_empty_key(&state)) {

            /* Get the XL XU LL UL line */

            if (!havename) {
                rval = ILLmps_error (&state, "BASIS data before NAME\n");
                ILL_CLEANUP;
            }

            if (!strcmp(state.field, "XL")) {
                rowtype = READ_BASIS_XL;
            } else if (!strcmp(state.field, "XU")) {
                rowtype = READ_BASIS_XU;
            } else if (!strcmp(state.field, "LL")) {
                rowtype = READ_BASIS_LL;
            } else if (!strcmp(state.field, "UL")) {
                rowtype = READ_BASIS_UL;
            } else {
                rval = ILLmps_error(&state, "BASIS \"%s\" is invalid\n",
                                    state.field);
                ILL_CLEANUP;
            }

            if (ILLmps_next_field(&state) == 0) {

                rval = ILLlib_colindex (lp, (const char *) state.field, &col);
                ILL_CLEANUP_IF (rval);
                if (col == -1) {
                    rval = ILLmps_error(&state, "BASIS col not in LP\n");
                    ILL_CLEANUP;
                }

                if (rowtype == READ_BASIS_XL || rowtype == READ_BASIS_XU) {
                    if (ILLmps_next_field(&state) == 0) {
                        rval = ILLlib_rowindex (lp, (const char *) state.field,
                                                &row);
                        ILL_CLEANUP_IF (rval);
                        if (row == -1) {
                            rval = ILLmps_error(&state,
                                         "BASIS row not in LP\n");
                            ILL_CLEANUP;
                        }
                        if (rowtype == READ_BASIS_XL) {
                            B->cstat[col] = QS_COL_BSTAT_BASIC;
                            B->rstat[row] = QS_ROW_BSTAT_LOWER;

                        } else {
                            B->cstat[col] = QS_COL_BSTAT_BASIC;
                            B->rstat[row] = QS_ROW_BSTAT_UPPER;
                        }
                    } else {
                        rval = ILLmps_error(&state,
                                     "BASIS line needs row and column\n");
                        ILL_CLEANUP;
                    }
                } else {
                    if (rowtype == READ_BASIS_LL) {
                        B->cstat[col] = QS_COL_BSTAT_LOWER;
                    } else {
                        B->cstat[col] = QS_COL_BSTAT_UPPER;
                    }
                }
            } else {
                rval = ILLmps_error(&state, "BASIS line has no row/column\n");
                ILL_CLEANUP;
            }
        } else {
            /* found a section indicator in col 1 */
            if (!strcmp(state.key, ILLmps_section_name[ILL_MPS_ENDATA])) {
                end = 1;
                break;  /* done reading */
            }

            sec = ILLutil_index(ILLmps_section_name, state.key);
            if (sec < 0 || sec != ILL_MPS_NAME) {
                rval = ILLmps_error(&state, "BASIS \"%s\" is not a key\n",
                                    state.key);
                ILL_CLEANUP;
            }

            if (havename) {
                rval = ILLmps_error (&state, "BASIS two name sections\n");
                ILL_CLEANUP;
            }

            havename = 1;

            if (ILLmps_empty_field(&state)) {
                ILLmps_warn(&state, "BASIS blank NAME.");
            } else {
                ILL_UTIL_STR(bname, state.field);
                printf ("Basis Name: %s\n", bname);
                fflush (stdout);
                if (strcmp (bname, qslp->probname)) {
                    ILLmps_warn(&state, "BASIS name does not match LP.");
                }
            } 
        }
    }

    if (!end) {
        ILLmps_warn(&state, "Missing ENDATA in basis file.");
    }
    if (!ILLmps_next_line(&state)) {
        ILLmps_warn(&state, "Ignoring text after ENDATA.");
    }

    if (!havename) {
        rval = ILLmps_error (&state, "BASIS no name section\n");
        ILL_CLEANUP;
    }

    /* Correct the free variables */

    for (j = 0; j < nstruct; j++) {
        col = lp->O->structmap[j];
        if (qslp->lower[col] == -ILL_MAXDOUBLE &&
            qslp->upper[col] == ILL_MAXDOUBLE &&
                 B->cstat[j] == QS_COL_BSTAT_LOWER) {
            B->cstat[j] = QS_COL_BSTAT_FREE;
        }
    }

CLEANUP:

    if (file_in) fclose (file_in);
	ILLline_reader_free(in); 

    if (rval) {
        ILLlp_basis_free (B);
    }
    ILL_IFFREE (bname, char);

    ILL_RETURN (rval, "ILLlib_readbasis");
}

int ILLlib_writebasis (lpinfo *lp, ILLlp_basis *B, const char *fname)
{
    int rval = 0;
    FILE *out = (FILE *) NULL;
    char *cstat = (char *) NULL;
    char *rstat = (char *) NULL;
    ILLlpdata *qslp;
    int i, j, nstruct, nrows; 

    /* NOTE: non-basic free variables are encoded as non-basic at lower */

    if (!lp) {
        fprintf (stderr, "ILLlib_writebasis called without an LP\n");
        rval = 1; ILL_CLEANUP;
    }
    if (!B && lp->basisid == -1) {
        fprintf (stderr, "ILLlib_writebasis called with unsolved LP\n");
        rval = 1; ILL_CLEANUP;
    }

    qslp  = lp->O;
    nstruct = qslp->nstruct;
    nrows   = qslp->nrows;

    out = fopen (fname, "w");
    if (out == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for writing\n", fname);
        rval = 1; ILL_CLEANUP;
    }

    if (B) {
        cstat = B->cstat;
        rstat = B->rstat;
    } else {
        ILL_SAFE_MALLOC (cstat, nstruct, char);
        ILL_SAFE_MALLOC (rstat, nrows, char);

        rval = ILLlib_getbasis (lp, cstat, rstat);
        ILL_CLEANUP_IF (rval);
    }

    fprintf (out, "NAME    %s\n", qslp->probname);

    /* Pick out the non-basic rows and find a matching basic column */

    i = 0;
    j = 0;
    do {
        while (i < nrows && rstat[i] == QS_ROW_BSTAT_BASIC) {
            i++;
        }
        if (i < nrows) {
            while (j < nstruct && cstat[j] != QS_COL_BSTAT_BASIC) {
                j++;
            }
            if (j == nstruct) {
                /* No basic column to match the non-basic row */
                fprintf (stderr, "No basic column to match non-basic row %d\n",
                                  i);
                rval = 1; goto CLEANUP;
            }

            if (rstat[i] == QS_ROW_BSTAT_LOWER) {
                fprintf (out, " XL %s %s\n", qslp->colnames[j],
                                             qslp->rownames[i]);
            } else {
                fprintf (out, " XU %s %s\n", qslp->colnames[j],
                                             qslp->rownames[i]);
            }
            i++;
            j++;
        }
    } while (i < nrows);

    /* Now go through and output the non-basic cols at upper bound */

    for (j = 0; j < nstruct; j++) {
        if (cstat[j] == QS_COL_BSTAT_UPPER) {
            fprintf (out, " UL %s\n", qslp->colnames[j]);
        }
    }

    fprintf (out, "ENDATA\n");
    
CLEANUP:

    if (out) fclose (out);
    if (!B) {
        ILL_IFFREE (cstat, char);
        ILL_IFFREE (rstat, char);
    }
    ILL_RETURN (rval, "ILLlib_writebasis");
}

int ILLlib_getrownorms (lpinfo *lp, price_info *pinf, double *rownorms)
{
    int rval = 0;
    int i, j, basic = 0;
    ILLlpdata *qslp = lp->O;
    int nstruct = lp->O->nstruct;
    int nrows   = lp->O->nrows;

    check_pinf (pinf, &rval);
    if (rval) {
/*
        fprintf (stderr, "dual steepest edge norms not available\n");
*/
        ILL_CLEANUP;
    }

    for (i = 0; i < nstruct; i++) {
        j = qslp->structmap[i];
        if (lp->vstat[j] == STAT_BASIC) {
            rownorms[basic++] = pinf->dsinfo.norms[lp->vindex[j]];
        }
    }
    for (i = 0; i < nrows; i++) {
        j = qslp->rowmap[i];
        if (lp->vstat[j] == STAT_BASIC) {
            rownorms[basic++] = pinf->dsinfo.norms[lp->vindex[j]];
        }
    }

    if (basic != nrows) {
        fprintf (stderr, "error in ILLlib_getrownorms\n");
        rval = 1; ILL_CLEANUP;
    }

CLEANUP:

/*
    ILL_RETURN (rval, "ILLlib_getrownorms");
*/
    return rval;   /* Don't want error message */
}

int ILLlib_loadrownorms (lpinfo *lp, price_info *pinf, double *rownorms)
{
    int rval = 0;

    rval = ILLprice_load_rownorms (lp, rownorms, pinf);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_loadrownorms");
}

int ILLlib_recompute_rownorms (lpinfo *lp, price_info *pinf)
{
    int rval = 0;

    rval = ILLprice_build_pricing_info (lp, pinf, DUAL_PHASEII);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlib_recompute_rownorms");
}

int ILLlib_iter (lpinfo *lp)
{
    int iter = 0;

    if (lp && lp->cnts) {
        iter = lp->cnts->pI_iter + lp->cnts->pII_iter +
               lp->cnts->dI_iter + lp->cnts->dII_iter;
    }

    return iter;
}

#define PRINT_TOL 0.000001

int ILLlib_print_x (FILE *fd, lpinfo *lp, ILLlp_cache *C, double *x,
        int nonZerosOnly)
{
    int rval = 0;
    int j;
    ILLlpdata *qslp = lp->O;
    double *dx, *myx = (double *) NULL;

    /* If x is not specified, grab the LP solution */

    if (!x) {
        ILL_SAFE_MALLOC (myx, lp->ncols, double);
        rval = ILLlib_get_x (lp, C, myx);
        ILL_CLEANUP_IF(rval);
        dx = myx;
    } else {
        dx = x;
    }

    fprintf (fd, "Solution Values\n");
    for (j = 0; j < qslp->nstruct; j++) {
        if (!nonZerosOnly || dx[j] > PRINT_TOL || dx[j] < -PRINT_TOL) {
            ILL_FAILfalse(qslp->colnames[j] != NULL,
                          "no NULL names PLEASE!");
            fprintf (fd, "%s = %.6f\n", qslp->colnames[j], dx[j]);
            fflush (fd);
        }
    }

 CLEANUP:

    ILL_IFFREE (myx, double);
    ILL_RETURN(rval, "ILLlib_print_x");
}

int ILLlib_findName(ILLlpdata *qslp, int forRow, const char *name, int id,
        char buf[ILL_namebufsize]) 
{
    ILLsymboltab *tab;
    const char* mode;
    const char *p1, *p2; 
    int sind, rval = 0; 

    id ++; 
    tab =  (forRow) ? &qslp->rowtab : &qslp->coltab;
    if (tab->tablesize == 0) ILLsymboltab_create(tab, 100); 
    p1  =  (forRow) ? "c"  : "x"; 
    p2  =  (forRow) ? "c_" : "x_"; 
    mode = (forRow) ? "row" : "column"; 
    if (name == (char *) NULL) {
        ILLsymboltab_unique_name(tab, id, p1, buf); 
        /*
        fprintf(stderr, "Generating %s name \"%s\".\n", mode, buf); 
        */
    }  else {
       strcpy(buf, name);
    } 
    if (!ILLsymboltab_lookup(tab, buf, &sind)) {
	rval = ILLsymboltab_uname(&qslp->rowtab, buf, p1, p2); 
	if (name != NULL) {
	    fprintf(stderr, "Changing %s name \"%s\" to \"%s\".\n", mode, name, buf); 
	} 
        ILL_CLEANUP_IF(rval); 
    }
CLEANUP: 
    ILL_RETURN(rval, "findName"); 
}

int ILLwrite_lp_file (ILLlpdata *lp, FILE *out, qserror_collector *c)
{
	int rval = 0; 
	qsstring_reporter rep; 

	ILLstring_reporter_copy(&rep, &lp->reporter); 
	ILLstring_reporter_init(&lp->reporter, 
                                (qsreport_string_fct) fprintf, out); 
	rval = ILLwrite_lp(lp, c);
	ILLstring_reporter_copy(&lp->reporter, &rep); 
	return rval;
}

static void check_pinf (price_info *pinf, int *it_exists)
{
    if (!pinf || pinf->dI_price != QS_PRICE_DSTEEP ||
         pinf->dII_price != QS_PRICE_DSTEEP ||
         pinf->dsinfo.norms == (double *) NULL) {
        *it_exists = 1;
    } else {
        *it_exists = 0;
    }
}

#if 0
static int test_matrix (ILLmatrix *A)
{
    int rval = 0;
    int i, j, k;
    int ncols = A->matcols;
    int nrows = A->matrows;
    int matsize = A->matsize;
    int *mbeg = A->matbeg;
    int *mcnt = A->matcnt;
    int *mind = A->matind;
    int *tempi = (int *) NULL;


    if (matsize == 0) ILL_CLEANUP;

    ILL_SAFE_MALLOC (tempi, matsize, int);
    for (i = 0; i < matsize; i++) tempi[i] = 0;

    for (i = 0; i < ncols; i++) {
        for (j = 0; j < mcnt[i]; j++) {
            k = mind[mbeg[i] + j];
            if (k < 0 || k >= nrows) {
                printf ("ERROR IN MATRIX: %d\n", k);
                printf ("ncols = %d, bad col = %d\n", ncols, i);
                printf ("bad cnt = %d, bad index = %d\n", mcnt[i], mbeg[i] + j);
                printf ("matcolsize = %d, matsize = %d\n", A->matcolsize,
                                                         A->matsize);
                rval = 1; ILL_CLEANUP;
            }
            if (tempi[mbeg[i] + j] != 0) {
                printf ("ERROR: over written matrix\n");
                printf ("ncols = %d, bad col = %d\n", ncols, i);
                printf ("nrows = %d\n", nrows);
                printf ("bad cnt = %d, bad index = %d\n", mcnt[i], mbeg[i] + j);
                rval = 1; ILL_CLEANUP;
            } else {
                tempi[mbeg[i] + j] = 1;
            }
        }
    }

    for (i = A->matsize - A->matfree; i < A->matsize; i++) {
        if (tempi[i] != 0) {
            printf ("ERROR: free space is being used\n");
            rval = 1; ILL_CLEANUP;
        }
    }

CLEANUP:

    ILL_IFFREE (tempi, int);
    ILL_RETURN (rval, "test_matrix");
}
#endif
