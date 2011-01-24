/* RCS_INFO = "$RCSfile: binary.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

/****************************************************************************/
/*                                                                          */
/*                     Simple MIP Code to test LP Solver                    */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLmip_bfs (lpinfo *lp, double *val, double *x)                   */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"
#include "simplex.h"
#include "binary.h"
#include "price.h"
#include "lib.h"
#include "qstruct.h"
#include "qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define  ILL_INTTOL (0.000001)

#define  STRONG_PIVOTS     (50)
#define  STRONG_CANDIDATES (10)

#define ILL_BRANCH_STRONG_WEIGHT (10.0)
#define ILL_BRANCH_STRONG_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (ILL_BRANCH_STRONG_WEIGHT * (v0) + (v1))         \
                  : (ILL_BRANCH_STRONG_WEIGHT * (v1) + (v0)))        \
                    / (ILL_BRANCH_STRONG_WEIGHT + 1.0))

#define ILL_BRANCH_PENALTY_WEIGHT (2.0)
#define ILL_BRANCH_PENALTY_VAL(v0,v1,f)                              \
    (((v0)*(f) < (v1)*(1.0-(f)) ?                                    \
        (ILL_BRANCH_PENALTY_WEIGHT * (v0)*(f) + (v1)*(1.0-(f)))    \
      : (ILL_BRANCH_PENALTY_WEIGHT * (v1)*(1.0-(f)) + (v0)*(f)))    \
                    / (ILL_BRANCH_PENALTY_WEIGHT + 1.0))



#define FIRSTBRANCH  1
#define MIDDLEBRANCH 2
#define STRONGBRANCH 3
#define PENALTYBRANCH 4


typedef struct bbnode {
    struct bbnode *next;
    struct bbnode *prev;
    int            id;
    int            depth;
    int            handle;
    double         bound;
    char          *cstat;
    char          *rstat;
    double        *rownorms;
    int            bound_cnt;
    int           *bound_indx;
    char          *lu;
    double        *bounds;
} bbnode;

typedef struct mipinfo {
    int         branching_rule;
    int         watch;
    int         depth;
    int         totalnodes;
    int         activenodes;
    int         totalpivots;
    int         lastpivots;
    int         objsense;
    double      objectivebound;
    double      value;
    double     *downpen;
    double     *uppen;
    double     *x;
    double     *bestx;
    double     *orig_lower;
    double     *orig_upper;
    double     *lower;
    double     *upper;
    lpinfo     *lp;
    price_info *pinf;
    bbnode      head_bbnode;
    ILLpriority *que;
    ILLptrworld ptrworld;
} mipinfo;


ILL_PTRWORLD_ROUTINES(bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE(bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (bbnode, bbnode_check_leaks, depth, int)


static void
    cleanup_mip (mipinfo *minf),
    choose_initial_price (price_info *pinf),
    best_bbnode (mipinfo *minf, bbnode **best),
    put_bbnode (mipinfo *minf, bbnode *b),
    remove_bbnode (bbnode *b),
    find_first_branch (lpinfo *lp, double *x, int *bvar),
    find_middle_branch (lpinfo *lp, double *x, int *bvar),
    check_integral (lpinfo *lp, double *x, int *yesno),
    copy_x (int nstruct, double *from_x, double *to_x),
    init_mipinfo (mipinfo *minf),
    free_mipinfo (mipinfo *minf),
    init_bbnode (bbnode *b),
    free_bbnode (bbnode *b);

static int
    startup_mip (mipinfo *minf, lpinfo *lp, price_info *pinf, double *lpval),
    run_bfs (mipinfo *minf),
    process_bfs_bbnode (mipinfo *minf, bbnode *b),
    child_work (mipinfo *minf, bbnode *active, int bvar, char bdir,
        double *cval, int *cp),
    fix_variables (lpinfo *lp, double bestval, bbnode *b, double *wupper,
        double *wlower, int *hit),
    find_branch (mipinfo *minf, double *x, double lpval, int *bvar),
    find_penalty_branch (lpinfo *lp, price_info *pinf, double *x,
        double *downpen, double *uppen, double lpval, int *bvar),
    find_strong_branch (lpinfo *lp, price_info *pinf, double *x, int *bvar),
    plunge (mipinfo *minf),
    plunge_work (mipinfo *minf, int depth),
    round_variables (mipinfo *minf, int *count, double tol);


static void choose_initial_price (price_info *pinf)
{
   pinf->pI_price  = QS_PRICE_PSTEEP;
   pinf->pII_price = QS_PRICE_PSTEEP;
   pinf->dI_price  = QS_PRICE_DSTEEP;
   pinf->dII_price = QS_PRICE_DSTEEP;
}

int ILLmip_bfs (lpinfo *lp, double *val, double *x)
{
    int tval, rval = 0;
    price_info pinf;
    mipinfo minf;
    bbnode *b;
    double lpval;
    double szeit = ILLutil_zeit ();

    ILLprice_init_pricing_info (&pinf);
    init_mipinfo (&minf);

    if (!lp) {
        fprintf (stderr, "ILLmip_bfs called without an LP\n");
        rval = 1; goto CLEANUP;
    }

    rval = startup_mip (&minf, lp, &pinf, &lpval);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (minf.que, 1, ILLpriority);
    rval = ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
    ILL_CLEANUP_IF (rval);

    b = bbnodealloc (&minf.ptrworld);
    init_bbnode (b);
    b->depth = 0;
    b->id = minf.totalnodes++;
    b->bound = lpval;
    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
    rval = ILLlib_getbasis (lp, b->cstat, b->rstat);
    ILL_CLEANUP_IF (rval);

    if (pinf.dII_price == QS_PRICE_DSTEEP) {
        ILL_SAFE_MALLOC (b->rownorms, lp->nrows, double);
        tval = ILLlib_getrownorms (lp, &pinf, b->rownorms);
        if (tval) {
            printf ("Row norms not available\n"); fflush (stdout);
            ILL_IFFREE (b->rownorms, double);
        }
    }

    rval = ILLutil_priority_insert (minf.que, (void *) b, lpval, &(b->handle));
    ILL_CLEANUP_IF (rval);

    b->prev = &(minf.head_bbnode);
    b->next = (bbnode *) NULL;
    minf.head_bbnode.next = b;
    minf.activenodes++;

    minf.branching_rule = PENALTYBRANCH;

    rval = run_bfs (&minf);
    ILL_CLEANUP_IF (rval);

    printf ("Total Number of Nodes: %d\n", minf.totalnodes);
    printf ("Total Number of Pivots: %d\n", minf.totalpivots);
    printf ("BFS MIP Runing Time: %.2f seconds\n", ILLutil_zeit () - szeit);
    fflush (stdout);

    if (minf.objsense == ILL_MAX) *val = -minf.value;
    else                          *val = minf.value;

    if (x && minf.value != ILL_MAXDOUBLE) {
        copy_x (lp->O->nstruct, minf.bestx, x);
    }

CLEANUP:

    if (minf.que) {
        ILLutil_priority_free (minf.que);
        ILL_IFFREE (minf.que, ILLpriority);
    }
    cleanup_mip (&minf);
    free_mipinfo (&minf);
    ILLprice_free_pricing_info (&pinf);
    ILL_RETURN (rval, "ILLmip_bfs");
}

static int startup_mip (mipinfo *minf, lpinfo *lp, price_info *pinf,
        double *lpval)
{
    int rval = 0;
    int i, col, status, intcount = 0;
    double val;
    ILLlpdata *qlp;

    choose_initial_price (pinf);

    qlp = lp->O;

    rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, pinf, DUAL_SIMPLEX,
                           &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += ILLlib_iter (lp);

    rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &val);
    ILL_CLEANUP_IF (rval);

    printf ("LP Value: %.6f\n", val); fflush (stdout);
    if (lpval) *lpval = val;

    if (qlp->intmarker) {
        for (i = 0; i < qlp->nstruct; i++) {
            if (qlp->intmarker[i]) {
                col = qlp->structmap[i];
                intcount++;
                if (qlp->lower[col] == -ILL_MAXDOUBLE ||
                    qlp->upper[col] ==  ILL_MAXDOUBLE) {
                    printf ("Instance has unbounded integer variable\n");
                    fflush (stdout);
                    rval = 1; goto CLEANUP;
                }
            }
        }
    }

    if (intcount == 0) {
        printf ("No integer variables\n"); fflush (stdout);
        rval = 1; goto CLEANUP;
    } else {
        printf ("%d integer variables\n", intcount); fflush (stdout);
    }

    if (qlp->sinfo) {   /* Free the presolve LP and work with orginal */
        ILLlp_sinfo_free (qlp->sinfo);
        ILL_IFFREE (qlp->sinfo, ILLlp_sinfo);
    }


    minf->lp = lp;
    minf->pinf = pinf;
    minf->objsense = qlp->objsense;
    if (qlp->objsense == ILL_MAX) {    /* MIP codes work with min */
        for (i = 0; i < lp->ncols; i++) {
            qlp->obj[i] = -qlp->obj[i];
        }
        qlp->objsense = ILL_MIN;
    }

    ILL_SAFE_MALLOC (minf->x, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->bestx, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->lower, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->upper, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->orig_lower, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->orig_upper, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->downpen, qlp->nstruct, double);
    ILL_SAFE_MALLOC (minf->uppen, qlp->nstruct, double);

    rval = ILLlib_get_x (lp, (ILLlp_cache *) NULL, minf->x);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
        minf->lower[i] = qlp->lower[i];
        minf->upper[i] = qlp->upper[i];
        minf->orig_lower[i] = qlp->lower[i];
        minf->orig_upper[i] = qlp->upper[i];
        minf->downpen[i] = -1.0;
        minf->uppen[i] = -1.0;
    }


CLEANUP:

    ILL_RETURN (rval, "startup_mip");
}

static void cleanup_mip (mipinfo *minf)
{
    int i;
    ILLlpdata *qslp;
    
    if (minf && minf->lp) {
        qslp = minf->lp->O;
        if (minf->objsense == ILL_MAX) {
            for (i = 0; i < minf->lp->ncols; i++) {
                qslp->obj[i] = -qslp->obj[i];
            }
            qslp->objsense = ILL_MIN;
        }
    }
}

static int run_bfs (mipinfo *minf)
{
    int rval = 0;
    bbnode *b;

    while (minf->head_bbnode.next) {
        best_bbnode (minf, &b);
        rval = process_bfs_bbnode (minf, b);
        ILL_CLEANUP_IF (rval);
        remove_bbnode (b);
        free_bbnode (b); 
        bbnodefree (&minf->ptrworld, b);
        minf->activenodes--;
    }

CLEANUP:

    ILL_RETURN (rval, "run_bfs");
}

static int process_bfs_bbnode (mipinfo *minf, bbnode *active)
{
    lpinfo *lp = minf->lp;
    ILLlp_basis B;
    int status, bvar;
    int i, j, hit, dnp = 0, upp = 0;
    int nstruct = lp->O->nstruct;
    double t, lpval, dnval, upval;
    double *wupper = (double *) NULL;
    double *wlower = (double *) NULL;
    int rval = 0;

    ILLlp_basis_init (&B);

    if (minf->watch > 1) {
        printf ("Node %4d: %.3f", active->id, active->bound);
        if (minf->value != ILL_MAXDOUBLE) printf (" %.3f", minf->value);
        else                              printf ("  None");
        printf (", Active %d ", minf->activenodes);
        fflush (stdout);
    } else if (minf->watch == 1) {
        if (minf->lastpivots > 1000) {
            minf->lastpivots = 0;
            printf ("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
                   minf->totalpivots, minf->activenodes, active->bound);
            if (minf->value != ILL_MAXDOUBLE) printf ("%.3f", minf->value);
            else                              printf ("None\n");
        }
    }

    if (active->bound >= minf->objectivebound) {
        if (minf->watch > 1) {
            printf ("  Node can be purged\n"); fflush (stdout);
        }
        goto CLEANUP;
    }

    /*  Set the LP bounds for the node. */

    ILL_SAFE_MALLOC (wlower, nstruct, double);
    ILL_SAFE_MALLOC (wupper, nstruct, double);

    for (i = 0; i < nstruct; i++) {
        wlower[i] = minf->orig_lower[i];
        wupper[i] = minf->orig_upper[i];
    }
    for (i = 0; i < active->bound_cnt; i++) {
        j = active->bound_indx[i];
        if (active->lu[i] == 'L') wlower[j] = active->bounds[i];
        else                      wupper[j] = active->bounds[i];
    }

    if (active->bound_cnt > 0) {
        rval = ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
                                   active->lu, active->bounds);
        ILL_CLEANUP_IF (rval);
    }

    /*  Solve the LP. */

    rval = ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
                             active->rstat);
    ILL_CLEANUP_IF (rval);
    if (active->rownorms) {
        ILL_SAFE_MALLOC (B.rownorms, lp->nrows, double);
        for (i = 0; i < lp->nrows; i++) {
            B.rownorms[i] = active->rownorms[i];
        }
    }

    rval = ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += ILLlib_iter (lp);
    minf->lastpivots += ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
        printf ("Simplex did not solve the LP\n"); fflush (stdout);
        rval = 1; ILL_CLEANUP;
    }

    if (status == QS_LP_INFEASIBLE) {
        printf ("  Infeasible LP, should have been purged earlier\n");
        fflush (stdout);
        rval = 1; ILL_CLEANUP;
    }

    if (active->depth < 0) {
        for (i = 0; i < nstruct; i++) {
            minf->lower[i] = wlower[i];
            minf->upper[i] = wupper[i];
        }
        rval = plunge (minf);
        ILL_CLEANUP_IF (rval);
    }

    /*  Fix variables. */

    if (minf->value < ILL_MAXDOUBLE) {
        rval = fix_variables (lp, minf->value, active, wupper, wlower, &hit);
        ILL_CLEANUP_IF (rval);

        if (hit) {
            rval = ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status,
                                    0);
            ILL_CLEANUP_IF (rval);

            minf->totalpivots += ILLlib_iter (lp);
            minf->lastpivots += ILLlib_iter (lp);

            if (status == QS_LP_UNSOLVED) {
                printf ("Simplex did not solve the LP\n"); fflush (stdout);
                rval = 1; ILL_CLEANUP;
            }

            if (status == QS_LP_INFEASIBLE) {
                printf ("  Infeasible LP after fixing\n");
                fflush (stdout);
                rval = 1; ILL_CLEANUP;
            }
        }
    }


    /*  Branch. */

    rval = ILLlib_get_x (lp, (ILLlp_cache *) NULL, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &lpval);
    ILL_CLEANUP_IF (rval);

    rval = find_branch (minf, minf->x, lpval, &bvar);
    ILL_CLEANUP_IF (rval);

    if (bvar == -1) {
        printf ("Found integral solution: %f\n", lpval);
        if (lpval < minf->value) {
            minf->value = lpval;
            minf->objectivebound = lpval - ILL_INTTOL;
            copy_x (nstruct, minf->x, minf->bestx);
        }
    } else {
        /* Create down child */

        rval = child_work (minf, active, bvar, 'D', &dnval, &dnp);
        ILL_CLEANUP_IF (rval);

        /* Restore parent basis */

        rval = ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
                             active->rstat);
        ILL_CLEANUP_IF (rval);
        if (active->rownorms) {
            ILL_SAFE_MALLOC (B.rownorms, lp->nrows, double);
            for (i = 0; i < lp->nrows; i++) {
                B.rownorms[i] = active->rownorms[i];
            }
        }

        /* Create up child */

        rval = child_work (minf, active, bvar, 'U', &upval, &upp);
        ILL_CLEANUP_IF (rval);

        if (minf->watch > 1) {
            if (dnval == ILL_MAXDOUBLE) {
                printf ("DN->XXX");
            } else {
                printf ("DN->%.3f%c", dnval, dnp ? 'X' : ' ');
            }
            if (upval == ILL_MAXDOUBLE) {
                printf ("UP->XXX\n");
            } else {
                printf ("UP->%.3f%c\n", upval, upp ? 'X' : ' ');
            }
            fflush (stdout);
        }
    }
    
    /* Set the LP bounds back to original values */

    for (i = 0; i < active->bound_cnt; i++) {
        if (active->lu[i] == 'L') t = minf->orig_lower[active->bound_indx[i]];
        else                      t = minf->orig_upper[active->bound_indx[i]];

        rval = ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_IFFREE (wlower, double);
    ILL_IFFREE (wupper, double);
    ILLlp_basis_free (&B);
    ILL_RETURN(rval, "process_bfs_bbnode");
}

static int child_work (mipinfo *minf, bbnode *active, int bvar, char bdir,
        double *cval, int *cp)
{
    int tval, rval = 0;
    int i, status, intsol;
    double t, oldt, lpval;
    double xi = minf->x[bvar];
    lpinfo *lp = minf->lp;
    bbnode *b;

    *cp = 0;

    if (bdir == 'D') {
        rval = ILLlib_getbnd (lp, bvar, 'U', &oldt);
        ILL_CLEANUP_IF (rval);
        t = ILLutil_our_floor (xi);
        rval = ILLlib_chgbnd (lp, bvar, 'U', t);
        ILL_CLEANUP_IF (rval);
    } else {
        rval = ILLlib_getbnd (lp, bvar, 'L', &oldt);
        ILL_CLEANUP_IF (rval);
        t = ILLutil_our_ceil (xi);
        rval = ILLlib_chgbnd (lp, bvar, 'L', t);
        ILL_CLEANUP_IF (rval);
    }

    rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf, DUAL_SIMPLEX,
                            &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += ILLlib_iter (lp);
    minf->lastpivots += ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
        printf ("Simplex did not solve Child LP\n"); fflush (stdout);
        rval = 1; ILL_CLEANUP;
    }

    if (status == QS_LP_INFEASIBLE) {
        *cval = ILL_MAXDOUBLE; *cp = 1;
    } else {
        rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &lpval);
        ILL_CLEANUP_IF (rval);
        *cval = lpval;

        /* What about the x vector?  Bico - 020531 */

        check_integral (lp, minf->x, &intsol);
        if (intsol) {
            if (lpval < minf->value) {
                printf ("Found integral solution: %f\n", lpval);
                minf->value = lpval;
                minf->objectivebound = lpval - ILL_INTTOL;
                copy_x (lp->O->nstruct, minf->x, minf->bestx);
            }
        }

        if (lpval >= minf->objectivebound) {
            *cp = 1;
        } else {
            b = bbnodealloc (&minf->ptrworld);
            init_bbnode (b);
            b->depth = active->depth + 1;
            b->id = minf->totalnodes;
            b->bound = lpval;
            ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
            ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
            rval = ILLlib_getbasis (lp, b->cstat, b->rstat);
            ILL_CLEANUP_IF (rval);
            if (minf->pinf->dII_price == QS_PRICE_DSTEEP) {
                ILL_SAFE_MALLOC (b->rownorms, lp->nrows, double);
                tval = ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
                if (tval) {
                    printf ("Row norms not available\n"); fflush (stdout);
                    printf ("A\n");
                    exit (1);
                    ILL_IFFREE (b->rownorms, double);
                }
            }
            ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
            ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);
            ILL_SAFE_MALLOC (b->bounds, active->bound_cnt + 1, double);
            for (i = 0; i < active->bound_cnt; i++) {
                b->bound_indx[i] = active->bound_indx[i];
                b->lu[i] = active->lu[i];
                b->bounds[i] = active->bounds[i];
            }
            b->bound_indx[active->bound_cnt] = bvar;
            if (bdir == 'D') b->lu[active->bound_cnt] = 'U';
            else             b->lu[active->bound_cnt] = 'L';
            b->bounds[active->bound_cnt] = t;
            b->bound_cnt = active->bound_cnt + 1;

            rval = ILLutil_priority_insert (minf->que, (void *) b, lpval,
                                            &(b->handle));
            ILL_CLEANUP_IF (rval);
    
            put_bbnode (minf, b);
            minf->activenodes++;
        }
    }
    minf->totalnodes++;

    if (bdir == 'D') {
        rval = ILLlib_chgbnd (lp, bvar, 'U', oldt);
        ILL_CLEANUP_IF (rval);
    } else {
        rval = ILLlib_chgbnd (lp, bvar, 'L', oldt);
        ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    return rval;
}

static int fix_variables (lpinfo *lp, double bestval, bbnode *b, double *wupper,
        double *wlower, int *hit)
{
    int rval = 0;
    int i, nnew = 0;
    int nstruct = lp->O->nstruct;
    double delta, lpval;
    int *new_indx = (int *) NULL;
    char *new_lu = (char *) NULL;
    double *new_bounds = (double *) NULL;
    double *dj = (double *) NULL;

    *hit = 0;

    if (bestval < ILL_MAXDOUBLE) {
        rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &lpval);
        ILL_CLEANUP_IF (rval);
        delta = bestval - lpval + ILL_INTTOL;

        ILL_SAFE_MALLOC (dj, nstruct, double);
        ILL_SAFE_MALLOC (new_indx, nstruct, int);
        ILL_SAFE_MALLOC (new_lu, nstruct, char);
        ILL_SAFE_MALLOC (new_bounds, nstruct, double);

        rval = ILLlib_solution (lp, (ILLlp_cache *) NULL, (double *) NULL,
                   (double *) NULL, (double *) NULL, (double *) NULL, dj);
        ILL_CLEANUP_IF (rval);

        for (i = 0; i < nstruct; i++) {
            if (lp->O->intmarker[i]) {
                if (wlower[i] != wupper[i]) {
                    if (dj[i] > delta) {
                        wupper[i] -= 1.0;
                        rval = ILLlib_chgbnd (lp, i, 'U', wupper[i]);
                        ILL_CLEANUP_IF (rval);
                        new_indx[nnew] = i;
                        new_lu[nnew] = 'U';
                        new_bounds[nnew] = wupper[i];
                        nnew++;
                    } else if (-dj[i] > delta) {
                        wlower[i] += 1.0;
                        rval = ILLlib_chgbnd (lp, i, 'L', wlower[i]);
                        ILL_CLEANUP_IF (rval);
                        new_indx[nnew] = i;
                        new_lu[nnew] = 'L';
                        new_bounds[nnew] = wlower[i];
                        nnew++;
                    }
                }
            }
        }

        if (nnew) {
            rval = ILLutil_reallocrus_count ((void **) &(b->bound_indx),
                          b->bound_cnt + nnew, sizeof (int));
            ILL_CLEANUP_IF (rval);
            rval = ILLutil_reallocrus_count ((void **) &(b->lu),
                          b->bound_cnt + nnew, sizeof (char));
            ILL_CLEANUP_IF (rval);
            rval = ILLutil_reallocrus_count ((void **) &(b->bounds),
                          b->bound_cnt + nnew, sizeof (double));
            ILL_CLEANUP_IF (rval);
            for (i = 0; i < nnew; i++) {
                b->bound_indx[b->bound_cnt + i] = new_indx[i];
                b->lu[b->bound_cnt + i]         = new_lu[i];
                b->bounds[b->bound_cnt + i]     = new_bounds[i];
            }
            b->bound_cnt += nnew;
        }
    }
  
    *hit = nnew;

CLEANUP:

    ILL_IFFREE (dj, double);
    ILL_IFFREE (new_indx, int);
    ILL_IFFREE (new_lu, char);
    ILL_IFFREE (new_bounds, double);

    return rval;
}

static void best_bbnode (mipinfo *minf, bbnode **best)
{
#if 0
    bbnode *b;
    double bestval = ILL_MAXDOUBLE;

    for (b = minf->head_bbnode.next; b; b = b->next) {
        if (b->bound < bestval) {
            *best = b;
            bestval = b->bound;
        }
    }
#endif

    double val;
    ILLutil_priority_deletemin (minf->que, &val, (void **) best);
}

static void put_bbnode (mipinfo *minf, bbnode *b)
{
    b->next = minf->head_bbnode.next;
    b->prev = &(minf->head_bbnode);
    if (b->next) b->next->prev = b;
    minf->head_bbnode.next = b;
}

static void remove_bbnode (bbnode *b)
{
    b->prev->next = b->next;
    if (b->next) b->next->prev = b->prev;
}

static int find_branch (mipinfo *minf, double *x, double lpval, int *bvar)
{
    lpinfo *lp = minf->lp;
    int rval = 0;

    switch (minf->branching_rule) {
    case PENALTYBRANCH:
        rval = find_penalty_branch (lp, minf->pinf, x, minf->downpen,
                                    minf->uppen, lpval, bvar);
        ILL_CLEANUP_IF (rval);
        break;
    case FIRSTBRANCH:
        find_first_branch (lp, x, bvar);
        break;
    case MIDDLEBRANCH:
        find_middle_branch (lp, x, bvar);
        break;
    case STRONGBRANCH:
        rval = find_strong_branch (lp, minf->pinf, x, bvar);
        ILL_CLEANUP_IF (rval);
        break;
    default:
        fprintf (stderr, "Unknown branching rule.\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "find_branch");
}

static void find_first_branch (lpinfo *lp, double *x, int *bvar)
{
    int i, ibest = -1;
    ILLlpdata *qslp = lp->O;
    double t;

    for (i = 0; i < qslp->nstruct; i++) {
        if (qslp->intmarker[i]) {
            t = ILLutil_our_frac (x[i]);
            if (t > ILL_INTTOL && t < 1.0 - ILL_INTTOL) {
                ibest = i; 
                break;
            }
        }
    }
    *bvar = ibest;
}

static void find_middle_branch (lpinfo *lp, double *x, int *bvar)
{
    int i, ibest = -1;
    double t, tbest = 0.5;
    ILLlpdata *qlp = lp->O;

    for (i = 0; i < qlp->nstruct; i++) {
        if (qlp->intmarker[i]) {
            t = ILLutil_our_frac (x[i]) - 0.5;
            if (t < 0.0)
                t = -t;
            if (t < tbest) {
                tbest = t;
                ibest = i;
            }
        }
    }

    if (tbest < (0.5 - ILL_INTTOL)) {
        *bvar = ibest;
    } else {
        *bvar = -1;
    }
}

static int find_penalty_branch (lpinfo *lp, price_info *pinf, double *x,
        double *downpen, double *uppen, double lpval, int *bvar)
{
    int rval = 0;
    int i, k, ibest = -1, ncand = 0, nneed = 0;
    ILLlpdata *qslp = lp->O;
    int *candidatelist = (int *) NULL;
    int *needlist = (int *) NULL;
    double *fval = (double *) NULL;
    double *xlist = (double *) NULL;
    double *newdown = (double *) NULL;
    double *newup = (double *) NULL;
    double t, tbest = -ILL_MAXDOUBLE;

    ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
    ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);
    ILL_SAFE_MALLOC (fval, qslp->nstruct, double);
    ILL_SAFE_MALLOC (xlist, qslp->nstruct, double);
    for (i = 0; i < qslp->nstruct; i++) {
        if (qslp->intmarker[i]) {
            fval[i] = ILLutil_our_frac (x[i]);
            if (fval[i] > ILL_INTTOL && fval[i] < 1.0 - ILL_INTTOL) {
                candidatelist[ncand++] = i;
                if (downpen[i] == -1.0) {
                    xlist[nneed] = x[i];
                    needlist[nneed++] = i;
                }
            }
        }
    }

    if (nneed > 0) {
        ILL_SAFE_MALLOC (newdown, nneed, double);
        ILL_SAFE_MALLOC (newup, nneed, double);
        rval = ILLlib_strongbranch (lp, pinf, needlist, nneed,
                   (double *) NULL, newdown, newup, 5 * STRONG_PIVOTS,
                   ILL_MAXDOUBLE / 10.0);
        ILL_CLEANUP_IF (rval);

        for (i = 0; i < nneed; i++) {
            k = needlist[i];
            downpen[k] = (newdown[i] - lpval) / fval[k];
            uppen[k]   = (newup[i] - lpval) / (1.0 - fval[k]);
            
        }
    }

    for (i = 0; i < ncand; i++) {
        k = candidatelist[i];
        t = ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]);
        if (t > tbest) {
            tbest = t;
            ibest = k;
        }
    }

    *bvar = ibest;

CLEANUP:

    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (needlist, int);
    ILL_IFFREE (newdown, double);
    ILL_IFFREE (newup, double);
    ILL_IFFREE (fval, double);
    ILL_IFFREE (xlist, double);
    return rval;
}

static int find_strong_branch (lpinfo *lp, price_info *pinf, double *x,
        int *bvar)
{
    int rval = 0;
    int i, ibest = -1, ncand = 0;
    int maxtrys = STRONG_CANDIDATES;
    double t, tbest = -ILL_MAXDOUBLE;
    ILLlpdata *qlp = lp->O;
    int *candidatelist = (int *) NULL;
    int *newlist       = (int *) NULL;
    int *perm          = (int *) NULL;
    double *tval       = (double *) NULL;
    double *xlist      = (double *) NULL;
    double *downpen    = (double *) NULL;
    double *uppen      = (double *) NULL;
    ILLrandstate rstate;

    ILLutil_sprand (999, &rstate);

    ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);

    ILL_SAFE_MALLOC (tval, qlp->nstruct, double);

    for (i = 0; i < qlp->nstruct; i++) {
        if (qlp->intmarker[i]) {
            t = ILLutil_our_frac (x[i]) - 0.5;
            if (t < 0.0)
                t = -t;
            if (t < (0.5 - ILL_INTTOL)) {
                candidatelist[ncand] = i;
                tval[ncand++]        = t;
            }
        }
    }

    if (ncand > 0) {
        if (ncand > maxtrys) {
            ILL_SAFE_MALLOC (perm, ncand, int);

            for (i = 0; i < ncand; i++) {
                perm[i] = i;
            }
            ILLutil_rselect (perm, 0, ncand-1, maxtrys, tval, &rstate);

            ILL_SAFE_MALLOC (newlist, maxtrys, int);

            for (i = 0; i < maxtrys; i++) {
                newlist[i] = candidatelist[perm[i]];
            }
            ILL_IFFREE (candidatelist, int);
            candidatelist = newlist; 
            newlist = (int *) NULL;
            ncand = maxtrys;
        }

        ILL_SAFE_MALLOC (downpen, ncand, double);
        ILL_SAFE_MALLOC (uppen, ncand, double);
        ILL_SAFE_MALLOC (xlist, ncand, double);

        for (i = 0; i < ncand; i++) {
            xlist[i] = x[candidatelist[i]];
        }
 
        rval = ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
                    (double *) NULL, downpen, uppen, STRONG_PIVOTS,
                    ILL_MAXDOUBLE / 10.0);
        ILL_CLEANUP_IF (rval);

        for (i = 0; i < ncand; i++) {
            t = ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]);
            if (t > tbest) {
                tbest = t;
                ibest = candidatelist[i];
            }
        }
    }

    *bvar = ibest;


CLEANUP:

    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (newlist, int);
    ILL_IFFREE (perm, int);
    ILL_IFFREE (tval, double);
    ILL_IFFREE (downpen, double);
    ILL_IFFREE (uppen, double);
    ILL_IFFREE (xlist, double);

    ILL_RETURN (rval, "find_strong_branch");
}

static void check_integral (lpinfo *lp, double *x, int *yesno)
{
    int i;
    double t;
    ILLlpdata *qlp = lp->O;

    for (i = 0; i < qlp->nstruct; i++) {
        if (qlp->intmarker[i]) {
            t = ILLutil_our_frac (x[i]);
            if (t > ILL_INTTOL && t < 1.0 - ILL_INTTOL) {
                *yesno = 0;
                return;
            }
        }
    }

    *yesno = 1;
}

static int plunge (mipinfo *minf)
{
    int rval = 0;
    int i, status;
    lpinfo *lp = minf->lp;
    ILLlpdata *qlp = minf->lp->O;
    double *oldlower = (double *) NULL;
    double *oldupper = (double *) NULL;

    if (minf->watch) {
        printf ("Plunging ...\n"); fflush (stdout);
    }

    ILL_SAFE_MALLOC (oldlower, qlp->nstruct, double);
    ILL_SAFE_MALLOC (oldupper, qlp->nstruct, double);

    for (i = 0; i < qlp->nstruct; i++) {
        oldlower[i] = minf->lower[i];
        oldupper[i] = minf->upper[i];
    }

    rval = plunge_work (minf, 0);
    ILL_CLEANUP_IF (rval);
    
    for (i = 0; i < qlp->nstruct; i++) {
        rval = ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
        ILL_CLEANUP_IF (rval);
        rval = ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
        ILL_CLEANUP_IF (rval);
        minf->lower[i] = oldlower[i];
        minf->upper[i] = oldupper[i];
    }

    rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf, DUAL_SIMPLEX,
                           &status, 0);
    ILL_CLEANUP_IF (rval);
    

CLEANUP:

    ILL_IFFREE (oldlower, double);
    ILL_IFFREE (oldupper, double);

    ILL_RETURN (rval, "plunge");
}

static int plunge_work (mipinfo *minf, int depth)
{
    int rval = 0;
    int bvar, status, count;
    double lpval, val0, val1;
    lpinfo *lp = minf->lp;

    rval = ILLlib_get_x (lp, (ILLlp_cache *) NULL, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = round_variables (minf, &count, 0.001 /* 0.001 */);
    ILL_CLEANUP_IF (rval);
    if (count) {
        rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf,
                               DUAL_SIMPLEX, &status, 0);
        ILL_CLEANUP_IF (rval);
        if (status != QS_LP_OPTIMAL) {
            goto CLEANUP;
        }
        rval = ILLlib_get_x (lp, (ILLlp_cache *) NULL, minf->x);
        ILL_CLEANUP_IF (rval);
    }
    
    find_middle_branch (lp, minf->x, &bvar);
    if (bvar == -1) {
        rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &lpval);
        ILL_CLEANUP_IF (rval);

        if (lpval < minf->value) {
            printf ("Plunge Integral Solution: %.6f (Depth: %d)\n",
                                                 lpval, depth);
            fflush (stdout);

            minf->value = lpval;
            minf->objectivebound = lpval - ILL_INTTOL;
            copy_x (lp->O->nstruct, minf->x, minf->bestx);
        }
        goto CLEANUP;
    }

    minf->lower[bvar] = 1.0;
    rval = ILLlib_chgbnd (lp, bvar, 'L', 1.0);
    ILL_CLEANUP_IF (rval);
    rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf, DUAL_SIMPLEX,
                           &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
        printf ("Simplex did not solve the plunge LP\n"); fflush (stdout);
        rval = 1; ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
        val1 = ILL_MAXDOUBLE;
    } else if (status == QS_LP_OPTIMAL) {
        rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &val1);
        ILL_CLEANUP_IF (rval);
    } else {
        ILL_CLEANUP; 
    }

    rval = ILLlib_chgbnd (lp, bvar, 'L', 0.0);
    ILL_CLEANUP_IF (rval);
    minf->lower[bvar] = 0.0;

    minf->upper[bvar] = 0.0;
    rval = ILLlib_chgbnd (lp, bvar, 'U', 0.0);
    ILL_CLEANUP_IF (rval);
    rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf, DUAL_SIMPLEX,
                           &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
        printf ("Simplex did not solve the plunge LP\n"); fflush (stdout);
        rval = 1; ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
        val0 = ILL_MAXDOUBLE;
    } else if (status == QS_LP_OPTIMAL) {
        rval = ILLlib_objval (lp, (ILLlp_cache *) NULL, &val0);
        ILL_CLEANUP_IF (rval);
    } else {
        ILL_CLEANUP; 
    }

    rval = ILLlib_chgbnd (lp, bvar, 'U', 1.0);
    ILL_CLEANUP_IF (rval);
    minf->upper[bvar] = 1.0;

    if (val0 == ILL_MAXDOUBLE && val1 == ILL_MAXDOUBLE) {
        ILL_CLEANUP; 
    }

    if (val0 < val1) {
        minf->upper[bvar] = 0.0;
        rval = ILLlib_chgbnd (lp, bvar, 'U', 0.0);
        ILL_CLEANUP_IF (rval);
        rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf,
                               DUAL_SIMPLEX, &status, 0);
        ILL_CLEANUP_IF (rval);
        rval = plunge_work (minf, depth+1);
        ILL_CLEANUP_IF (rval);
        rval = ILLlib_chgbnd (lp, bvar, 'U', 1.0);
        ILL_CLEANUP_IF (rval);
        minf->upper[bvar] = 1.0;
    } else {
        minf->lower[bvar] = 1.0;
        rval = ILLlib_chgbnd (lp, bvar, 'L', 1.0);
        ILL_CLEANUP_IF (rval);
        rval = ILLlib_optimize (lp, (ILLlp_basis *) NULL, minf->pinf,
                               DUAL_SIMPLEX, &status, 0);
        ILL_CLEANUP_IF (rval);
        rval = plunge_work (minf, depth+1);
        ILL_CLEANUP_IF (rval);
        rval = ILLlib_chgbnd (lp, bvar, 'L', 0.0);
        ILL_CLEANUP_IF (rval);
        minf->lower[bvar] = 0.0;
    }

CLEANUP:

    ILL_RETURN (rval, "plunge_work");
}

static int round_variables (mipinfo *minf, int *count, double tol)
{
    int rval = 0;
    int i, hit = 0;
    lpinfo *lp = minf->lp;
    ILLlpdata *qlp = lp->O;

    *count = 0;

    for (i = 0; i < qlp->nstruct; i++) {
        if (qlp->intmarker[i]) {
            if (minf->lower[i] != minf->upper[i]) {
                if (minf->x[i] < tol) {
                    minf->upper[i] = 0.0;
                    rval = ILLlib_chgbnd (lp, i, 'U', 0.0);
                    ILL_CLEANUP_IF (rval);
                    hit++;
                } else if (minf->x[i] > 1.0 - tol) {
                    minf->lower[i] = 1.0;
                    rval = ILLlib_chgbnd (lp, i, 'L', 1.0);
                    ILL_CLEANUP_IF (rval);
                    hit++;
                }
            }
        }
    }
    *count = hit;

CLEANUP:

    ILL_RETURN (rval, "round_variables");
}

static void copy_x (int nstruct, double *from_x, double *to_x) 
{
    int j;

    for (j = 0; j < nstruct; j++) {
        to_x[j] = from_x[j];
    }
}

static void init_mipinfo (mipinfo *minf)
{
    if (minf) {
        minf->branching_rule = /* MIDDLEBRANCH */ STRONGBRANCH;
        minf->watch = 1;
        minf->depth = 0;
        minf->totalnodes = 0;
        minf->activenodes = 0;
        minf->totalpivots = 0;
        minf->lastpivots = 0;
        minf->objectivebound = ILL_MAXDOUBLE;
        minf->value = ILL_MAXDOUBLE;
        minf->downpen = (double *) NULL;
        minf->uppen   = (double *) NULL;
        minf->x     = (double *) NULL;
        minf->bestx = (double *) NULL;
        minf->lower = (double *) NULL;
        minf->upper = (double *) NULL;
        minf->lp   = (lpinfo *) NULL;
        minf->pinf = (price_info *) NULL;
        minf->head_bbnode.prev = (bbnode *) NULL;
        minf->head_bbnode.next = (bbnode *) NULL;
        minf->que = (ILLpriority *) NULL;
        ILLptrworld_init (&minf->ptrworld);
    }
}

static void free_mipinfo (mipinfo *minf)
{
    int total, onlist;

    if (minf) {
        ILL_IFFREE (minf->downpen, double);
        ILL_IFFREE (minf->uppen, double);
        ILL_IFFREE (minf->x, double);
        ILL_IFFREE (minf->bestx, double);
        ILL_IFFREE (minf->lower, double);
        ILL_IFFREE (minf->upper, double);
        bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
        if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist)) {
            fprintf (stderr, "WARNING: %d outstanding bbnodes\n",
                     total - onlist);
        }
        ILLptrworld_delete (&minf->ptrworld);
        init_mipinfo (minf);
    }
}

static void init_bbnode (bbnode *b)
{
    if (b) {
        b->next       = (bbnode *) NULL;
        b->prev       = (bbnode *) NULL;
        b->id         = 0;
        b->depth      = 0;
        b->handle     = 0;
        b->bound      = -ILL_MAXDOUBLE;
        b->cstat      = (char *) NULL;
        b->rstat      = (char *) NULL;
        b->rownorms   = (double *) NULL;
        b->bound_cnt  = 0;
        b->bound_indx = (int *) NULL;
        b->lu         = (char *) NULL;
        b->bounds     = (double *) NULL;
    }
} 

static void free_bbnode (bbnode *b)
{
    if (b) {
        ILL_IFFREE (b->cstat, char);
        ILL_IFFREE (b->rstat, char);
        ILL_IFFREE (b->rownorms, double);
        ILL_IFFREE (b->bound_indx, int);
        ILL_IFFREE (b->lu, char);
        ILL_IFFREE (b->bounds, double);
        init_bbnode (b);
    }
}
