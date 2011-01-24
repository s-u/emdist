/* RCS_INFO = "$RCSfile: presolve.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

/****************************************************************************/
/*                                                                          */
/*                 Presolve Routine for Simplex Method                      */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLlp_add_logicals (ILLlpata *lp)                                 */
/*    int ILLlp_presolve (ILLlpdata *lp)                                    */
/*    int ILLlp_scale (ILLlpdata *lp)                                       */
/*    void ILLlp_sinfo_init (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_sinfo_free (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_predata_init (ILLlp_predata *pre)                          */
/*    void ILLlp_predata_free (ILLlp_predata *pre)                          */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*    presolve will assume that logicals have been added.                   */
/*                                                                          */
/****************************************************************************/

#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"

#define ILL_LP_STATUS_OK (0)
#define ILL_PRE_FEAS_TOL (0.000001)
#define ILL_PRE_ZERO_TOL (0.0000000001)

#define ILL_PRE_DELETE_EMPTY_ROW               (1)
#define ILL_PRE_DELETE_SINGLETON_ROW           (2)
#define ILL_PRE_DELETE_FIXED_VARIABLE          (3)
#define ILL_PRE_DELETE_FORCED_VARIABLE         (4)
#define ILL_PRE_DELETE_SINGLETON_VARIABLE      (5)
#define ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE (6)
#define ILL_PRE_DELETE_EMPTY_COLUMN            (7)

#define ILL_PRE_COL_STRUC                      (0)
#define ILL_PRE_COL_LOGICAL                    (1)

static int debug = 0;

typedef struct edge {
    int               row;
    int               col;
    char              coltype;
    char              mark;
    char              del;
    double            coef;
} edge;

typedef struct node {
    struct edge     **adj;
    double            obj;
    double            lower;
    double            upper;
    double            rhs;
    int               deg;
    char              mark;
    char              del;
    char              coltype;
    char              rowsense;
} node;

typedef struct intptr {
    int              this;
    struct intptr   *next;
} intptr;

typedef struct graph {
    struct edge      *edgelist;
    struct node      *rows;
    struct node      *cols;
    int               ecount;
    int               nrows;
    int               ncols;
    struct edge     **adjspace;
    ILLptrworld        intptrworld;
    int               objsense;
} graph;

void
    ILLlp_sinfo_init (ILLlp_sinfo *sinfo),
    ILLlp_sinfo_free (ILLlp_sinfo *sinfo),
    ILLlp_predata_init (ILLlp_predata *pre),
    ILLlp_predata_free (ILLlp_predata *pre),
    ILLlp_preop_init (ILLlp_preop *op),
    ILLlp_preop_free (ILLlp_preop *op),
    ILLlp_preline_init (ILLlp_preline *line),
    ILLlp_preline_free (ILLlp_preline *line);

int
    ILLlp_sinfo_print (ILLlp_sinfo *s);

static void
    set_fixed_variable (graph *G, int j, double val),
    get_implied_rhs_bounds (graph *G, int i, double *lb, double *ub),
    get_implied_variable_bounds (graph *G, int j, edge *a_ij,
        double *lb, double *ub),
    dump_line (ILLlp_preline *line),
    init_graph (graph *G),
    free_graph (graph *G),
    dump_graph (graph *G);
    
static int
    simple_presolve (ILLlpdata *lp, ILLlp_predata *pre, ILLlp_sinfo *info,
        int pre_types, int *status),
    grab_lp_line (graph *G, int indx, ILLlp_preline *line, int row_or_col),
    grab_lp_info (graph *G, char **colnames, ILLlp_sinfo *info),
    fixed_variables (graph *G, ILLlp_predata *pre),
    empty_columns (graph *G, ILLlp_predata *pre),
    singleton_rows (graph *G, ILLlp_predata *pre, int *hit),
    forcing_constraints (graph *G, ILLlp_predata *pre, int *hit),
    singleton_columns (graph *G, ILLlp_predata *pre, int *hit),
    duplicate_rows (graph *G, int *hit),
    duplicate_cols (graph *G, int *hit),
    gather_dup_lists (int *s, int count, int *duptotal,
        int **dupcnt, int **dupind),
    get_next_preop (ILLlp_predata *pre, ILLlp_preop **op),
    add_to_list (ILLptrworld *world, intptr **list, int i),
    build_graph (ILLlpdata *lp, graph *G);


ILL_PTRWORLD_ROUTINES(intptr, intptralloc, intptr_bulkalloc, intptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE(intptr, intptr_listfree, intptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

int ILLlp_add_logicals (ILLlpdata *lp) 
{
    int rval = 0;
    int ncols, nrows, nzcount, i, aindex;
    char *sense;
    ILLmatrix *A;

    if (!lp) {
        fprintf (stderr, "ILLlp_add_logicals called with a NULL pointer\n");
        rval = 1; goto CLEANUP;
    }

    printf ("ILLlp_add_logicals ...\n"); fflush (stdout);

    A = &lp->A;
    sense = lp->sense;
    ncols = lp->ncols;
    nrows = lp->nrows;
    nzcount = lp->nzcount;

    if (nrows == 0) goto CLEANUP;

    rval = ILLutil_reallocrus_count ((void **) &(lp->obj),
                                    lp->colsize + nrows, sizeof (double));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(lp->upper),
                                    lp->colsize + nrows, sizeof (double));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(lp->lower),
                                    lp->colsize + nrows, sizeof (double));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(lp->colnames),
                                    lp->colsize + nrows, sizeof (char *));
    ILL_CLEANUP_IF (rval);
    for (i = ncols; i < ncols + nrows; i++) {
        lp->colnames[i] = (char *) NULL;
    }

    ILL_SAFE_MALLOC (lp->rowmap, lp->rowsize, int);


    rval = ILLutil_reallocrus_count ((void **) &(A->matcnt),
                                    A->matcolsize + nrows, sizeof (int));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(A->matbeg),
                                    A->matcolsize + nrows, sizeof (int));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(A->matind),
                                    A->matsize + nrows, sizeof (int));
    ILL_CLEANUP_IF (rval);

    rval = ILLutil_reallocrus_count ((void **) &(A->matval),
                                    A->matsize + nrows, sizeof (double));
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < nrows; i++) {
        A->matind[A->matsize + i] = -1;
    }

    aindex = A->matsize - A->matfree;

    for (i = 0; i < nrows; i++) {
        lp->rowmap[i] = ncols;
        lp->obj[ncols] = 0.0;
        A->matcnt[ncols] = 1;
        A->matbeg[ncols] = aindex;
        A->matind[aindex] = i;
        switch (sense[i]) {
        case 'E':                         /* Arificial */
            lp->lower[ncols] = 0.0;       
            lp->upper[ncols] = 0.0;
            A->matval[aindex] = 1.0;
            break;
        case 'G':                         /* Surplus   */
            lp->lower[ncols] = 0.0;       
            lp->upper[ncols] = ILL_MAXDOUBLE;
            A->matval[aindex] = -1.0;
            break;
        case 'L':                         /* Slack     */
            lp->lower[ncols] = 0.0;       
            lp->upper[ncols] = ILL_MAXDOUBLE;
            A->matval[aindex] = 1.0;
            break;
        case 'R':                         /* Range     */
            lp->lower[ncols] = 0.0;
            lp->upper[ncols] = lp->rangeval[i];
            A->matval[aindex] = -1.0;
            break;
        default:
            fprintf (stderr, "unknown sense %c in ILLlp_add_logicals\n",
                                                              sense[i]);
            rval = 1; goto CLEANUP;
        }
        ncols++;
        nzcount++;
        aindex++;
    }

    lp->ncols   = ncols;
    lp->nzcount = nzcount;
    A->matcols  = ncols;

    lp->colsize   += nrows;
    A->matsize    += nrows;
    A->matcolsize += nrows;

    if (lp->rA) {
        ILLlp_rows_clear (lp->rA);
    } else {
        ILL_SAFE_MALLOC (lp->rA, 1, ILLlp_rows);
    }

    rval = ILLlp_rows_init (lp->rA, lp, 1);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "ILLlp_add_logicals");
}

int ILLlp_scale (ILLlpdata *lp) 
{
    int rval = 0;
    int i, j, k, col, row, nstruct, start, stop;
    ILLmatrix *A;
    double rho;
    double *gamma = (double *) NULL;

    /* Columns - divide by largest absolute value */

    if (!lp) {
        ILL_ERROR (rval, "ILLlp_scale called with a NULL pointer");
    }

    if (lp->nrows == 0 || lp->ncols == 0) goto CLEANUP;

    A = &lp->A;
    nstruct = lp->nstruct;

    for (j = 0; j < nstruct; j++) {
        col = lp->structmap[j];
        rho = 0.0;

        start = A->matbeg[col];
        stop = start + A->matcnt[col];

        for (k = start; k < stop; k++) {
            if (rho < fabs (A->matval[k])) {
                rho = fabs (A->matval[k]);
            }
        }

        if (rho > 0.0) {
            for (k = start; k < stop; k++) {
                A->matval[k] /= rho;
            }
            lp->obj[col] /= rho;
            if (lp->lower[col] != -ILL_MAXDOUBLE) lp->lower[col] *= rho;
            if (lp->upper[col] !=  ILL_MAXDOUBLE) lp->upper[col] *= rho;
        }
    }

    ILL_SAFE_MALLOC (gamma, lp->nrows, double);
    for (i = 0; i < lp->nrows; i++) {
        gamma[i] = 0.0;
    }

    for (j = 0; j < nstruct; j++) {
        col = lp->structmap[j];
        start = A->matbeg[col];
        stop = start + A->matcnt[col];

        for (k = start; k < stop; k++) {
            row = A->matind[k];
            if (gamma[row] < fabs (A->matval[k])) {
                gamma[row] = fabs (A->matval[k]);
            }
        }
    }

    for (j = 0; j < nstruct; j++) {
        col = lp->structmap[j];
        start = A->matbeg[col];
        stop = start + A->matcnt[col];

        for (k = start; k < stop; k++) {
            row = A->matind[k];
            if (gamma[row] > 0.0) {
                A->matval[k] /= gamma[row];
            }
        }
    }

    for (i = 0; i < lp->nrows; i++) {
        if (gamma[i] > 0.0) {
            lp->rhs[i] /= gamma[i];
            col = lp->rowmap[i];
            if (lp->upper[col] !=  ILL_MAXDOUBLE) {
                lp->upper[col] /= gamma[i];    /* Ranged row */
            }
        }
    }

    if (lp->rA) {  /* Need to clear the row version of data */
        ILLlp_rows_clear (lp->rA);
        ILL_IFFREE (lp->rA, ILLlp_rows);
    }


CLEANUP:

    ILL_IFFREE (gamma, double);
    ILL_RETURN (rval, "ILLlp_scale");
}

int ILLlp_presolve (ILLlpdata *lp, int pre_types) 
{
    int rval = 0;
    int status        = ILL_LP_STATUS_OK;
    ILLlp_predata *pre = (ILLlp_predata *) NULL;
    ILLlp_sinfo *info  = (ILLlp_sinfo *) NULL;

    if (!lp) {
        fprintf (stderr, "ILLlp_presolve called with a NULL pointer\n");
        rval = 1; goto CLEANUP;
    }


/*
    ILLlpdata_writelp (lp, (char *) NULL);
    printf ("\n"); fflush (stdout);
*/

    ILL_SAFE_MALLOC (pre, 1, ILLlp_predata);
    ILLlp_predata_init (pre);

    ILL_SAFE_MALLOC (info, 1, ILLlp_sinfo);
    ILLlp_sinfo_init (info);

    rval = simple_presolve (lp, pre, info, pre_types, &status);
    ILL_CLEANUP_IF (rval);
    if (status != ILL_LP_STATUS_OK) {
        printf ("simple_presolve returned with bad status\n");
        rval = 1; goto CLEANUP;
    }

/*
    rval = ILLlp_sinfo_print (info);
    ILL_CLEANUP_IF (rval);
*/

CLEANUP:

    if (rval) {
        if (pre) {
            ILLlp_predata_free (pre);
            ILL_IFFREE (pre, ILLlp_predata);
        }

        if (info) {
            ILLlp_sinfo_free (info);
            ILL_IFFREE (info, ILLlp_sinfo);
        }
    } else {
        lp->presolve = pre;
        lp->sinfo    = info;
    }

    ILL_RETURN (rval, "ILLlp_presolve");
}


#if 0
int ILLlp_presolve_addrow (lpinfo *lp, int cnt, int *ind, double *val,
        double rhs)
{
    int rval = 0;
    ILLlpdata *qslp;
    ILLlp_sinfo *S;
    ILLmatrix *A;

    /* This will need to evolve into a function that handles the task */
    /* of working through the presolve data to determine the new LP   */
    /* created when a row is added to the original LP.                */

    /* The copies of the obj and bound used in the simplex code are   */
    /* also updated in this function.                                 */

    if (!lp) {
        fprintf (stderr, "ILLlp_presolve_addrow is called without an LP\n");
        rval = 1; goto CLEANUP;
    }

    if (lp->presolve != (ILLlp_predata *) NULL) {
        fprintf (stderr, "Not yet set up to handle addrows after presolve\n");
        rval = 1; goto CLEANUP;
    }

    qslp = lp->O;
    S = qslp->sinfo;
    A = S->A; 

 
    rval = ILLlib_matrix_addrow (A, cnt, ind, val, rhs);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_RETURN (rval, "ILLlp_presolve_addrow");
}
#endif


static int simple_presolve (ILLlpdata *lp, ILLlp_predata *pre, ILLlp_sinfo *info,
        int pre_types, int *status)
{
    int rval = 0;
    int i, hit, newhit;
    graph G;

    if (status) *status = ILL_LP_STATUS_OK;
    init_graph (&G);
 
    if (!lp) {
        fprintf (stderr, "simple_presolve called with a NULL pointer\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Initial Rows = %d, Cols = %d, Nzcount = %d\n",
               lp->nrows, lp->ncols, lp->nzcount);
    fflush (stdout);

    rval = build_graph (lp, &G);
    ILL_CLEANUP_IF (rval);
    if (debug) dump_graph (&G);

    if (pre_types & ILL_PRE_FIXED) {
        rval = fixed_variables (&G, pre);
        ILL_CLEANUP_IF (rval);
    }

    do {
        hit = 0;
        if (pre_types & ILL_PRE_SINGLE_ROW) {
            rval = singleton_rows (&G, pre, &newhit);
            ILL_CLEANUP_IF (rval);
            hit += newhit;
        }

        if (pre_types & ILL_PRE_FORCING) {
            rval = forcing_constraints (&G, pre, &newhit);
            ILL_CLEANUP_IF (rval);
            hit += newhit;
        }

        if (pre_types & ILL_PRE_SINGLE_COL) {
            rval = singleton_columns (&G, pre, &newhit);
            ILL_CLEANUP_IF (rval);
            hit += newhit;
        }

        if (pre_types & ILL_PRE_DUPLICATE_ROW) {
            rval = duplicate_rows (&G, &newhit);
            ILL_CLEANUP_IF (rval);
            hit += newhit;
        }

        if (pre_types & ILL_PRE_DUPLICATE_COL) {
            rval = duplicate_cols (&G, &newhit);
            ILL_CLEANUP_IF (rval);
            hit += newhit;
        }


/*
        {
            int k, cnt = 0;
            for (i = 0; i < G.ncols; i++) {
                if (G.cols[i].del == 0) {
                    for (k = 0; k < G.cols[i].deg; k++)  {
                        if (G.cols[i].adj[k]->del == 0) {
                            cnt++;
                        }
                    }
                }
            }
            printf ("Current NZCOUNT = %d\n", cnt); fflush (stdout);
        }
*/
    } while (hit);

    if (ILL_PRE_EMPTY_COL) {
        rval = empty_columns (&G, pre);
        ILL_CLEANUP_IF (rval);
    }

    if (debug) {
        printf ("Operations\n");
        for (i = 0; i < pre->opcount; i++) {
            switch (pre->oplist[i].ptype) {
            case ILL_PRE_DELETE_EMPTY_ROW:
                printf ("Delete Empty Row: %d\n", pre->oplist[i].rowindex);
                fflush (stdout);
                break;
            case ILL_PRE_DELETE_SINGLETON_ROW:
                printf ("Delete Singleton Row: %d (col %d)\n",
                    pre->oplist[i].rowindex, pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            case ILL_PRE_DELETE_FIXED_VARIABLE:
                printf ("Delete Fixed Variable: %d\n",
                                                pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            case ILL_PRE_DELETE_FORCED_VARIABLE:
                printf ("Delete Forced Variable: %d\n",
                                                pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            case ILL_PRE_DELETE_SINGLETON_VARIABLE:
                printf ("Delete Singleton Variable: %d\n",
                                                pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            case ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE:
                printf ("Delete Free Singleton Variable: %d\n",
                                                pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            case ILL_PRE_DELETE_EMPTY_COLUMN:
                printf ("Delete Empty Column: %d\n",
                                                pre->oplist[i].colindex);
                fflush (stdout);
                dump_line (&pre->oplist[i].line);
                break;
            default:
                fprintf (stderr, "unknon presolve operation\n");
                rval = 1; goto CLEANUP;
            }
        }
        printf ("\n");
    }

    rval = grab_lp_info (&G, lp->colnames, info);
    ILL_CLEANUP_IF (rval);

/*
    printf ("Final Rows = %d, Cols = %d, Nzcount = %d\n",
               info->nrows, info->ncols, info->nzcount);
    fflush (stdout);
*/


CLEANUP:

    free_graph (&G);
    ILL_RETURN (rval, "simple_presolve");
}

static int grab_lp_line (graph *G, int indx, ILLlp_preline *line,
        int row_or_col)
{
    int rval = 0;
    int k, cnt;
    node *n;

    if (row_or_col == 0) n = &G->rows[indx];
    else                 n = &G->cols[indx];

    line->count = 0;

    for (k = 0; k < n->deg; k++) {
        if (n->adj[k]->del == 0) {
            line->count++;
        }
    }

    if (line->count) {
        ILL_SAFE_MALLOC (line->ind, line->count, int);
        ILL_SAFE_MALLOC (line->val, line->count, double);
        if (!line->ind || !line->val) {
            fprintf (stderr, "out of memory in grab_lp_line\n");
            rval = 1; goto CLEANUP;
        }
        for (k = 0, cnt = 0; k < n->deg; k++) {
            if (n->adj[k]->del == 0) {
                line->ind[cnt] = n->adj[k]->row;
                line->val[cnt] = n->adj[k]->coef;
                cnt++;
            }
        }
    }

    if (row_or_col == 0) {
        line->rhs   = n->rhs;
    } else {
        line->obj   = n->obj;
        line->lower = n->lower;
        line->upper = n->upper;
    }

    line->row_or_col = row_or_col;

CLEANUP:

    ILL_RETURN (rval, "grab_lp_line");
}

static void dump_line (ILLlp_preline *line)
{
    int k;

    printf (" ");
    if (line->row_or_col == 0) {
        for (k = 0; k < line->count; k++) {
            printf (" C%d->%g", line->ind[k], line->val[k]);
        }
        printf (" RHS->%g\n", line->rhs);
    } else {
        for (k = 0; k < line->count; k++) {
            printf (" R%d->%g", line->ind[k], line->val[k]);
        }
        printf (" Obj->%g  LB->%g  UB->%g\n", line->obj, line->lower,
                                                         line->upper);
    }
    fflush (stdout);
}

static int grab_lp_info (graph *G, char **colnames, ILLlp_sinfo *info)
{
    int rval = 0;
    int ncols = 0, nrows = 0, nzcount = 0;
    int i, j, k, cnt, len;
    node *grows = G->rows;
    node *gcols = G->cols;
    int *tdeg = (int *) NULL;
    int *map  = (int *) NULL;
    char *buf = (char *) NULL;
    ILLmatrix *A = &info->A;

    ILL_SAFE_MALLOC (tdeg, G->ncols, int);
    ILL_SAFE_MALLOC (map, G->nrows, int);
    if (!tdeg || !map) {
        fprintf (stderr, "out of memory in grab_lp_info\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < G->nrows; i++) {
        if (grows[i].del == 0) {
            map[i] = nrows;
            nrows++;
        }
    }

    for (j = 0; j < G->ncols; j++) {
        if (gcols[j].del == 0) {
            tdeg[ncols] = 0;
            for (k = 0; k < gcols[j].deg; k++) {
                if (gcols[j].adj[k]->del == 0) {
                    tdeg[ncols]++;
                    nzcount++;
                }
            }
            ncols++;
        }
    }

    info->ncols   = ncols;
    info->nrows   = nrows;
    info->nzcount = nzcount;

    info->rowsize = nrows; 
    info->colsize = ncols;

    ILL_SAFE_MALLOC (info->rhs, info->rowsize, double);
    ILL_SAFE_MALLOC (info->obj, info->colsize, double);
    ILL_SAFE_MALLOC (info->lower, info->colsize, double);
    ILL_SAFE_MALLOC (info->upper, info->colsize, double);
    ILL_SAFE_MALLOC (A->matval, info->nzcount + 1, double);
    ILL_SAFE_MALLOC (A->matind, info->nzcount + 1, int);
    ILL_SAFE_MALLOC (A->matcnt, info->colsize, int);
    ILL_SAFE_MALLOC (A->matbeg, info->colsize, int);

    if (!info->rhs || !info->obj || !info->lower || !info->upper ||
        !A->matval || !A->matind || !A->matcnt   || !A->matbeg) { 
        fprintf (stderr, "out of memory in grab_lp\n");
        rval = 1; goto CLEANUP;
    }

    A->matind[info->nzcount] = -1;
    A->matsize = info->nzcount + 1;
    A->matcolsize = info->colsize;
    A->matfree = 1;
    A->matcols = ncols;
    A->matrows = nrows;


    nrows = 0;
    for (i = 0; i < G->nrows; i++) {
        if (grows[i].del == 0) {
            info->rhs[nrows] = grows[i].rhs;
            nrows++;
        }
    }

    ncols = 0;
    cnt   = 0;
    for (j = 0; j < G->ncols; j++) {
        if (gcols[j].del == 0) {
            info->obj[ncols]   = gcols[j].obj;
            info->lower[ncols] = gcols[j].lower;
            info->upper[ncols] = gcols[j].upper;
            A->matcnt[ncols] = tdeg[ncols];
            A->matbeg[ncols] = cnt;
            for (k = 0; k < gcols[j].deg; k++) {
                if (gcols[j].adj[k]->del == 0) {
                    A->matval[cnt] = gcols[j].adj[k]->coef;
                    A->matind[cnt] = map[gcols[j].adj[k]->row];
                    cnt++;
                }
            }
            ncols++;
        }
    }

    if (colnames) {
        ILL_SAFE_MALLOC (info->colnames, info->colsize, char *);
        if (!info->colnames) {
            fprintf (stderr, "out of memory in grab_lp\n");
            rval = 1; goto CLEANUP;
        }
        for (j = 0; j < info->colsize; j++) {
            info->colnames[j] = (char *) NULL;
        }

        ILL_SAFE_MALLOC (buf, ILL_namebufsize, char);
        if (!buf) {
            fprintf (stderr, "out of memory in grab_lp\n");
            rval = 1; goto CLEANUP;
        }
        ncols = 0;
        for (j = 0; j < G->ncols; j++) {
            if (gcols[j].del == 0) {
                if (gcols[j].coltype == ILL_PRE_COL_STRUC) {
                    len = strlen (colnames[j]) + 1;
                    ILL_SAFE_MALLOC (info->colnames[ncols], len, char);
                    if (!info->colnames[ncols]) {
                        fprintf (stderr, "out of memory in grab_lp\n");
                        rval = 1; goto CLEANUP;
                    }
                    strcpy (info->colnames[ncols], colnames[j]);
                } else {
                    for (k = 0; k < gcols[j].deg; k++) {
                        if (gcols[j].adj[k]->del == 0) {
                            i = gcols[j].adj[k]->row;
                            break;
                        }
                    }
                    if (k == gcols[j].deg) {
                        fprintf (stderr, "problem with graph in grab_lp\n");
                        rval = 1; goto CLEANUP;
                    }
                    sprintf (buf, "s%d", i);
                    len = strlen (buf) + 1;
                    ILL_SAFE_MALLOC (info->colnames[ncols], len, char);
                    if (!info->colnames[ncols]) {
                        fprintf (stderr, "out of memory in grab_lp\n");
                        rval = 1; goto CLEANUP;
                    }
                    strcpy (info->colnames[ncols], buf);
                }
                ncols++;
            }
        }
    }

/* ADD STRUCT VARIABLE STUFF */
    

CLEANUP:

    if (rval) {
        ILLlp_sinfo_free (info);
    }
    ILL_IFFREE (tdeg, int);
    ILL_IFFREE (map,  int);
    ILL_IFFREE (buf, char);

    ILL_RETURN (rval, "grab_lp_info");
}

static int fixed_variables (graph *G, ILLlp_predata *pre)
{
    int rval = 0;
    int j;
    int ncols  = G->ncols;
    node *cols = G->cols;
    ILLlp_preop *op;

    for (j = 0; j < ncols; j++) {
        if (cols[j].del == 0) {
            if (cols[j].lower == cols[j].upper) {
                rval = get_next_preop (pre, &op);
                ILL_CLEANUP_IF (rval);

                op->colindex = j;
                op->rowindex = -1;
                op->ptype = ILL_PRE_DELETE_FIXED_VARIABLE;

                rval = grab_lp_line (G, op->colindex, &op->line, 1);
                ILL_CLEANUP_IF (rval);
                pre->opcount++;

                set_fixed_variable (G, j, cols[j].lower);
            }
        }
    }

CLEANUP:

    ILL_RETURN (rval, "fixed_variables");
}

static int empty_columns (graph *G, ILLlp_predata *pre)
{
    int rval = 0;
    int j, k;
    int ncols  = G->ncols;
    node *cols = G->cols;
    ILLlp_preop *op;

    for (j = 0; j < ncols; j++) {
        if (cols[j].del == 0) {
            for (k = 0; k < cols[j].deg; k++) {
                if (cols[j].adj[k]->del == 0) break;
            }
            if (k == cols[j].deg) {
                rval = get_next_preop (pre, &op);
                ILL_CLEANUP_IF (rval);

                op->colindex = j;
                op->rowindex = -1;
                op->ptype = ILL_PRE_DELETE_EMPTY_COLUMN;

                rval = grab_lp_line (G, op->colindex, &op->line, 1);
                ILL_CLEANUP_IF (rval);
                pre->opcount++;

                if (cols[j].obj * G->objsense > ILL_PRE_FEAS_TOL) {
                    if (cols[j].lower == -ILL_MAXDOUBLE) {
                        printf ("unbounded prob detected in empty_columns\n");
                        printf ("col %d, obj %g\n", j, cols[j].obj);
                        fflush (stdout);
                        rval = 1; goto CLEANUP;
                    } else {
                        set_fixed_variable (G, j, cols[j].lower);
                    }
                } else if (cols[j].obj * G->objsense < -ILL_PRE_FEAS_TOL) {
                    if (cols[j].upper == ILL_MAXDOUBLE) {
                        printf ("unbounded prob detected in empty_columns\n");
                        printf ("col %d, obj %g\n", j, cols[j].obj);
                        fflush (stdout);
                        rval = 1; goto CLEANUP;
                    } else {
                        set_fixed_variable (G, j, cols[j].upper);
                    }
                } else {
                    set_fixed_variable (G, j, cols[j].lower);
                }
            }
        }
    }

CLEANUP:

    ILL_RETURN (rval, "empty_columns");
}

static int singleton_rows (graph *G, ILLlp_predata *pre, int *hit)
{
    int rval = 0;
    int rowindex, i, k, h;
    int nrows  = G->nrows;
    node *rows = G->rows;
    node *cols = G->cols;
    node *r, *c;
    edge *pivot, *f;
    intptr *next, *list = (intptr *) NULL;
    int *tdeg =  (int *) NULL;
    double val;
    ILLlp_preop *op;

    *hit = 0;
    if (G->nrows == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (tdeg, G->nrows, int);
    if (!tdeg) {
        fprintf (stderr, "out of memory in singleton_rows\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < nrows; i++) {
        if (rows[i].del == 0) {
            tdeg[i] = 0;
            for (k = 0; k < rows[i].deg; k++) {
                if (rows[i].adj[k]->del == 0) {
                    tdeg[i]++;
                }
            }
            if (tdeg[i] <= 1) {
                rval = add_to_list (&G->intptrworld, &list, i);
                ILL_CLEANUP_IF (rval);
            }
        }
    } 

    while (list) {
        (*hit)++;
        rowindex = list->this;
        next     = list->next;
        intptrfree (&G->intptrworld, list);
        list = next;

        rval = get_next_preop (pre, &op);
        ILL_CLEANUP_IF (rval);

        r = &rows[rowindex];

        if (tdeg[rowindex] == 0) {
            if (r->rhs < -ILL_PRE_FEAS_TOL || r->rhs > ILL_PRE_FEAS_TOL) {
                printf ("infeasible row detected in singleton_row\n");
                printf ("empty row with rhs = %g\n", r->rhs);
                fflush (stdout);
                rval = 1; goto CLEANUP;
            }
            op->ptype = ILL_PRE_DELETE_EMPTY_ROW;
            op->rowindex = rowindex;
        } else {
            /*  Find the "pivot" entry and colum */

            for (k = 0; k < r->deg; k++) {
                if (r->adj[k]->del == 0) break;
            }
            if (k == r->deg) {
                fprintf (stderr, "lost an edge in singleton_rows\n");
                rval = 1; goto CLEANUP;
            }

            pivot = r->adj[k];
            c = &cols[pivot->col];

            /*  Store data from operation (incluing the col coefs) */

            op->ptype    = ILL_PRE_DELETE_SINGLETON_ROW;
            op->rowindex = rowindex;
            op->colindex = c - cols; 
            op->line.rhs = r->rhs;
            rval = grab_lp_line (G, op->colindex, &op->line, 1);
            ILL_CLEANUP_IF (rval);

            /*  Fix the x[c] to its rhs value */

            val = r->rhs / pivot->coef;
            if (val < c->lower - ILL_PRE_FEAS_TOL ||
                val > c->upper + ILL_PRE_FEAS_TOL) {
                printf ("infeasible bounds detected in singleton_row %d\n",
                                 rowindex);
                printf ("lower->%g  upper->%g  val = %g\n",
                                 c->lower, c->upper, val);
                fflush (stdout);
                rval = 1; goto CLEANUP;
            } else {
                c->lower = val;
                c->upper = val;
            }
 
            /*  Delete x[c] from other rows (and adjust their rhs) */ 
 
            c->del = 1;

            for (h = 0; h < c->deg; h++) {
                f = c->adj[h];
                if (f->del == 0) {
                    rows[f->row].rhs -= (f->coef * c->lower);
                    tdeg[f->row]--;
                    if (tdeg[f->row] == 1) {
                        if (f == pivot) {
                            fprintf (stderr, "bad pivot element\n");
                            rval = 1; goto CLEANUP;
                        }
                        rval = add_to_list (&G->intptrworld, &list, f->row);
                        ILL_CLEANUP_IF (rval);
                    }
                    f->del = 1;
                }
            }
        }

        r->del = 1;
        pre->opcount++;
    }

CLEANUP:

    ILL_IFFREE (tdeg, int);
    intptr_listfree (&G->intptrworld, list);

    ILL_RETURN (rval, "singleton_rows");
}

static int forcing_constraints (graph *G, ILLlp_predata *pre, int *hit)
{
    int rval = 0;
    int i, j, k, ts;
    node *rows = G->rows;
    node *cols = G->cols;
    edge *e;
    int nrows  = G->nrows;
    double ub, lb;
    ILLlp_preop *op;

    *hit = 0;

    for (i = 0; i < nrows; i++) {
        if (rows[i].del == 0) {
            get_implied_rhs_bounds (G, i, &lb, &ub);
            if (lb > rows[i].rhs + ILL_PRE_FEAS_TOL ||
                ub < rows[i].rhs - ILL_PRE_FEAS_TOL) {
                printf ("infeasible row detected in forcing_constraints\n");
                printf ("Row %d:  RHS->%g  LBnd->%g  UBnd->%g\n",
                            i, rows[i].rhs, lb, ub);
                fflush (stdout);
                rval = 1; goto CLEANUP;
            } else if (lb > rows[i].rhs - ILL_PRE_FEAS_TOL ||
                       ub < rows[i].rhs + ILL_PRE_FEAS_TOL) {
                (*hit)++;
                ts = (lb > rows[i].rhs - ILL_PRE_FEAS_TOL ? 0 : 1);
                for (k = 0; k < rows[i].deg; k++) {
                    e = rows[i].adj[k];
                    if (e->del == 0) {
                        j = e->col;

                        rval = get_next_preop (pre, &op);
                        ILL_CLEANUP_IF (rval);

                        op->colindex = j;
                        op->rowindex = i;
                        op->ptype = ILL_PRE_DELETE_FORCED_VARIABLE;

                        rval = grab_lp_line (G, j, &op->line, 1);
                        ILL_CLEANUP_IF (rval);
                        pre->opcount++;

                        if ((ts == 0 && e->coef < 0.0) ||
                            (ts == 1 && e->coef > 0.0)) {
                            set_fixed_variable (G, j, cols[j].upper);
                        } else {
                            set_fixed_variable (G, j, cols[j].lower);
                        }
                    }
                }
            }
        }
    }

CLEANUP:

    ILL_RETURN (rval, "forcing_constraints");
}

static int singleton_columns (graph *G, ILLlp_predata *pre, int *hit)
{
    int rval = 0;
    int ncols = G->ncols;
    int j, k, deg, rdeg, single = 0, irow;
    double lb, ub, b;
    node *cols = G->cols;
    node *rows = G->rows;
    edge *b_edge;
    ILLlp_preop *op;

    *hit = 0;
    if (G->ncols == 0) goto CLEANUP;

    for (j = 0; j < ncols; j++) {
        if (cols[j].del == 0) {
            deg = 0;
            for (k = 0; k < cols[j].deg && deg <= 1; k++) {
                if (cols[j].adj[k]->del == 0) {
                    single = k;
                    deg++;
                }
            }
            if (deg == 1) {
                irow   = cols[j].adj[single]->row;
                b      = cols[j].adj[single]->coef;
                b_edge = cols[j].adj[single];

                get_implied_variable_bounds (G, j, b_edge, &lb, &ub);

                if (lb >= cols[j].lower && ub <= cols[j].upper) {
                    double eb;
                    edge  *a_edge;

                    /*  The jth variable can be substituted out of problem */
                    /*        x = (c/b) - (a/b)y                           */


                    rval = get_next_preop (pre, &op);
                    ILL_CLEANUP_IF (rval);

                    op->colindex = j;
                    op->rowindex = irow;
                    op->ptype = ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE;

                    rval = grab_lp_line (G, irow, &op->line, 0);
                    ILL_CLEANUP_IF (rval);
                    pre->opcount++;

                    /*  Adjust the objective function                      */
                    /*     dy ==> (d - (e/b))ay   (e is obj coef of y)     */

                    eb = cols[j].obj / b;

                    for (k = 0; k < rows[irow].deg; k++) {
                        a_edge = rows[irow].adj[k];
                        if (a_edge->del == 0 && a_edge != b_edge) {
                            cols[a_edge->col].obj -= (eb * a_edge->coef);
                        }
                    }


                    /*  Delete y from graph */

                    cols[j].del = 1;

                    /*  Delete equation ay + bx = c */

                    rows[irow].del = 1;
                    for (k = 0; k < rows[irow].deg; k++) {
                        rows[irow].adj[k]->del = 1;
                    }

                } else {
                    rdeg = 0;
                    for (k = 0; k < rows[irow].deg && rdeg <= 2; k++) {
                        if (rows[irow].adj[k]->del == 0) {
                            rdeg++;
                        }
                    }
                    if (rdeg == 2) {
                        double newub = ILL_MAXDOUBLE, newlb = -ILL_MAXDOUBLE;
                        double a = 0.0, c, l, u;
                        edge *a_edge = (edge *) NULL;
                        int col2 = 0;

                        /*    ay + bx = c                                */
                        /*    l <= x <= u                                */
                        /*      x - is column singleton                  */
                        /*      derive bounds on y and substitute out x  */

                        c      = rows[irow].rhs;
                        l      = cols[j].lower;
                        u      = cols[j].upper;

                        /* Find the ay term */

                        for (k = 0; k < rows[irow].deg; k++) {
                            if (rows[irow].adj[k]->del == 0 &&
                                rows[irow].adj[k]->col != j) {
                                a_edge = rows[irow].adj[k];
                                a      = rows[irow].adj[k]->coef;
                                col2   = rows[irow].adj[k]->col;
                                break;
                            }
                        }
                        if (k == rows[irow].deg) {
                            fprintf (stderr, "graph error in singleton_col\n");
                            rval = 1; goto CLEANUP;
                        }

                        /*  Record the operation             */
                        /*  x is column j,  y is column col2 */

                        rval = get_next_preop (pre, &op);
                        ILL_CLEANUP_IF (rval);

                        op->colindex = j;
                        op->rowindex = irow;
                        op->ptype = ILL_PRE_DELETE_SINGLETON_VARIABLE;

                        rval = grab_lp_line (G, irow, &op->line, 0);
                        ILL_CLEANUP_IF (rval);
                        pre->opcount++;

                        /*  Adjust the bounds on y           */
                        /*  Using x = c/b - (a/b)y            */

                        if (a/b > 0) {
                            if (l > -ILL_MAXDOUBLE) newub = (c/a) - (l*b)/a;
                            if (u <  ILL_MAXDOUBLE) newlb = (c/a) - (u*b)/a;
                        } else {
                            if (l > -ILL_MAXDOUBLE) newlb = (c/a) - (l*b)/a;
                            if (u <  ILL_MAXDOUBLE) newub = (c/a) - (u*b)/a;
                        }

                        if (newlb > cols[col2].lower) cols[col2].lower = newlb;
                        if (newub < cols[col2].upper) cols[col2].upper = newub;
                        cols[col2].obj -= a/b;

                        /*  Delete x (and the bx term) from graph */

                        cols[j].del = 1;
                        b_edge->del = 1;

                        /*  Delete equation ay + bx = c (and the ax term) */

                        rows[irow].del = 1;
                        a_edge->del = 1;
                    }
                }
            }
        }
    }


CLEANUP:

    ILL_RETURN (rval, "singleton_columns");
}

static int duplicate_rows (graph *G, int *hit)
{
    int rval = 0;
    node *cols = G->cols;
    node *rows = G->rows;
    int  ncols = G->ncols;
    int  nrows = G->nrows;
    int     *s = (int *) NULL;
    double  *f = (double *) NULL;
    double szeit = ILLutil_zeit ();
    double q;
    int i, j, k, k2, ri, r0 = 0, n, nu = 0, got, t0, t = 1;
    node *c;


    /*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (nrows == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (s, nrows, int);

    ILL_SAFE_MALLOC (f, nrows, double);

    for (i = 0; i < nrows; i++) {
        if (rows[i].del || rows[i].rowsense != 'E') {
            s[i] = ILL_MAXINT;   /* ILL_MAXINT means no longer eligible    */
        } else {
            s[i] = 0;           /* 0 means eligible, >0 means in a group */ 
            nu++;               /* Tracks the number of eligible rows    */
        }
    }

    for (j = 0; j < ncols; j++) {
        c = &cols[j];
        if (c->del) continue;
        if (c->coltype != ILL_PRE_COL_STRUC) continue;

        n = 0;
        t0 = t++;

        for (k = 0; k < c->deg; k++) {
            if (c->adj[k]->del) continue;
            
            ri = c->adj[k]->row;
            if (s[ri] == 0) {
                s[ri] = t0;
                f[ri] = c->adj[k]->coef;
                r0    = ri;
                n++;
            } else if (s[ri] < t0) {
                got = 0;
                for (k2 = k+1; k2 < c->deg; k2++) {
                    if (c->adj[k2]->del) continue;

                    i = c->adj[k2]->row;
                    if (s[i] == s[ri]) {
                        q = (c->adj[k]->coef *(f[i])) /
                            (f[ri] * (c->adj[k2]->coef));
                        if (q >= 1 - ILL_PRE_ZERO_TOL &&
                            q <= 1 + ILL_PRE_ZERO_TOL) {
                            s[ri] = t;
                            s[i]  = t;
                            got++;
                        }
                    }
                }
                if (got) {
                    t++;
                } else {
                    s[ri] = ILL_MAXINT;
                    if (--nu == 0) goto DONE;
                }
            }
        }

        if (n == 1) {
            s[r0] = ILL_MAXINT;
            if (--nu == 0) goto DONE;
        }
    }

DONE:

    {
        int idup = 0;

        for (i = 0; i < nrows; i++)  {
            if (s[i] > 0 && s[i] < ILL_MAXINT) {
                printf ("Row %d: %d\n", i, s[i]);
                idup++;
            }
        }
        printf ("Number of duplicate rows: %d\n", idup);
    }

    printf ("Time in duplicate_rows: %.2f (seconds)\n", ILLutil_zeit() - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    ILL_IFFREE (f, double);

    ILL_RETURN (rval, "duplicate_rows");
}

static int duplicate_cols (graph *G, int *hit)
{
    int rval = 0;
    node *cols = G->cols;
    node *rows = G->rows;
    int  ncols = G->ncols;
    int  nrows = G->nrows;
    int     *s = (int *) NULL;
    double  *f = (double *) NULL;
    double szeit = ILLutil_zeit ();
    double q;
    int i, j, k, k2, ci, c0 = 0, n, nu = 0, got, t0, t = 1;
    node *r;


    /*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (ncols == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (s, ncols, int);

    ILL_SAFE_MALLOC (f, ncols, double);

    for (j = 0; j < ncols; j++) {
        if (cols[j].del || cols[j].coltype != ILL_PRE_COL_STRUC) {
            s[j] = ILL_MAXINT;   /* ILL_MAXINT means no longer eligible    */
        } else {
            s[j] = 0;           /* 0 means eligible, >0 means in a group */ 
            nu++;               /* Tracks the number of eligible rows    */
        }
    }

    for (i = 0; i < nrows; i++) {
        r = &rows[i];
        if (r->del) continue;

        n = 0;
        t0 = t++;

        for (k = 0; k < r->deg; k++) {
            if (r->adj[k]->del) continue;
            
            ci = r->adj[k]->col;
            if (s[ci] == 0) {
                s[ci] = t0;
                f[ci] = r->adj[k]->coef;
                c0    = ci;
                n++;
            } else if (s[ci] < t0) {
                got = 0;
                for (k2 = k+1; k2 < r->deg; k2++) {
                    if (r->adj[k2]->del) continue;

                    j = r->adj[k2]->col;
                    if (s[j] == s[ci]) {
                        q = (r->adj[k]->coef *(f[j])) /
                            (f[ci] * (r->adj[k2]->coef));
                        if (q >= 1 - ILL_PRE_ZERO_TOL &&
                            q <= 1 + ILL_PRE_ZERO_TOL) {
                            s[ci] = t;
                            s[j]  = t;
                            got++;
                        }
                    }
                }
                if (got) {
                    t++;
                } else {
                    s[ci] = ILL_MAXINT;
                    if (--nu == 0) goto DONE;
                }
            }
        }

        if (n == 1) {
            s[c0] = ILL_MAXINT;
            if (--nu == 0) goto DONE;
        }
    }

DONE:

    {
        int dcount;
        int *dcnt;
        int *dlist;

        rval = gather_dup_lists ( s, ncols, &dcount, &dcnt, &dlist);
        ILL_CLEANUP_IF (rval);
    }

    printf ("Time in duplicate_cols: %.2f (seconds)\n", ILLutil_zeit() - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    ILL_IFFREE (f, double);

    ILL_RETURN (rval, "duplicate_cols");
}

static int gather_dup_lists (/* graph *G, */ int *s, /* double *f, */ int count,
    int *duptotal, int **dupcnt, int **dupind)
{
    int rval = 0;
    int *cnt   = (int *) NULL;
    int *ind   = (int *) NULL;
    int *beg   = (int *) NULL;
    int i, smax = 0, ndup = 0, total = 0;

    *duptotal = 0;
    *dupcnt = (int *) NULL;
    *dupind = (int *) NULL;


    for (i = 0; i < count; i++) {
        if (s[i] < ILL_MAXINT && s[i] > smax) smax = s[i];
    }
    if (smax == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (cnt, smax+1, int);

    ILL_SAFE_MALLOC (ind, smax+1, int);

    for (i = 0; i < smax+1; i++) {
        cnt[i] = 0;
    }

    for (i = 0; i < count; i++) {
        if (s[i] < ILL_MAXINT) {
            cnt[s[i]]++;
        }
    }

    if (cnt[0] > 0) printf ("%d Empty Lines\n", cnt[0]);

    printf ("Duplicate Classes:"); fflush (stdout);
    for (i = 1; i < smax+1; i++) {
        if (cnt[i] > 1) {
            ndup++;
            printf (" %d", cnt[i]);
        }
    }
    printf ("  Number %d\n", ndup);
    fflush (stdout);

    if (ndup == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (beg, ndup, int);

    for (i = 1, ndup = 0; i < smax+1; i++) {
        if (cnt[i] > 1) {
            beg[ndup] = total;
            total += cnt[i];
            ind[i] = ndup;
            ndup++;
        }
    }

    if (total == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (*dupcnt, ndup, int);

    ILL_SAFE_MALLOC (*dupind, total, int);

    for (i = 0; i < ndup; i++) {
        (*dupcnt)[i] = 0;
    }

    for (i = 0; i < count; i++) {
        if (s[i] < ILL_MAXINT && s[i] > 0) {
            if (cnt[s[i]] > 1) {
                (*dupind)[beg[ind[s[i]]] + (*dupcnt)[ind[s[i]]]] = i;
                (*dupcnt)[ind[s[i]]]++;
            }
        }
    }

    for (i = 0;  i < ndup; i++) {
        int j;
        for (j = beg[i]; j < beg[i] + (*dupcnt)[i]; j++) {
            printf (" %d", (*dupind)[j]);
        }
        printf (" | "); 
        fflush (stdout);
    }

    *duptotal = ndup;

CLEANUP:

    ILL_IFFREE (cnt, int);
    ILL_IFFREE (ind, int);
    ILL_IFFREE (beg, int);

    ILL_RETURN (rval, "gather_dup_lists");
}

static void set_fixed_variable (graph *G, int j, double val)
{
    int k;
    edge *e;

    G->cols[j].del = 1;
    for (k = 0; k < G->cols[j].deg; k++) {
        e = G->cols[j].adj[k];
        if (e->del == 0) {
            G->rows[e->row].rhs -= (e->coef * val);
            e->del = 1;
        }
    }
}

static void get_implied_rhs_bounds (graph *G, int i, double *lb, double *ub)
{
    int k;
    double l, u;
    node *cols = G->cols;
    node *rows = G->rows;
    edge *e;

    l = 0.0;
    for (k = 0; k < rows[i].deg; k++) {
        e = rows[i].adj[k];
        if (e->del == 0) {
            if (e->coef < 0.0) {
                if (cols[e->col].upper == ILL_MAXDOUBLE) {
                    l = -ILL_MAXDOUBLE;
                    break;
                } else {
                    l += (e->coef * cols[e->col].upper);
                }
            } else if (e->coef > 0.0) {
                if (cols[e->col].lower == -ILL_MAXDOUBLE) {
                    l = -ILL_MAXDOUBLE;
                    break;
                } else {
                    l += (e->coef * cols[e->col].lower);
                }
            }
        }
    }
    u = 0.0;
    for (k = 0; k < rows[i].deg; k++) {
        e = rows[i].adj[k];
        if (e->del == 0) {
            if (e->coef < 0.0) {
                if (cols[e->col].lower == -ILL_MAXDOUBLE) {
                    u = ILL_MAXDOUBLE;
                } else {
                    u += (e->coef * cols[e->col].lower);
                }
            } else if (e->coef > 0.0) {
                if (cols[e->col].upper == ILL_MAXDOUBLE) {
                    u = ILL_MAXDOUBLE;
                } else {
                    u += (e->coef * cols[e->col].upper);
                }
            }
        }
    }

    *lb = l;
    *ub = u;
}

static void get_implied_variable_bounds (graph *G, int j, edge *a_ij,
        double *lb, double *ub)
{
    int i = a_ij->row;
    double l, u;

    get_implied_rhs_bounds (G, i, &l, &u);

    *lb = -ILL_MAXDOUBLE;
    *ub =  ILL_MAXDOUBLE;

    if (a_ij->coef > ILL_PRE_FEAS_TOL) {
       if (u <  ILL_MAXDOUBLE) {
           *lb = (G->rows[i].rhs - u)/a_ij->coef + G->cols[j].upper;
       }
       if (l > -ILL_MAXDOUBLE) {
           *ub = (G->rows[i].rhs - l)/a_ij->coef + G->cols[j].lower;
       }
    } else if (a_ij->coef < ILL_PRE_FEAS_TOL) {
       if (l > -ILL_MAXDOUBLE) {
           *lb = (G->rows[i].rhs - l)/a_ij->coef + G->cols[j].upper;
       }
       if (u <  ILL_MAXDOUBLE) {
           *ub = (G->rows[i].rhs - u)/a_ij->coef + G->cols[j].lower;
       }
    }
}

static int get_next_preop (ILLlp_predata *pre, ILLlp_preop **op)
{
    int rval = 0;

    if (pre->opcount >= pre->opsize) {
        rval = ILLutil_reallocrus_scale ((void **) &pre->oplist,
                     &pre->opsize, pre->opcount + 1, 1.3,
                      sizeof (ILLlp_preop));
        ILL_CLEANUP_IF (rval);
    }
    *op = &pre->oplist[pre->opcount];
    ILLlp_preop_init (*op);

CLEANUP:

    ILL_RETURN (rval, "get_next_preop");
}

static int add_to_list (ILLptrworld *world, intptr **list, int i)
{
    int rval = 0;
    intptr *ip;

    ip = intptralloc (world);
    if (!ip) { rval = 1; goto CLEANUP; }
    ip->this = i;
    ip->next = *list;
    *list = ip;

CLEANUP:

    ILL_RETURN (rval, "add_to_list");
}

static int build_graph (ILLlpdata *lp, graph *G)
{
    int rval = 0;
    int ncols   = lp->ncols;
    int nrows   = lp->nrows;
    int nzcount = lp->nzcount;
    int i, j, k, stop, count;
    edge *edgelist;
    node *rows, *cols;
    ILLmatrix *A = &lp->A;

    G->objsense = lp->objsense;

    ILL_SAFE_MALLOC (G->rows, nrows, node);
    if (!G->rows) {
        fprintf (stderr, "out of memory in build_graph\n");
        rval = 1; goto CLEANUP;
    }
    rows = G->rows;

    for (i = 0; i < nrows; i++) {
        rows[i].rowsense = lp->sense[i];
        rows[i].deg = 0;
    }

    ILL_SAFE_MALLOC (G->cols, ncols, node);
    ILL_SAFE_MALLOC (G->edgelist, nzcount, edge);
    ILL_SAFE_MALLOC (G->adjspace, 2 * nzcount, edge *);

    if (!G->cols || !G->edgelist || !G->adjspace) {
        fprintf (stderr, "out of memory in build_graph\n");
        rval = 1; goto CLEANUP;
    }

    cols     = G->cols;
    edgelist = G->edgelist;

    for (j = 0; j < ncols; j++) {
        stop = A->matbeg[j] + A->matcnt[j];
        for (k = A->matbeg[j]; k < stop; k++) {
            rows[A->matind[k]].deg++;
        }
    }

    for (i = 0, count = 0; i < nrows; i++) {
        rows[i].adj = G->adjspace + count;
        count += rows[i].deg;
        rows[i].deg = 0;
    }

    for (j = 0; j < ncols; j++) {
        cols[j].adj = G->adjspace + count;
        count += A->matcnt[j];
        cols[j].deg = 0;
        cols[j].coltype = ILL_PRE_COL_STRUC;
    }
    for (i = 0; i < nrows; i++) {
        cols[lp->rowmap[i]].coltype = ILL_PRE_COL_LOGICAL;
    }

    for (j = 0, count = 0; j < ncols; j++) {
        cols[j].obj   = lp->obj[j];
        cols[j].lower = lp->lower[j];
        cols[j].upper = lp->upper[j];
        stop = A->matbeg[j] + A->matcnt[j];
        for (k = A->matbeg[j]; k < stop; k++) {
            i = A->matind[k];
            rows[i].adj[rows[i].deg++] = &(edgelist[count]); 
            cols[j].adj[cols[j].deg++] = &(edgelist[count]);
            edgelist[count].row     = i;
            edgelist[count].col     = j;
            edgelist[count].coef    = A->matval[k];
            edgelist[count].mark    = 0;
            edgelist[count].del     = 0;
            edgelist[count].coltype = cols[j].coltype;
            count++;
        }
    }
    if (count != nzcount) {
        fprintf (stderr, "counts are off in build_graph\n");
        rval = 1; goto CLEANUP;
    }

    G->ecount = count;
    G->nrows  = nrows;
    G->ncols = ncols;

    for (i = 0; i < G->nrows; i++) {
        G->rows[i].del = 0;
        G->rows[i].rhs = lp->rhs[i];
    }
    for (j = 0; j < G->ncols; j++) {
        G->cols[j].del = 0;
    }

CLEANUP:

    ILL_RETURN (rval, "build_graph");
}

static void dump_graph (graph *G)
{
    int i, j, k;

    printf ("ecount = %d, nrows = %d, ncols = %d\n",
               G->ecount, G->nrows, G->ncols);
    fflush (stdout);

    for (i = 0; i < G->nrows; i++) {
        printf ("Row %d:", i);
        for (k = 0; k < G->rows[i].deg; k++) {
            printf (" %d", G->rows[i].adj[k]->col);
            if (G->rows[i].adj[k]->coltype == ILL_PRE_COL_LOGICAL) printf ("S");
            printf ("(%g)", G->rows[i].adj[k]->coef);
        }
        printf ("  rhs: %g", G->rows[i].rhs);
        if (G->rows[i].del) {
            printf (" (deleted)\n");
        } else {
            printf ("\n");
        }
    }

    for (j = 0; j < G->ncols; j++) {
        if (G->cols[j].coltype == ILL_PRE_COL_LOGICAL) {
            printf ("Slk %d:", j);
        } else {
            printf ("Col %d:", j);
        }
        for (k = 0; k < G->cols[j].deg; k++) {
            printf (" %d", G->cols[j].adj[k]->row);
        }
        printf ("  obj: %g  bnd: (%g, %g)", G->cols[j].obj, G->cols[j].lower,
                                                 G->cols[j].upper);
        if (G->cols[j].del) {
            printf (" (deleted)\n");
        } else {
            printf ("\n");
        }
    }
}

static void init_graph (graph *G)
{
    if (G) {
        G->edgelist = (edge *) NULL;
        G->rows     = (node *) NULL;
        G->cols     = (node *) NULL;
        G->ecount   = 0;
        G->nrows    = 0;
        G->ncols    = 0;
        G->adjspace = (edge **) NULL;
        ILLptrworld_init (&G->intptrworld);
    }
}

static void free_graph (graph *G)
{
    if (G) {
        int total, onlist;
        ILL_IFFREE (G->edgelist, edge);
        ILL_IFFREE (G->rows, node);
        ILL_IFFREE (G->cols, node);
        ILL_IFFREE (G->adjspace, edge *);
        if (intptr_check_leaks (&G->intptrworld, &total, &onlist)) {
            fprintf (stderr, "WARNING: %d outstanding intptrs\n",
                     total - onlist);
        }
        ILLptrworld_delete (&G->intptrworld);
        init_graph (G);
    }
}

int ILLlp_sinfo_print (ILLlp_sinfo *s)
{
    int rval = 0;
    int i;
    ILLlpdata lp;
    char *sense = (char *) NULL;

    ILLlpdata_init (&lp);

    lp.nrows     = s->nrows;
    lp.ncols     = s->ncols;
    lp.nzcount   = s->nzcount;
    lp.objsense  = s->objsense;
    lp.obj       = s->obj;
    lp.rhs       = s->rhs;
    lp.lower     = s->lower;
    lp.upper     = s->upper;
    lp.A.matval  = s->A.matval;
    lp.A.matcnt  = s->A.matcnt;
    lp.A.matbeg  = s->A.matbeg;
    lp.A.matind  = s->A.matind;
    lp.rownames  = (char **) NULL;
    lp.colnames  = s->colnames;
    lp.objname   = (char *) NULL;
    lp.probname  = (char *) NULL;
    lp.intmarker = (char *) NULL;

    ILL_SAFE_MALLOC (sense, s->nrows, char);
    if (!sense) {
        fprintf (stderr, "out of memory in ILLlp_sinfo_print\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < s->nrows; i++) {
        sense[i] = 'E';
    }
    lp.sense = sense;

/*
    rval = ILLlpdata_writelp (&lp, (char *) NULL);
    ILL_CLEANUP_IF (rval);
*/

CLEANUP:

    ILL_IFFREE (sense, char);
    ILL_RETURN (rval, "ILLlp_sinfo_print");
}

void ILLlp_sinfo_init (ILLlp_sinfo *sinfo)
{
    if (sinfo) {
        sinfo->ncols = 0;
        sinfo->nrows = 0;
        sinfo->nzcount = 0;
        sinfo->rowsize = 0;
        sinfo->colsize = 0;
        sinfo->objsense = ILL_MIN;
        sinfo->obj      = (double *) NULL;
        sinfo->rhs      = (double *) NULL;
        sinfo->lower    = (double *) NULL;
        sinfo->upper    = (double *) NULL;

        ILLmatrix_init (&sinfo->A);

        sinfo->colnames = (char **) NULL;
    }
}

void ILLlp_sinfo_free (ILLlp_sinfo *sinfo)
{
    if (sinfo) {
        ILL_IFFREE (sinfo->obj, double);
        ILL_IFFREE (sinfo->rhs, double);
        ILL_IFFREE (sinfo->lower, double);
        ILL_IFFREE (sinfo->upper, double);
        ILLmatrix_free (&sinfo->A);
        if (sinfo->colnames) {
            int i;
            for (i = 0; i < sinfo->ncols; i++) {
                ILL_IFFREE (sinfo->colnames[i], char); 
            }
            ILL_IFFREE (sinfo->colnames, char *);
        }
        ILLlp_sinfo_init (sinfo);
    }
}

void ILLlp_predata_init (ILLlp_predata *pre)
{
    if (pre) {
        pre->opcount = 0;
        pre->opsize = 0;
        pre->oplist = (ILLlp_preop *) NULL;
        pre->r_nrows = 0;
        pre->r_ncols = 0;
        pre->colmap = (int *) NULL;
        pre->rowmap = (int *) NULL;
        pre->colscale  = (double *) NULL;
        pre->rowscale  = (double *) NULL;
        pre->colfixval = (double *) NULL;
        pre->rowfixval = (double *) NULL;
    }
}

void ILLlp_predata_free (ILLlp_predata *pre)
{
    if (pre) {
        int i;

        for (i = 0; i < pre->opcount; i++) {
            ILLlp_preop_free (&pre->oplist[i]);
        }
        ILL_IFFREE (pre->oplist, ILLlp_preop);
        ILL_IFFREE (pre->colmap, int);
        ILL_IFFREE (pre->rowmap, int);
        ILL_IFFREE (pre->colscale, double);
        ILL_IFFREE (pre->rowscale, double);
        ILL_IFFREE (pre->colfixval, double);
        ILL_IFFREE (pre->rowfixval, double);
        ILLlp_predata_init (pre);
    }
}

void ILLlp_preop_init (ILLlp_preop *op)
{
    if (op) {
        op->ptype = 0;
        op->rowindex = -1;
        op->colindex = -1;
        ILLlp_preline_init (&op->line);
    }
}

void ILLlp_preop_free (ILLlp_preop *op)
{
    if (op) {
        ILLlp_preline_free (&op->line);
        ILLlp_preop_init (op);
    }
}

void ILLlp_preline_init (ILLlp_preline *line)
{
    if (line) {
        line->rhs   = 0.0;
        line->obj   = 0.0;
        line->upper = 0.0;
        line->lower = 0.0;
        line->count = 0;
        line->ind   = (int *) NULL;
        line->val   = (double *) NULL;
    }
}

void ILLlp_preline_free (ILLlp_preline *line)
{
    if (line) {
        ILL_IFFREE (line->ind, int);
        ILL_IFFREE (line->val, double);
        ILLlp_preline_init (line);
    }
}
