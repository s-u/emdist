/* RCS_INFO = "$RCSfile: factor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#include <stdio.h>
#include <stdlib.h>
#include "iqsutil.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include "lpdefs.h"
#include "factor.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define SPARSE_FACTOR 0.05

#undef  RECORD
#undef  DEBUG
#undef  SOLVE_DEBUG

#undef  FACTOR_DEBUG
#undef  UPDATE_DEBUG

#define TRACK_FACTOR
#define NOTICE_BLOWUP

#undef  FACTOR_STATS
#undef  UPDATE_STATS
#undef  GROWTH_STATS

#undef  UPDATE_STUDY

#undef  SORT_RESULTS

#ifdef UPDATE_STUDY
int nupdate = 0;
double colspiketot = 0.0;
double rowspiketot = 0.0;
double permshifttot = 0.0;
double leftetatot = 0.0;
#endif

void ILLfactor_init_factor_work (factor_work *f)
{
    f->max_k        = 1000; /* must be less than 46340 (2^15.5) */
    f->fzero_tol    = 1e-15;
    f->szero_tol    = 1e-15;
    f->partial_tol  = 0.01;
    f->ur_space_mul = 2.0;
    f->uc_space_mul = 1.1;
    f->lc_space_mul = 1.1;
    f->er_space_mul = 1000.0;
    f->grow_mul     = 1.5;
    f->p            = 4;
    f->etamax       = 100;
    f->minmult      = 1e3;
    f->maxmult      = 1e5;
    f->updmaxmult   = 1e7;
    f->dense_fract  = 0.25;
    f->dense_min    = 25;

    f->partial_cur  = f->partial_tol;

    f->work_coef = (double *) NULL;
    f->work_indx = (int *) NULL;
    f->uc_inf    = (uc_info *) NULL;
    f->ur_inf    = (ur_info *) NULL;
    f->lc_inf    = (lc_info *) NULL;
    f->lr_inf    = (lr_info *) NULL;
    f->er_inf    = (er_info *) NULL;
    f->ucindx    = (int *) NULL;
    f->ucrind    = (int *) NULL;
    f->uccoef    = (double *) NULL;
    f->urindx    = (int *) NULL;
    f->urcind    = (int *) NULL;
    f->urcoef    = (double *) NULL;
    f->lcindx    = (int *) NULL;
    f->lccoef    = (double *) NULL;
    f->lrindx    = (int *) NULL;
    f->lrcoef    = (double *) NULL;
    f->erindx    = (int *) NULL;
    f->ercoef    = (double *) NULL;
    f->rperm     = (int *) NULL;
    f->rrank     = (int *) NULL;
    f->cperm     = (int *) NULL;
    f->crank     = (int *) NULL;
    f->dmat      = (double *) NULL;
    ILLsvector_init (&f->xtmp);
}

void ILLfactor_free_factor_work (factor_work *f)
{
#ifdef UPDATE_STUDY
    if (nupdate) {
        printf ("UPDATE STUDY: avg %d upd: %.2f col, %.2f row, %.2f lefteta, %.2f perm\n",
                nupdate, colspiketot/nupdate, rowspiketot/nupdate,
                leftetatot/nupdate, permshifttot/nupdate);
    }
#endif
    ILL_IFFREE (f->work_coef, double);
    ILL_IFFREE (f->work_indx, int);
    ILL_IFFREE (f->uc_inf, uc_info);
    ILL_IFFREE (f->ur_inf, ur_info);
    ILL_IFFREE (f->lc_inf, lc_info);
    ILL_IFFREE (f->lr_inf, lr_info);
    ILL_IFFREE (f->er_inf, er_info);
    ILL_IFFREE (f->ucindx, int);
    ILL_IFFREE (f->ucrind, int);
    ILL_IFFREE (f->uccoef, double);
    ILL_IFFREE (f->urindx, int);
    ILL_IFFREE (f->urcind, int);
    ILL_IFFREE (f->urcoef, double);
    ILL_IFFREE (f->lcindx, int);
    ILL_IFFREE (f->lccoef, double);
    ILL_IFFREE (f->lrindx, int);
    ILL_IFFREE (f->lrcoef, double);
    ILL_IFFREE (f->erindx, int);
    ILL_IFFREE (f->ercoef, double);
    ILL_IFFREE (f->rperm, int);
    ILL_IFFREE (f->rrank, int);
    ILL_IFFREE (f->cperm, int);
    ILL_IFFREE (f->crank, int);
    ILL_IFFREE (f->dmat, double);
    ILLsvector_free (&f->xtmp);
}

int ILLfactor_set_factor_iparam (factor_work *f, int param, int val)
{
    switch (param) {
    case QS_FACTOR_MAX_K:     f->max_k     = val; break;
    case QS_FACTOR_P:         f->p         = val; break;
    case QS_FACTOR_ETAMAX:    f->etamax    = val; break;
    case QS_FACTOR_DENSE_MIN: f->dense_min = val; break;
    default:
        fprintf (stderr, "Invalid param %d in ILLfactor_set_factor_iparam\n", param);
        return 1;
    }
    return 0;
}

int ILLfactor_set_factor_dparam (factor_work *f, int param, double val)
{
    switch (param) {
    case QS_FACTOR_FZERO_TOL:    f->fzero_tol    = val; break;
    case QS_FACTOR_SZERO_TOL:    f->szero_tol    = val; break;
    case QS_FACTOR_UR_SPACE_MUL: f->ur_space_mul = val; break;
    case QS_FACTOR_UC_SPACE_MUL: f->uc_space_mul = val; break;
    case QS_FACTOR_LC_SPACE_MUL: f->lc_space_mul = val; break;
    case QS_FACTOR_LR_SPACE_MUL: f->lr_space_mul = val; break;
    case QS_FACTOR_ER_SPACE_MUL: f->er_space_mul = val; break;
    case QS_FACTOR_GROW_MUL:     f->grow_mul     = val; break;
    case QS_FACTOR_MAXMULT:      f->maxmult      = val; break;      
    case QS_FACTOR_UPDMAXMULT:   f->updmaxmult   = val; break;
    case QS_FACTOR_DENSE_FRACT:  f->dense_fract  = val; break;
    case QS_FACTOR_PARTIAL_TOL:
        f->partial_tol  = val;
        f->partial_cur  = val;
        break;
    default:
        fprintf (stderr, "Invalid param %d in ILLfactor_set_factor_dparam\n", param);
        return 1;
    }
    return 0;
}

int ILLfactor_create_factor_work (factor_work *f, int dim)
{
    int i;
    int rval;
    
    f->dim       = dim;
    f->etacnt    = 0;
    ILL_SAFE_MALLOC (f->work_coef, dim, double);
    ILL_SAFE_MALLOC (f->work_indx, dim, int);
    ILL_SAFE_MALLOC (f->uc_inf, dim + (f->max_k + 1), uc_info);
    ILL_SAFE_MALLOC (f->ur_inf, dim + (f->max_k + 1), ur_info);
    ILL_SAFE_MALLOC (f->lc_inf, dim, lc_info);
    ILL_SAFE_MALLOC (f->lr_inf, dim, lr_info);
    ILL_SAFE_MALLOC (f->rperm, dim, int);
    ILL_SAFE_MALLOC (f->rrank, dim, int);
    ILL_SAFE_MALLOC (f->cperm, dim, int);
    ILL_SAFE_MALLOC (f->crank, dim, int);
    
    for (i=0; i<dim; i++) {
        f->work_coef[i] = 0.0;
        f->work_indx[i] = 0;
        f->uc_inf[i].nzcnt = 0;
        f->ur_inf[i].nzcnt = 0;
        f->lc_inf[i].nzcnt = 0;
        f->lr_inf[i].nzcnt = 0;
        f->rperm[i] = i;
        f->rrank[i] = i;
        f->cperm[i] = i;
        f->crank[i] = i;
    }
    for (i=0; i<=f->max_k; i++) {
        f->uc_inf[dim+i].nzcnt = i;
        f->uc_inf[dim+i].next  = dim + i;
        f->uc_inf[dim+i].prev  = dim + i;
        f->ur_inf[dim+i].nzcnt = i;
        f->ur_inf[dim+i].next  = dim + i;
        f->ur_inf[dim+i].prev  = dim + i;
    }

    rval = ILLsvector_alloc (&f->xtmp, dim);
    ILL_CLEANUP_IF (rval);
    
    rval = 0;

 CLEANUP:
    if (rval) {
        ILLfactor_free_factor_work (f);
    }
    ILL_RETURN (rval, "ILLfactor_create_factor_work");
}

static void dump_matrix (factor_work *f, int remaining)
{
    int       dim    = f->dim;
    ur_info  *ur_inf = f->ur_inf;
    uc_info  *uc_inf = f->uc_inf;
    lc_info  *lc_inf = f->lc_inf;
    lr_info  *lr_inf = f->lr_inf;
    er_info  *er_inf = f->er_inf;
    int       nzcnt;
    int       beg;
    
    int i;
    int j;
    
    for (i=0; i<dim; i++) {
        if (!remaining || ur_inf[i].next >= 0) {
            printf ("Row %d %d (max %.3f):", i, f->rrank[i], ur_inf[i].max);
            nzcnt = ur_inf[i].nzcnt;
            beg   = ur_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                if (j == ur_inf[i].pivcnt) {
                    printf (" |");
                }
                printf (" %.3f*%d", f->urcoef[beg + j], f->urindx[beg + j]);
                if (f->urcind) printf ("@%d", f->urcind[beg + j]);
            }
            printf ("\n");
        }
    }
    if (f->dmat) {
        int start = 0;
        if (remaining) start = f->stage - f->dense_base;
        printf ("Dcols at %d %d - %d    :", f->stage - f->dense_base,
                f->dense_base+start, f->nstages);
        for (j=start; j<f->dcols; j++) {
            printf (" %5d", f->cperm[j+f->dense_base]);
        }
        printf ("\n");
        for (i=start; i<f->drows; i++) {
            printf ("DRow %d %d (max %.3f):", i,
                    f->rperm[i + f->dense_base],
                    ur_inf[f->rperm[i + f->dense_base]].max);
            for (j=start; j<f->dcols; j++) {
                if (j == f->drows) {
                    printf (" |");
                }
                printf (" %.3f", f->dmat[i*f->dcols+j]);
            }
            printf ("\n");
        }
    }

    if (!remaining) {
        for (i=0; i<f->stage; i++) {
            printf ("L col %d:", lc_inf[i].c);
            nzcnt = lc_inf[i].nzcnt;
            beg   = lc_inf[i].cbeg;
            for (j=0; j<nzcnt; j++) {
                printf (" %.3f*%d", f->lccoef[beg + j], f->lcindx[beg + j]);
            }
            printf ("\n");
        }
        for (i=f->nstages; i<f->dim; i++) {
            printf ("L col %d:", lc_inf[i].c);
            nzcnt = lc_inf[i].nzcnt;
            beg   = lc_inf[i].cbeg;
            for (j=0; j<nzcnt; j++) {
                printf (" %.3f*%d", f->lccoef[beg + j], f->lcindx[beg + j]);
            }
            printf ("\n");
        }
        for (i=0; i<f->dim; i++) {
            if (!lr_inf[i].nzcnt) continue;
            printf ("L row %d:", lr_inf[i].r);
            nzcnt = lr_inf[i].nzcnt;
            beg   = lr_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                printf (" %.3f*%d", f->lrcoef[beg + j], f->lrindx[beg + j]);
            }
            printf ("\n");
        }
    }

    if (!remaining) {
        for (i=0; i<f->etacnt; i++) {
            printf ("Eta row %d:", f->er_inf[i].r);
            nzcnt = er_inf[i].nzcnt;
            beg   = er_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                printf (" %.3f*%d", f->ercoef[beg + j], f->erindx[beg + j]);
            }
            printf ("\n");
        }
    }
    
    for (i=0; i<dim; i++) {
        if (!remaining || uc_inf[i].next >= 0) {
            printf ("Col %d %d:", i, f->crank[i]);
            nzcnt = uc_inf[i].nzcnt;
            beg   = uc_inf[i].cbeg;
            for (j=0; j<nzcnt; j++) {
                if (f->uccoef != (double *) NULL) {
                    printf (" %.3f*%d", f->uccoef[beg + j], f->ucindx[beg + j]);
                    if (f->ucrind) printf ("@%d", f->ucrind[beg + j]);
                } else {
                    printf (" %d", f->ucindx[beg + j]);
                }
            }
            printf ("\n");
        }
    }

    if (!remaining) {
        printf ("rperm:");
        for (i=0; i<dim; i++) {
            if (i == f->nstages) printf ("|");
            if (i == f->stage) printf ("|");
            printf (" %d", f->rperm[i]);
        }
        printf ("\n");
        
        printf ("cperm:");
        for (i=0; i<dim; i++) {
            if (i == f->nstages) printf ("|");
            if (i == f->stage) printf ("|");
            printf (" %d", f->cperm[i]);
        }
        printf ("\n");
    }
    
    printf ("Rows by nzcnt:\n");
    for (i=0; i<=f->max_k; i++) {
        if (ur_inf[dim + i].next != dim + i) {
            printf ("%d:", i);
            for (j = ur_inf[dim + i].next; j != dim + i;
                 j = ur_inf[j].next) {
                printf (" %d", j);
            }
            printf ("\n");
        }
    }
    
    printf ("Cols by nzcnt:\n");
    for (i=0; i<=f->max_k; i++) {
        if (uc_inf[dim + i].next != dim + i) {
            printf ("%d:", i);
            for (j = uc_inf[dim + i].next; j != dim + i;
                 j = uc_inf[j].next) {
                printf (" %d", j);
            }
            printf ("\n");
        }
    }
    
    printf ("\n");
    fflush (stdout);
}

#ifdef SORT_RESULTS
static void sort_vector2 (int nzcnt, int *indx, double *coef)
{
  int i;
  int j;
  int itmp;
  double ctmp;

  for (i=1; i<nzcnt; i++) {
    itmp = indx[i];
    ctmp = coef[i];
    for (j=i; j>=1 && indx[j-1] > itmp; j--) {
      indx[j] = indx[j-1];
      coef[j] = coef[j-1];
    }
    indx[j] = itmp;
    coef[j] = ctmp;
  }
}

static void sort_vector (svector *x)
{
  sort_vector2 (x->nzcnt, x->indx, x->coef);
}
#endif

#ifdef DEBUG
static int check_matrix (factor_work *f)
{
    ur_info *ur_inf = f->ur_inf;
    uc_info *uc_inf = f->uc_inf;
    int rbeg;
    int nzcnt;
    int cbeg;
    int c;
    int r;
    int j;
    int nerr = 0;

    for (r=0; r<f->dim; r++) {
        nzcnt = ur_inf[r].nzcnt;
        rbeg  = ur_inf[r].rbeg;
        for (j=0; j<nzcnt; j++) {
            c = f->urindx[rbeg + j];
            cbeg = uc_inf[c].cbeg;
            if (f->ucindx[cbeg + f->urcind[rbeg + j]] != r) {
                printf ("index mismatch, row %d column %d\n", r, c);
                nerr++;
            }
            if (f->uccoef[cbeg + f->urcind[rbeg + j]] != f->urcoef[rbeg+j]) {
                printf ("coef mismatch, row %d column %d\n", r, c);
                nerr++;
            }
        }
    }
    if (f->urindx[f->ur_space] != 0) {
        printf ("last urindx entry %d != 0\n", f->urindx[f->ur_space]);
        nerr++;
    }

    for (c=0; c<f->dim; c++) {
        nzcnt = uc_inf[c].nzcnt;
        cbeg  = uc_inf[c].cbeg;
        for (j=0; j<nzcnt; j++) {
            r = f->ucindx[cbeg + j];
            rbeg = ur_inf[r].rbeg;
            if (f->urindx[rbeg + f->ucrind[cbeg + j]] != c) {
                printf ("index mismatch, column %d row %d\n", c, r);
                nerr++;
            }
            if (f->urcoef[rbeg + f->ucrind[cbeg + j]] != f->uccoef[cbeg+j]) {
                printf ("coef mismatch, column %d row %d\n", c, r);
                nerr++;
            }
        }
    }
    if (f->ucindx[f->uc_space] != 0) {
        printf ("last ucindx entry %d != 0\n", f->ucindx[f->uc_space]);
    }
    if (nerr) {
        fflush (stdout);
        dump_matrix (f, 0);
        return E_CHECK_FAILED;
    }
    return 0;
}
#endif

#ifdef FACTOR_STATS
static void dump_factor_stats (factor_work *f)
{
    int dim = f->dim;
    int ecnt = f->etacnt;
    ur_info *ur_inf = f->ur_inf;
    lc_info *lc_inf = f->lc_inf;
    er_info *er_inf = f->er_inf;
    double   *urcoef = f->urcoef;
    double   *lccoef = f->lccoef;
    double   *ercoef = f->ercoef;
    int lnzcnt = 0;
    int unzcnt = 0;
    int enzcnt = 0;
    int nzcnt;
    int beg;
    double umax;
    double lmax;
    double emax;
    int i;
    int j;

    umax = 0.0;
    for (i=0; i<dim; i++) {
        nzcnt = ur_inf[i].nzcnt;
        beg = ur_inf[i].rbeg;
        unzcnt += nzcnt;
        for (j=0; j<nzcnt; j++) {
            if (fabs(urcoef[beg+j]) > umax) umax = fabs(urcoef[beg+j]);
        }
    }
    lmax = 0.0;
    for (i=0; i<dim; i++) {
        nzcnt = lc_inf[i].nzcnt;
        beg = lc_inf[i].cbeg;
        lnzcnt += nzcnt;
        for (j=0; j<nzcnt; j++) {
            if (fabs(lccoef[beg+j]) > lmax) lmax = fabs(lccoef[beg+j]);
        }
    }
    emax = 0.0;
    for (i=0; i<ecnt; i++) {
        nzcnt = er_inf[i].nzcnt;
        beg = er_inf[i].rbeg;
        enzcnt += nzcnt;
        for (j=0; j<nzcnt; j++) {
            if (fabs(ercoef[beg+j]) > emax) emax = fabs(ercoef[beg+j]);
        }
    }
    printf ("factor U %d nzs %.3e max L %d nzs %.3e max E %d nzs %.3e max\n",
            unzcnt, umax, lnzcnt, lmax, enzcnt, emax);
    fflush (stdout);
}
#endif

static void clear_work (factor_work *f)
{
    int i;
    int dim = f->dim;
    double *work_coef = f->work_coef;

    for (i=0; i<dim; i++) {
        work_coef[i] = 0.0;
    }
}

static void load_row (factor_work *f, int r)
{
    double *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int    *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int     prow_nzcnt  = f->ur_inf[r].nzcnt;
    double *work_coef   = f->work_coef;
    int    *work_indx   = f->work_indx;
    int i;
    int j;
    
    for (i=0; i<prow_nzcnt; i++) {
        j = prow_urindx[i];
        work_coef[j] = prow_urcoef[i];
        work_indx[j] = 1;
    }
}

static void clear_row (factor_work *f, int r)
{
    int    *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int     prow_nzcnt  = f->ur_inf[r].nzcnt;
    double *work_coef   = f->work_coef;
    int    *work_indx   = f->work_indx;
    int i;
    int j;
    
    for (i=0; i<prow_nzcnt; i++) {
        j = prow_urindx[i];
        work_coef[j] = 0.0;
        work_indx[j] = 0;
    }
}

static int make_ur_space (factor_work *f, int space)
{
    double   *new_urcoef = (double *) NULL;
    int      *new_urindx = (int *) NULL;
    int      *new_urcind = (int *) NULL;
    double   *urcoef     = f->urcoef;
    int      *urindx     = f->urindx;
    int      *urcind     = f->urcind;
    int       minspace;
    ur_info  *ur_inf     = f->ur_inf;
    int       dim        = f->dim;
    int       new_nzcnt  = 0;
    int       rbeg;
    int       nzcnt;
    int i;
    int j;
    int rval;

    minspace = f->ur_space;
    nzcnt = space;
    for (i=0; i<dim; i++) nzcnt += ur_inf[i].nzcnt;
    while (nzcnt * 2 >= minspace) {
        minspace *= f->grow_mul;
    }

#ifdef GROWTH_STATS
    printf ("make_ur_space growing from %d to %d...", f->ur_space, minspace);
    fflush (stdout);
#endif
    
    ILL_SAFE_MALLOC (new_urcoef, minspace, double);
    ILL_SAFE_MALLOC (new_urindx, minspace+1, int);

    if (urcind) {
        ILL_SAFE_MALLOC (new_urcind, minspace, int);
    }

    if (urcind) {
        for (j=0; j<dim; j++) {
            rbeg = ur_inf[j].rbeg;
            nzcnt = ur_inf[j].nzcnt;
            ur_inf[j].rbeg = new_nzcnt;
            for (i=0; i<nzcnt; i++) {
                new_urindx[new_nzcnt] = urindx[rbeg+i];
                new_urcoef[new_nzcnt] = urcoef[rbeg+i];
                new_urcind[new_nzcnt] = urcind[rbeg+i];
                new_nzcnt++;
            }
        }
    } else {
        for (j=0; j<dim; j++) {
            rbeg = ur_inf[j].rbeg;
            nzcnt = ur_inf[j].nzcnt;
            ur_inf[j].rbeg = new_nzcnt;
            for (i=0; i<nzcnt; i++) {
                new_urindx[new_nzcnt] = urindx[rbeg+i];
                new_urcoef[new_nzcnt] = urcoef[rbeg+i];
                new_nzcnt++;
            }
        }
    }
    
    for (i=new_nzcnt; i<minspace; i++) {
        new_urindx[i] = -1;
    }
    new_urindx[minspace] = 0;

    ILL_IFFREE (f->urcoef, double);
    f->urcoef = new_urcoef;
    new_urcoef = (double *) NULL;

    ILL_IFFREE (f->urindx, int);
    f->urindx = new_urindx;
    new_urindx = (int *) NULL;

    ILL_IFFREE (f->urcind, int);
    f->urcind = new_urcind;
    new_urcind = (int *) NULL;

    f->ur_freebeg = new_nzcnt;
    f->ur_space = minspace;

#ifdef GROWTH_STATS
    printf ("%d nonzeros\n", new_nzcnt);
    fflush (stdout);
    dump_factor_stats (f);
#endif
    
    rval = 0;

 CLEANUP:
    ILL_IFFREE (new_urcoef, double);
    ILL_IFFREE (new_urindx, int);
    ILL_IFFREE (new_urcind, int);
    ILL_RETURN (rval, "make_ur_space");
}

static int make_uc_space (factor_work *f, int space)
{
    double   *new_uccoef = (double *) NULL;
    int      *new_ucindx = (int *) NULL;
    int      *new_ucrind = (int *) NULL;
    int       uc_freebeg = f->uc_freebeg;
    double   *uccoef     = f->uccoef;
    int      *ucindx     = f->ucindx;
    int      *ucrind     = f->ucrind;
    int       minspace   = uc_freebeg + space;
    uc_info  *uc_inf     = f->uc_inf;
    int       dim        = f->dim;
    int       new_nzcnt  = 0;
    int       cbeg;
    int       nzcnt;
    int i;
    int j;
    int rval;
    
    if (f->uc_space * f->grow_mul > minspace) {
        minspace = f->uc_space * f->grow_mul;
    }

#ifdef GROWTH_STATS
    printf ("make_uc_space growing from %d to %d...", f->uc_space, minspace);
    fflush (stdout);
#endif
    
    ILL_SAFE_MALLOC (new_ucindx, minspace+1, int);

    if (ucrind) {
        ILL_SAFE_MALLOC (new_uccoef, minspace, double);
        ILL_SAFE_MALLOC (new_ucrind, minspace, int);
    }

    if (ucrind) {
        for (j=0; j<dim; j++) {
            cbeg = uc_inf[j].cbeg;
            nzcnt = uc_inf[j].nzcnt;
            uc_inf[j].cbeg = new_nzcnt;
            for (i=0; i<nzcnt; i++) {
                new_ucindx[new_nzcnt] = ucindx[cbeg+i];
                new_uccoef[new_nzcnt] = uccoef[cbeg+i];
                new_ucrind[new_nzcnt] = ucrind[cbeg+i];
                new_nzcnt++;
            }
        }
    } else {
        for (j=0; j<dim; j++) {
            cbeg = uc_inf[j].cbeg;
            nzcnt = uc_inf[j].nzcnt;
            uc_inf[j].cbeg = new_nzcnt;
            for (i=0; i<nzcnt; i++) {
                new_ucindx[new_nzcnt] = ucindx[cbeg+i];
                new_nzcnt++;
            }
        }
    }

    for (i=new_nzcnt; i<minspace; i++) {
        new_ucindx[i] = -1;
    }
    new_ucindx[minspace] = 0;

    ILL_IFFREE (f->uccoef, double);
    f->uccoef = new_uccoef;
    new_uccoef = (double *) NULL;

    ILL_IFFREE (f->ucindx, int);
    f->ucindx = new_ucindx;
    new_ucindx = (int *) NULL;

    ILL_IFFREE (f->ucrind, int);
    f->ucrind = new_ucrind;
    new_ucrind = (int *) NULL;
    
    f->uc_freebeg = new_nzcnt;
    f->uc_space = minspace;

#ifdef GROWTH_STATS
    printf ("%d nonzeros\n", new_nzcnt);
    fflush (stdout);
    dump_factor_stats (f);
#endif

    rval = 0;

 CLEANUP:
    ILL_IFFREE (new_uccoef, double);
    ILL_IFFREE (new_ucindx, int);
    ILL_IFFREE (new_ucrind, int);
    ILL_RETURN (rval, "make_uc_space");
}

static int make_lc_space (factor_work *f, int space)
{
    double   *new_lccoef = (double *) NULL;
    int      *new_lcindx = (int *) NULL;
    int       lc_freebeg = f->lc_freebeg;
    double   *lccoef     = f->lccoef;
    int      *lcindx     = f->lcindx;
    int       minspace   = lc_freebeg + space;
    int i;
    int rval;
    
    if (f->lc_space * f->grow_mul > minspace) {
        minspace = f->lc_space * f->grow_mul;
    }

#ifdef GROWTH_STATS
    printf ("make_lc_space growing from %d to %d...", f->lc_space, minspace);
    fflush (stdout);
#endif
    
    ILL_SAFE_MALLOC (new_lccoef, minspace, double);
    ILL_SAFE_MALLOC (new_lcindx, minspace, int);

    for (i=0; i<lc_freebeg; i++) {
        new_lccoef[i] = lccoef[i];
        new_lcindx[i] = lcindx[i];
    }
    
    ILL_IFFREE (lccoef, double);
    f->lccoef = new_lccoef;
    new_lccoef = (double *) NULL;

    ILL_IFFREE (lcindx, int);
    f->lcindx = new_lcindx;
    new_lcindx = (int *) NULL;

    f->lc_space = minspace;
    
#ifdef GROWTH_STATS
    printf ("done\n");
    fflush (stdout);
    dump_factor_stats (f);
#endif
    
    rval = 0;

 CLEANUP:
    ILL_IFFREE (new_lccoef, double);
    ILL_IFFREE (new_lcindx,  int);
    ILL_RETURN (rval, "make_lc_space");
}

static void set_col_nz (factor_work *f, int c)
{
    uc_info  *uc_inf = f->uc_inf;
    int       nzcnt  = uc_inf[c].nzcnt;
    int       max_k  = f->max_k;
    int       dim    = f->dim;
    
    if (uc_inf[c].next >= 0) {
        uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
        uc_inf[uc_inf[c].prev].next = uc_inf[c].next;
        
        if (nzcnt >= max_k) nzcnt = max_k;
        uc_inf[c].next = uc_inf[dim + nzcnt].next;
        uc_inf[c].prev = dim + nzcnt;
        uc_inf[dim+nzcnt].next = c;
        uc_inf[uc_inf[c].next].prev = c;
    }
}

static void set_row_nz (factor_work *f, int r)
{
    ur_info  *ur_inf = f->ur_inf;
    int       nzcnt  = ur_inf[r].pivcnt;
    int       max_k  = f->max_k;
    int       dim    = f->dim;
    
    if (ur_inf[r].next >= 0) {
        ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
        ur_inf[ur_inf[r].prev].next = ur_inf[r].next;
        
        if (nzcnt >= max_k) nzcnt = max_k;
        ur_inf[r].next = ur_inf[dim + nzcnt].next;
        ur_inf[r].prev = dim + nzcnt;
        ur_inf[dim+nzcnt].next = r;
        ur_inf[ur_inf[r].next].prev = r;
    }
}

static void remove_col_nz (factor_work *f, int r, int c)
{
    uc_info  *uc_inf = f->uc_inf;
    int      *ucindx = f->ucindx + uc_inf[c].cbeg;
    int       nzcnt  = uc_inf[c].nzcnt;
    int i;

    for (i=0; i<nzcnt; i++) {
        if (ucindx[i] == r) {
            --nzcnt;
            ucindx[i] = ucindx[nzcnt];
            ucindx[nzcnt] = -1;
            break;
        }
    }
    uc_inf[c].nzcnt = nzcnt;

    set_col_nz (f, c);
}

static void remove_row_nz (factor_work *f, int r, int c)
{
    ur_info  *ur_inf = f->ur_inf;
    int      *urindx = f->urindx + ur_inf[r].rbeg;
    double   *urcoef = f->urcoef + ur_inf[r].rbeg;
    int       pivcnt = ur_inf[r].pivcnt;
    double max = 0.0;
    int tind;
    double tcoef;
    int i;

    for (i=0; i<pivcnt; i++) {
        if (urindx[i] == c) {
            --pivcnt;
            ILL_SWAP(urindx[i], urindx[pivcnt], tind);
            ILL_SWAP(urcoef[i], urcoef[pivcnt], tcoef);
            --i;
        } else {
            if (fabs(urcoef[i]) > max) max = fabs(urcoef[i]);
        }
    }
    ur_inf[r].pivcnt = pivcnt;
    ur_inf[r].max = max;

    set_row_nz (f, r);
}

static int add_col_nz (factor_work *f, int r, int c)
{
    uc_info  *uc_inf     = f->uc_inf;
    int       cbeg       = uc_inf[c].cbeg;
    int       nzcnt      = uc_inf[c].nzcnt;
    int       uc_freebeg = f->uc_freebeg;
    int      *ucindx     = f->ucindx;
    int i;
    int rval = 0;

    if (uc_inf[c].next == -1) {
        return 0;
    }
    
    if (ucindx[cbeg+nzcnt] == -1) {
        ucindx[cbeg+nzcnt] = r;
        uc_inf[c].nzcnt++;
        if (nzcnt + cbeg == uc_freebeg) {
            f->uc_freebeg = uc_freebeg+1;
        }
    } else {
        if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
            rval = make_uc_space (f, nzcnt+1);
            ILL_CLEANUP_IF (rval);
            uc_freebeg = f->uc_freebeg;
            cbeg = uc_inf[c].cbeg;
            ucindx = f->ucindx;
        }
        for (i=0; i<nzcnt; i++) {
            ucindx[uc_freebeg+i] = ucindx[cbeg+i];
            ucindx[cbeg+i] = -1;
        }
        ucindx[uc_freebeg+nzcnt] = r;
        uc_inf[c].cbeg = uc_freebeg;
        uc_inf[c].nzcnt++;
        f->uc_freebeg = uc_freebeg + nzcnt + 1;
    }

    set_col_nz (f, c);
 CLEANUP:
    ILL_RETURN (rval, "add_col_nz");
}

static void disable_col (factor_work *f, int c)
{
    uc_info *uc_inf = f->uc_inf;

    if (uc_inf[c].next >= 0) {
        uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
        uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

        uc_inf[c].next = -2;
        uc_inf[c].prev = -2;
    }
}

static void remove_col (factor_work *f, int c)
{
    uc_info  *uc_inf = f->uc_inf;
    int       cbeg   = uc_inf[c].cbeg;
    int       nzcnt  = uc_inf[c].nzcnt;
    int      *ucindx = f->ucindx;
    int i;

    for (i=0; i<nzcnt; i++) {
        ucindx[cbeg+i] = -1;
    }
    uc_inf[c].cbeg = 0;
    uc_inf[c].nzcnt = 0;

    if (uc_inf[c].next >= 0) {
        uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
        uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

        uc_inf[c].next = -1;
        uc_inf[c].prev = -1;
    }
}

static void remove_row (factor_work *f, int r)
{
    ur_info *ur_inf = f->ur_inf;

    if (ur_inf[r].next >= 0) {
        ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
        ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

        ur_inf[r].next = -1;
        ur_inf[r].prev = -1;
    }
}

static double find_coef (factor_work *f, int r, int c)
{
    double *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int    *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int i;
    int prow_nzcnt = f->ur_inf[r].nzcnt;

    for (i=0; i<prow_nzcnt; i++) {
        if (prow_urindx[i] == c) {
            return prow_urcoef[i];
        }
    }
    fprintf (stderr, "Coefficient not found\n");
    return 0.0;
}

static int elim_row (factor_work *f, int elim_r, int r, int c,
        double *p_pivot_coef)
{
    ur_info  *ur_inf      = f->ur_inf;
    double   *work_coef   = f->work_coef;
    int      *work_indx   = f->work_indx;
    double    elim_coef   = find_coef (f, r, c) / work_coef[c];
    double   *urcoef      = f->urcoef;
    int      *urindx      = f->urindx;
    int       prow_beg    = ur_inf[r].rbeg;
    int       prow_nzcnt  = ur_inf[r].nzcnt;
    int       prow_pivcnt = ur_inf[r].pivcnt;
    double    fzero_tol   = f->fzero_tol;
    int       fill        = ur_inf[elim_r].nzcnt;
    int       cancel      = 0;
    double    max         = 0.0;
    int       erow_beg;
    int       erow_nzcnt;
    int       erow_pivcnt;
    double x;
    int i;
    int j;
    int rval = 0;

    *p_pivot_coef = elim_coef;

    for (i=0; i<prow_nzcnt; i++) {
        j = urindx[prow_beg+i];
        if (work_indx[j] == 1) {
            x = urcoef[prow_beg+i] - elim_coef * work_coef[j];
            if ((x <= fzero_tol && x >= -fzero_tol) || j == c) {
                cancel++;
                if (j != c) {
                    remove_col_nz (f, r, j);
                }
                if (i < prow_pivcnt) {
                    prow_pivcnt--;
                    prow_nzcnt--;
                    urindx[prow_beg+i] = urindx[prow_beg+prow_pivcnt];
                    urcoef[prow_beg+i] = urcoef[prow_beg+prow_pivcnt];
                    if (prow_pivcnt != prow_nzcnt) {
                        urindx[prow_beg+prow_pivcnt] = urindx[prow_beg+prow_nzcnt];
                        urcoef[prow_beg+prow_pivcnt] = urcoef[prow_beg+prow_nzcnt];
                    }
                } else {
                    prow_nzcnt--;
                    urindx[prow_beg+i] = urindx[prow_beg+prow_nzcnt];
                    urcoef[prow_beg+i] = urcoef[prow_beg+prow_nzcnt];
                }
                urindx[prow_beg + prow_nzcnt] = -1;
                i--;
            } else {
                urcoef[prow_beg+i] = x;
                if (i < prow_pivcnt && fabs(x) > max) max = fabs(x);
            }
            work_indx[j] = 0;
            fill--;
        } else {
            if (i < prow_pivcnt && fabs(urcoef[prow_beg+i]) > max) {
                max = fabs(urcoef[prow_beg+i]);
            }
        }
    }

    if (fill > 0) {
        ur_inf[r].nzcnt = prow_nzcnt;
        ur_inf[r].pivcnt = prow_pivcnt;
        if (fill > cancel) {
            int ur_freebeg = f->ur_freebeg;

            if (ur_freebeg + prow_nzcnt + fill >= f->ur_space) {
                rval = make_ur_space (f, prow_nzcnt + fill);
                ILL_CLEANUP_IF (rval);
                urcoef = f->urcoef;
                urindx = f->urindx;
                ur_freebeg = f->ur_freebeg;
                prow_beg = f->ur_inf[r].rbeg;
            }
            for (i=0; i<prow_nzcnt; i++) {
                urindx[ur_freebeg+i] = urindx[prow_beg+i];
                urcoef[ur_freebeg+i] = urcoef[prow_beg+i];
                urindx[prow_beg+i] = -1;
            }
            ur_inf[r].rbeg = ur_freebeg;
            f->ur_freebeg = ur_freebeg + prow_nzcnt + fill;
            prow_beg = ur_freebeg;
        }

        erow_beg    = ur_inf[elim_r].rbeg;
        erow_nzcnt  = ur_inf[elim_r].nzcnt;
        erow_pivcnt = ur_inf[elim_r].pivcnt;

        for (i=0; i<erow_pivcnt; i++) {
            j = urindx[erow_beg + i];
            if (work_indx[j] == 1) {
                x = -urcoef[erow_beg+i] * elim_coef;
                if (x > fzero_tol || x < -fzero_tol) {
                    rval = add_col_nz (f, r, j);
                    ILL_CLEANUP_IF (rval);
                    if (prow_pivcnt != prow_nzcnt) {
                        urindx[prow_beg+prow_nzcnt] =
                                                urindx[prow_beg + prow_pivcnt];
                        urcoef[prow_beg+prow_nzcnt] =
                                                urcoef[prow_beg + prow_pivcnt];
                    }
                    urindx[prow_beg+prow_pivcnt] = j;
                    urcoef[prow_beg+prow_pivcnt] = x;
                    if (fabs(x) > max) max = fabs(x);
                    prow_pivcnt++;
                    prow_nzcnt++;
                }
            } else {
                work_indx[j] = 1;
            }
        }
        for (i=erow_pivcnt; i<erow_nzcnt; i++) {
            j = urindx[erow_beg + i];
            if (work_indx[j] == 1) {
                x = -urcoef[erow_beg+i] * elim_coef;
                if (x > fzero_tol || x < -fzero_tol) {
                    rval = add_col_nz (f, r, j);
                    ILL_CLEANUP_IF (rval);
                    urindx[prow_beg+prow_nzcnt] = j;
                    urcoef[prow_beg+prow_nzcnt] = x;
                    prow_nzcnt++;
                }
            } else {
                work_indx[j] = 1;
            }
        }
    } else {
        erow_nzcnt = ur_inf[elim_r].nzcnt;
        erow_beg   = ur_inf[elim_r].rbeg;
        for (i=0; i<erow_nzcnt; i++) {
            j = urindx[erow_beg+i];
            work_indx[j] = 1;
        }
    }
            
    ur_inf[r].nzcnt = prow_nzcnt;
    ur_inf[r].pivcnt = prow_pivcnt;
    ur_inf[r].max = max;

    set_row_nz (f, r);
 CLEANUP:
    ILL_RETURN (rval, "elim_row");
}

#define SETPERM(f,s,r,c) {                    \
        f->rperm[f->rrank[r]] = f->rperm[s];  \
        f->rrank[f->rperm[s]] = f->rrank[r];  \
        f->rperm[s] = r;                      \
        f->rrank[r] = s;                      \
                                              \
        f->cperm[f->crank[c]] = f->cperm[s];  \
        f->crank[f->cperm[s]] = f->crank[c];  \
        f->cperm[s] = c;                      \
        f->crank[c] = s;                      \
}

static int elim (factor_work *f, int r, int c)
{
    uc_info *uc_inf = f->uc_inf;
    ur_info *ur_inf = f->ur_inf;
    lc_info *lc_inf = f->lc_inf;
    int *urindx;
    int *ucindx;
    int *lcindx;
    double *urcoef;
    double *lccoef;
    double pivot_coef;
    int nzcnt;
    int lc_freebeg;
    int s = f->stage;
    int i;
    int j;
    int rval = 0;

    if (uc_inf[c].nzcnt == 1) {
        /* col singleton */
        SETPERM(f,s,r,c);

        lc_inf[s].cbeg = -1;
        lc_inf[s].c = r;
        lc_inf[s].nzcnt = 0;
        f->stage++;

        urindx = f->urindx + ur_inf[r].rbeg;
        urcoef = f->urcoef + ur_inf[r].rbeg;
        nzcnt = ur_inf[r].nzcnt;
        for (i=0; i<nzcnt; i++) {
            j = urindx[i];
            remove_col_nz (f, r, j);
            if (j == c) {
                urindx[i] = urindx[0];
                urindx[0] = c;
                ILL_SWAP(urcoef[0], urcoef[i], pivot_coef);
            }
        }
        remove_row (f, r);
        remove_col (f, c);
    } else if (ur_inf[r].nzcnt == 1) {
        /* row singleton */
        --(f->nstages);
        SETPERM(f,f->nstages,r,c);

        lc_inf[f->nstages].cbeg = -1;
        lc_inf[f->nstages].c = r;
        lc_inf[f->nstages].nzcnt = 0;

        ucindx = f->ucindx + uc_inf[c].cbeg;
        nzcnt = uc_inf[c].nzcnt;
        for (i=0; i<nzcnt; i++) {
            j = ucindx[i];
            remove_row_nz (f, j, c);
        }
        remove_row (f, r);
        remove_col (f, c);
    } else {
        SETPERM(f,s,r,c);
        f->stage++;

        nzcnt = uc_inf[c].nzcnt;
        if (f->lc_freebeg + nzcnt >= f->lc_space) {
            rval = make_lc_space (f, nzcnt);
            ILL_CLEANUP_IF (rval);
        }
        lc_freebeg = f->lc_freebeg;
        lc_inf[s].cbeg = lc_freebeg;
        lc_inf[s].c = r;
        lcindx = f->lcindx;
        lccoef = f->lccoef;
        load_row (f, r);
        ucindx = f->ucindx + uc_inf[c].cbeg;
        for (i=0; i<nzcnt; i++) {
            j = f->ucindx[uc_inf[c].cbeg + i];
            if (j != r) {
                rval = elim_row (f, r, j, c, &pivot_coef);
                ILL_CLEANUP_IF (rval);
                lcindx[lc_freebeg] = j;
                lccoef[lc_freebeg] = pivot_coef;
                lc_freebeg++;
#ifdef TRACK_FACTOR
                if (fabs(pivot_coef) > f->maxelem_factor) {
                    f->maxelem_factor = fabs(pivot_coef);
                }
                if (ur_inf[r].max > f->maxelem_factor) {
                    f->maxelem_factor = ur_inf[r].max;
                }
#endif /* TRACK_FACTOR */
            }
        }
        lc_inf[s].nzcnt = lc_freebeg - lc_inf[s].cbeg;
        f->lc_freebeg = lc_freebeg;

        clear_row (f, r);

        urindx = f->urindx + ur_inf[r].rbeg;
        urcoef = f->urcoef + ur_inf[r].rbeg;
        nzcnt = ur_inf[r].nzcnt;
        for (i=0; i<nzcnt; i++) {
            j = urindx[i];
            remove_col_nz (f, r, j);
            if (j == c) {
                urindx[i] = urindx[0];
                urindx[0] = c;
                ILL_SWAP(urcoef[0], urcoef[i], pivot_coef);
            }
        }
        remove_row (f, r);
        remove_col (f, c);
    }
 CLEANUP:
    ILL_RETURN (rval, "elim");
}

static void find_pivot_column (factor_work *f, int c, int *p_r)
{
    uc_info  *uc_inf = f->uc_inf;
    ur_info  *ur_inf = f->ur_inf;
    int      *ucindx = f->ucindx;
    int       nzcnt  = uc_inf[c].nzcnt;
    int       cbeg   = uc_inf[c].cbeg;
    double    partial_cur = f->partial_cur;
    int       bestnz = -1;
    int i;
    int r;

    *p_r = -1;
    for (i=0; i<nzcnt; i++) {
        r = ucindx[cbeg+i];
        if ((bestnz == -1 || ur_inf[r].pivcnt < bestnz) &&
            fabs(find_coef (f, r, c)) >= partial_cur * ur_inf[r].max) {
            bestnz = ur_inf[r].pivcnt;
            *p_r = r;
        }
    }
}

static void find_pivot_row (factor_work *f, int r, int *p_c)
{
    uc_info  *uc_inf = f->uc_inf;
    ur_info  *ur_inf = f->ur_inf;
    int      *urindx = f->urindx;
    double   *urcoef = f->urcoef;
    int       pivcnt = ur_inf[r].pivcnt;
    int       rbeg   = ur_inf[r].rbeg;
    double    thresh = f->partial_cur * ur_inf[r].max;
    int       bestnz = -1;
    int i;
    int c;

    *p_c = -1;
    for (i=0; i<pivcnt; i++) {
        c = urindx[rbeg+i];
        if ((bestnz == -1 || uc_inf[c].nzcnt < bestnz) &&
            fabs(urcoef[rbeg+i]) >= thresh) {
            bestnz = uc_inf[c].nzcnt;
            *p_c = c;
        }
    }
}

static int find_pivot (factor_work *f, int *p_r, int *p_c)
{
    uc_info  *uc_inf = f->uc_inf;
    ur_info  *ur_inf = f->ur_inf;
    int       dim    = f->dim;
    int       max_k  = f->max_k;
    int       p      = f->p;
    int       c;
    int       r;
    int       mm;
    int       n;
    int m;
    int k;

    if (uc_inf[dim+1].next != dim+1) {
        c = uc_inf[dim+1].next;
        r = f->ucindx[uc_inf[c].cbeg];
        *p_c = c;
        *p_r = r;
        return 0;
    } else if (ur_inf[dim+1].next != dim+1) {
        r = ur_inf[dim+1].next;
        c = f->urindx[ur_inf[r].rbeg];
        *p_c = c;
        *p_r = r;
        return 0;
    }
    mm = 0;
    n = 0;
    *p_r = -1;
    *p_c = -1;
    for (k=2; k<=max_k && (mm == 0 || mm > (k-1)*(k-1)); k++) {
        if (uc_inf[dim+k].next != dim+k) {
            for (c = uc_inf[dim+k].next; c != dim+k; c = uc_inf[c].next) {
                find_pivot_column (f, c, &r);
                if (r >= 0) {
                    m = (uc_inf[c].nzcnt-1) * (ur_inf[r].pivcnt-1);
                    if (mm == 0 || m < mm) {
                        mm = m;
                        *p_c = c;
                        *p_r = r;
                        if (mm <= (k-1)*(k-1)) {
                            return 0;
                        }
                    }
                } else {
                    c = uc_inf[c].prev;
                    disable_col (f, uc_inf[c].next);
                }
                n++;
                if (n >= p && mm != 0) {
                    return 0;
                }
            }
        }

        if (ur_inf[dim+k].next != dim+k) {
            for (r = ur_inf[dim+k].next; r != dim+k; r = ur_inf[r].next) {
                find_pivot_row (f, r, &c);
                if (c >= 0) {
                    m = (uc_inf[c].nzcnt-1) * (ur_inf[r].pivcnt-1);
                    if (mm == 0 || m < mm) {
                        mm = m;
                        *p_c = c;
                        *p_r = r;
                        if (mm <= k*(k-1)) {
                            return 0;
                        }
                    }
                }
                n++;
                if (n >= p && mm != 0) {
                    return 0;
                }
            }
        }
    }
    if (mm != 0) {
        return 0;
    } else {
        fprintf (stderr, "No acceptable pivot found\n");
        return E_NO_PIVOT;
    }
}

static int create_factor_space (factor_work *f)
{
    uc_info  *uc_inf = f->uc_inf;
    ur_info  *ur_inf = f->ur_inf;
    int       dim    = f->dim;
    int nzcnt;
    int i;
    int rval;

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        nzcnt += ur_inf[i].nzcnt;
    }
    
    if (f->ucindx == (int *) NULL) {
        f->uc_space = nzcnt * f->uc_space_mul;
        ILL_SAFE_MALLOC (f->ucindx, f->uc_space+1, int);
    }
    
    if (f->urindx == (int *) NULL ||
        f->urcoef == (double *) NULL) {
        ILL_IFFREE (f->urindx, int);
        ILL_IFFREE (f->urcoef, double);
        f->ur_space = nzcnt * f->ur_space_mul;
        ILL_SAFE_MALLOC (f->urindx, f->ur_space+1, int);
        ILL_SAFE_MALLOC (f->urcoef, f->ur_space, double);
    }

    if (f->lcindx == (int *) NULL ||
        f->lccoef == (double *) NULL) {
        ILL_IFFREE (f->lcindx, int);
        ILL_IFFREE (f->lccoef, double);
        f->lc_space = nzcnt * f->lc_space_mul;
        ILL_SAFE_MALLOC (f->lcindx, f->lc_space, int);
        ILL_SAFE_MALLOC (f->lccoef, f->lc_space, double);
    }

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        ur_inf[i].rbeg = nzcnt;
        nzcnt += ur_inf[i].nzcnt;
        ur_inf[i].nzcnt = ur_inf[i].rbeg;
    }
    f->ur_freebeg = nzcnt;

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        uc_inf[i].cbeg = nzcnt;
        nzcnt += uc_inf[i].nzcnt;
        uc_inf[i].nzcnt = uc_inf[i].cbeg;
    }
    f->uc_freebeg = nzcnt;

    f->lc_freebeg = 0;

    rval = 0;
 CLEANUP:
    ILL_RETURN (rval, "create_factor_space");
}
    
static int init_matrix (factor_work *f, int *basis, int *cbeg, int *clen,
        int *in_ucindx, double *in_uccoef)
{
    uc_info  *uc_inf = f->uc_inf;
    ur_info  *ur_inf = f->ur_inf;
    int       dim    = f->dim;
    int       max_k  = f->max_k;
    int      *ucindx;
    int      *urindx;
    double   *urcoef;
    int       nzcnt;
    int       beg;
    double    fzero_tol = f->fzero_tol;
    double    v;
    double max;
    int i;
    int j;
    int r;
    int rval = 0;

    for (i=0; i<dim; i++) {
        ur_inf[i].nzcnt = 0;
    }
    for (i=0; i<dim; i++) {
        nzcnt = clen[basis[i]];
        beg = cbeg[basis[i]];
        uc_inf[i].nzcnt = nzcnt;
        for (j=0; j<nzcnt; j++) {
            r = in_ucindx[beg+j];
            ur_inf[r].nzcnt++;
        }
    }

    rval = create_factor_space (f);
    ILL_CLEANUP_IF (rval);

    urindx = f->urindx;
    ucindx = f->ucindx;
    urcoef = f->urcoef;
    
    for (i=0; i<dim; i++) {
        nzcnt = clen[basis[i]];
        beg = cbeg[basis[i]];
        for (j=0; j<nzcnt; j++) {
            v = in_uccoef[beg+j];
            if (v <= fzero_tol && v >= -fzero_tol) continue;
            r = in_ucindx[beg+j];
            ucindx[uc_inf[i].nzcnt++] = r;
            urindx[ur_inf[r].nzcnt] = i;
            urcoef[ur_inf[r].nzcnt] = v;
            ur_inf[r].nzcnt++;
        }
    }

    for (i=0; i<dim; i++) {
        uc_inf[i].nzcnt -= uc_inf[i].cbeg;
        ur_inf[i].nzcnt -= ur_inf[i].rbeg;
    }

    j = f->uc_space;
    for (i=f->uc_freebeg; i < j; i++) {
        ucindx[i] = -1;
    }
    ucindx[j] = 0;

    j = f->ur_space;
    for (i=f->ur_freebeg; i<j; i++) {
        urindx[i] = -1;
    }
    urindx[j] = 0;
    
    for (i=0; i<dim; i++) {
        nzcnt = ur_inf[i].nzcnt;
        ur_inf[i].pivcnt = nzcnt;
        beg = ur_inf[i].rbeg;
        max = 0.0;
        for (j=0; j<nzcnt; j++) {
            if (fabs(urcoef[beg+j]) > max) max = fabs(urcoef[beg+j]);
        }
        ur_inf[i].max = max;
    }
    
    for (i=0; i<=max_k; i++) {
        ur_inf[dim+i].next = dim+i;
        ur_inf[dim+i].prev = dim+i;
        uc_inf[dim+i].next = dim+i;
        uc_inf[dim+i].prev = dim+i;
    }

    for (i=0; i<dim; i++) {
        nzcnt = uc_inf[i].nzcnt;
        if (nzcnt >= max_k) nzcnt = max_k;
        uc_inf[i].next = uc_inf[dim+nzcnt].next;
        uc_inf[i].prev = dim+nzcnt;
        uc_inf[dim+nzcnt].next = i;
        uc_inf[uc_inf[i].next].prev = i;

        nzcnt = ur_inf[i].pivcnt;
        if (nzcnt >= max_k) nzcnt = max_k;
        ur_inf[i].next = ur_inf[dim+nzcnt].next;
        ur_inf[i].prev = dim+nzcnt;
        ur_inf[dim+nzcnt].next = i;
        ur_inf[ur_inf[i].next].prev = i;
    }

#ifdef TRACK_FACTOR
    max = 0;
    nzcnt = 0;
    for (i=0; i<dim; i++) {
        if (ur_inf[i].max > max) max = ur_inf[i].max;
        nzcnt += ur_inf[i].nzcnt;
    }

    f->maxelem_orig = max;
    f->nzcnt_orig   = nzcnt;

    f->maxelem_factor = f->maxelem_orig;
    f->nzcnt_factor = f->nzcnt_orig;
#endif /* TRACK_FACTOR */
    
    /* sentinal for column space */
    ucindx[f->uc_space] = 0;

    clear_work (f);
    
 CLEANUP:
    ILL_RETURN (rval, "init_matrix");
}

static int build_iteration_u_data (factor_work *f)
{
    int       dim    = f->dim;
    ur_info  *ur_inf = f->ur_inf;
    uc_info  *uc_inf = f->uc_inf;
    double   *uccoef = (double *) NULL;
    int      *ucindx = (int *) NULL;
    int      *urindx = f->urindx;
    double   *urcoef = f->urcoef;
    int      *ucrind = (int *) NULL;
    int      *urcind = (int *) NULL;
    int       nzcnt;
    int       beg;
    int       cbeg;
    int       cnzcnt;
    int       uc_space = f->uc_space;
    int       er_space;
    int i;
    int j;
    int k;
    int rval;

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        nzcnt += ur_inf[i].nzcnt;
    }

#ifdef TRACK_FACTOR
    f->nzcnt_factor = nzcnt;
#endif /* TRACK_FACTOR */

    ILL_IFFREE (f->uccoef, double);
    ILL_SAFE_MALLOC (uccoef, nzcnt, double);
    f->uccoef = uccoef;

    ILL_IFFREE (f->ucrind, int);
    ILL_SAFE_MALLOC (ucrind, nzcnt, int);
    f->ucrind = ucrind;

    ILL_IFFREE (f->urcind, int);
    ILL_SAFE_MALLOC (urcind, f->ur_space, int);
    f->urcind = urcind;

    if (uc_space < nzcnt) {
        ILL_IFFREE (f->ucindx, int);
        ILL_SAFE_MALLOC (f->ucindx, nzcnt+1, int);
    }
    f->uc_space = nzcnt;
    uc_space = nzcnt;
    ucindx = f->ucindx;

    for (i=0; i<dim; i++) {
        uc_inf[i].nzcnt = 0;
    }

    for (i=0; i<dim; i++) {
        nzcnt = ur_inf[i].nzcnt;
        beg   = ur_inf[i].rbeg;
        for (j=0; j<nzcnt; j++) {
            uc_inf[urindx[beg+j]].nzcnt++;
        }
        ur_inf[i].delay = 0;
    }

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        uc_inf[i].cbeg = nzcnt;
        nzcnt += uc_inf[i].nzcnt;
        uc_inf[i].nzcnt = 0;
        uc_inf[i].delay = 0;
    }

    f->uc_freebeg = nzcnt;
    for (i=nzcnt; i<uc_space; i++) {
        ucindx[i] = -1;
    }
    ucindx[uc_space] = 0;
    
    for (i=0; i<dim; i++) {
        nzcnt = ur_inf[i].nzcnt;
        beg   = ur_inf[i].rbeg;
        k = urindx[beg];
        cbeg = uc_inf[k].cbeg;
        cnzcnt = uc_inf[k].nzcnt;
        if (cnzcnt != 0) {
            ucindx[cbeg + cnzcnt] = ucindx[cbeg];
            uccoef[cbeg + cnzcnt] = uccoef[cbeg];
            ucrind[cbeg + cnzcnt] = ucrind[cbeg];
            urcind[ur_inf[ucindx[cbeg]].rbeg + ucrind[cbeg]] = cnzcnt;
        }
        ucindx[cbeg] = i;
        uccoef[cbeg] = urcoef[beg];
        ucrind[cbeg] = 0;
        urcind[beg] = 0;
        uc_inf[k].nzcnt = cnzcnt+1;
        for (j=1; j<nzcnt; j++) {
            k = urindx[beg+j];
            cbeg = uc_inf[k].cbeg;
            cnzcnt = uc_inf[k].nzcnt;
            ucindx[cbeg + cnzcnt] = i;
            uccoef[cbeg + cnzcnt] = urcoef[beg+j];
            ucrind[cbeg + cnzcnt] = j;
            urcind[beg+j] = cnzcnt;
            uc_inf[k].nzcnt++;
        }
    }

    for (i=0; i<dim; i++) {
        f->rrank[f->rperm[i]] = i;
    }
    
    nzcnt = f->ur_space;

    for (i=f->ur_freebeg; i<nzcnt; i++) {
        urindx[i] = -1;
    }
    urindx[nzcnt] = 0;

    clear_work (f);
    
    er_space = f->er_space_mul * f->etamax;
    ILL_SAFE_MALLOC (f->er_inf, f->etamax, er_info);
    ILL_SAFE_MALLOC (f->erindx, er_space, int);
    ILL_SAFE_MALLOC (f->ercoef, er_space, double);
    f->etacnt = 0;
    f->er_freebeg = 0;
    f->er_space = er_space;

    rval = 0;

 CLEANUP:
    ILL_RETURN (rval, "build_iteration_u_data");
}

static int build_iteration_l_data (factor_work *f)
{
    int       dim    = f->dim;
    lc_info  *lc_inf = f->lc_inf;
    lr_info  *lr_inf = f->lr_inf;
    double   *lrcoef = (double *) NULL;
    int      *lrindx = (int *) NULL;
    double   *lccoef = f->lccoef;
    int      *lcindx = f->lcindx;
    int       nzcnt;
    int       beg;
    int       rnzcnt;
    int       rbeg;
    int i;
    int j;
    int k;
    int c;
    int rval;

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        nzcnt += lc_inf[i].nzcnt;
        lr_inf[i].nzcnt = 0;
        lr_inf[i].delay = 0;
        lc_inf[lc_inf[i].c].crank = i;
    }

    ILL_IFFREE (f->lrcoef, double);
    if (nzcnt) {
      ILL_SAFE_MALLOC (lrcoef, nzcnt, double);
      f->lrcoef = lrcoef;
    }

    ILL_IFFREE (f->lrindx, int);
    ILL_SAFE_MALLOC (lrindx, nzcnt+1, int);
    f->lrindx = lrindx;

    for (i=0; i<dim; i++) {
        nzcnt = lc_inf[i].nzcnt;
        beg   = lc_inf[i].cbeg;
        lc_inf[i].delay = 0;
        for (j=0; j<nzcnt; j++) {
            lr_inf[lc_inf[lcindx[beg+j]].crank].nzcnt++;
        }
    }

    nzcnt = 0;
    for (i=0; i<dim; i++) {
        lr_inf[i].rbeg = nzcnt;
        nzcnt += lr_inf[i].nzcnt;
        lr_inf[i].nzcnt = 0;
        lr_inf[i].r = lc_inf[i].c;
        lr_inf[lr_inf[i].r].rrank = i;
    }

    for (i=0; i<dim; i++) {
        nzcnt = lc_inf[i].nzcnt;
        beg   = lc_inf[i].cbeg;
        c     = lc_inf[i].c;
        for (j=0; j<nzcnt; j++) {
            k = lc_inf[lcindx[beg+j]].crank;
            rbeg = lr_inf[k].rbeg;
            rnzcnt = lr_inf[k].nzcnt;
            lrindx[rbeg + rnzcnt] = c;
            lrcoef[rbeg + rnzcnt] = lccoef[beg+j];
            lr_inf[k].nzcnt++;
        }
    }

#ifdef TRACK_FACTOR
    nzcnt = f->nzcnt_factor;
    for (i=0; i<dim; i++) {
        nzcnt += lc_inf[i].nzcnt;
    }
    f->nzcnt_factor = nzcnt;

    f->maxelem_cur = f->maxelem_factor;
    f->nzcnt_cur   = f->nzcnt_factor;

/*
    dump_factor_stats (f);
    printf ("orig max  %e nzcnt %d\n", f->maxelem_orig, f->nzcnt_orig);
    printf ("f maxelem %e nzcnt %d\n", f->maxelem_cur, f->nzcnt_cur);
*/
#endif /* TRACK_FACTOR */
    
    rval = 0;

 CLEANUP:
    ILL_RETURN (rval, "build_iteration_l_data");
}

static int handle_singularity (factor_work *f)
{
    int rval = 0;
    int nsing;
    int *singr = (int *) NULL;
    int *singc = (int *) NULL;
    int i;

    if (f->p_nsing == (int *) NULL ||
        f->p_singr == (int **) NULL ||
        f->p_singc == (int **) NULL) {
        fprintf (stderr, "singular basis, but no place for singularity data\n");
        return E_SING_NO_DATA;
    }
    
    nsing = f->nstages - f->stage;
    ILL_SAFE_MALLOC (singr, nsing, int);
    ILL_SAFE_MALLOC (singc, nsing, int);
    for (i=f->stage; i<f->nstages; i++) {
        singr[i - f->stage] = f->rperm[i];
        singc[i - f->stage] = f->cperm[i];
    }
    *f->p_nsing = nsing;
    *f->p_singr = singr;
    *f->p_singc = singc;
    singr = (int *) NULL;
    singc = (int *) NULL;

 CLEANUP:
    ILL_IFFREE (singr, int);
    ILL_IFFREE (singc, int);
    ILL_RETURN (rval, "handle_singularity");
}

static int dense_build_matrix (factor_work *f)
{
    double *dmat = (double *) NULL;
    int stage = f->stage;
    int drows = f->nstages - stage;
    int dcols = f->dim - stage;
    int dsize = drows * dcols;
    int    *crank  = f->crank;
    double *urcoef = f->urcoef;
    int    *urindx = f->urindx;
    int nzcnt;
    int beg;
    int i;
    int r;
    int j;
    int rval = 0;

    ILL_SAFE_MALLOC (dmat, dsize, double);

    for (i=0; i<dsize; i++) dmat[i] = 0.0;

    for (i=0; i<drows; i++) {
        r = f->rperm[i + stage];
        nzcnt = f->ur_inf[r].nzcnt;
        beg   = f->ur_inf[r].rbeg;
        for (j=0; j<nzcnt; j++) {
            dmat[i*dcols - stage + crank[urindx[beg+j]]] = urcoef[beg+j];
        }
    }

    f->drows = drows;
    f->dcols = dcols;
    f->dense_base = f->stage;
    f->dmat = dmat;
    dmat = (double *) NULL;
    
 CLEANUP:
    ILL_IFFREE (dmat, double);
    ILL_RETURN (rval, "dense_build_matrix");
}

static int dense_find_pivot (factor_work *f, int *p_r, int *p_c)
{
    int     dcols   = f->dcols;
    int     drows   = f->drows;
    double *dmat    = f->dmat;
    int     dense_base = f->dense_base;
    int     s       = f->stage - dense_base;
    ur_info *ur_inf = f->ur_inf;
    int     *rperm  = f->rperm;
    double  maxval;
    int     max_r;
    int     max_c;
    int i;

    maxval = 0.0;
    max_r  = -1;
    for (i=s; i<drows; i++) {
        if (ur_inf[rperm[dense_base+i]].max > maxval) {
            maxval = ur_inf[rperm[dense_base+i]].max;
            max_r = i;
        }
    }
    if (max_r == -1) {
        return E_NO_PIVOT;
    }

    maxval = 0.0;
    max_c  = -1;
    for (i=s; i<drows; i++) {
        if (fabs(dmat[max_r*dcols + i]) > maxval) {
            maxval = fabs(dmat[max_r*dcols + i]);
            max_c  = i;
        }
    }
    if (max_c == -1) {
        return E_NO_PIVOT;
    }
    *p_r = max_r;
    *p_c = max_c;
    
    return 0;
}

static void dense_swap (factor_work *f, int r, int c)
{
    int     dcols   = f->dcols;
    int     drows   = f->drows;
    double *dmat    = f->dmat;
    int     dense_base = f->dense_base;
    int     s       = f->stage - dense_base;
    double  v;
    int i;

    if (r != s) {
        ILL_SWAP(f->rperm[dense_base+s], f->rperm[dense_base+r], i);
        f->rrank[f->rperm[dense_base+s]] = dense_base+s;
        f->rrank[f->rperm[dense_base+r]] = dense_base+r;
        for (i=0; i<dcols; i++) {
            ILL_SWAP (dmat[s*dcols+i], dmat[r*dcols+i], v);
        }
    }
    if (c != s) {
        ILL_SWAP(f->cperm[dense_base+s], f->cperm[dense_base+c], i);
        f->crank[f->cperm[dense_base+s]] = dense_base+s;
        f->crank[f->cperm[dense_base+c]] = dense_base+c;
        for (i=0; i<drows; i++) {
            ILL_SWAP (dmat[i*dcols+s], dmat[i*dcols+c], v);
        }
    }
}    

static void dense_elim (factor_work *f, int r, int c)
{
    int     dcols   = f->dcols;
    int     drows   = f->drows;
    double *dmat    = f->dmat;
    int     dense_base = f->dense_base;
    int     s       = f->stage - dense_base;
    ur_info *ur_inf = f->ur_inf;
    int     *rperm  = f->rperm;
    double   fzero_tol = f->fzero_tol;
    double  maxelem_factor = f->maxelem_factor;
    double  pivval;
    double  max;
    double  v;
    double  w;
    int i;
    int j;

    dense_swap (f, r, c);
    f->stage++;
    pivval = 1.0/dmat[s*dcols+s];
    for (i=s+1; i<drows; i++) {
        v = dmat[i*dcols+s];
        if (v != 0.0) {
            v *= pivval;
            if (v > fzero_tol || v < -fzero_tol) {
                dmat[i*dcols+s] = v;
                if (fabs(v) > maxelem_factor) maxelem_factor = fabs(v);
                max = 0.0;
                for (j=s+1; j<drows; j++) {
                    w = dmat[i*dcols+j] - v*dmat[s*dcols+j];
                    dmat[i*dcols+j] = w;
                    if (fabs(w) > max) max = fabs(w);
                }
                for (j=drows; j<dcols; j++) {
                    w = dmat[i*dcols+j] - v*dmat[s*dcols+j];
                    dmat[i*dcols+j] = w;
                }
                ur_inf[rperm[dense_base+i]].max = max;
                if (max > maxelem_factor) maxelem_factor = max;
            } else {
                dmat[i*dcols+s] = 0.0;
            }
        }
    }
    f->maxelem_factor = maxelem_factor;
}

static int dense_replace_row (factor_work *f, int i)
{
    int      dcols      = f->dcols;
    int      dense_base = f->dense_base;
    double  *dmat       = f->dmat + i * dcols;
    double   fzero_tol  = f->fzero_tol;
    ur_info *ur_inf     = f->ur_inf;
    int     *cperm      = f->cperm;
    int      r          = f->rperm[dense_base + i];
    double *urcoef;
    int    *urindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j=i; j<dcols; j++) {
        if (dmat[j] > fzero_tol || dmat[j] < -fzero_tol) {
            nzcnt++;
        }
    }
    if (nzcnt > ur_inf[r].nzcnt) {
        if (ur_inf[r].rbeg + ur_inf[r].nzcnt == f->ur_freebeg) {
            f->ur_freebeg = ur_inf[r].rbeg;
        }
        ur_inf[r].nzcnt = 0;
        if (f->ur_freebeg + nzcnt > f->ur_space) {
            rval = make_ur_space (f, nzcnt);
            ILL_CLEANUP_IF (rval);
        }
        ur_inf[r].rbeg = f->ur_freebeg;
        f->ur_freebeg += nzcnt;
    }
    beg = ur_inf[r].rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    for (j=i; j<dcols; j++) {
        if (dmat[j] > fzero_tol || dmat[j] < -fzero_tol) {
            urcoef[beg] = dmat[j];
            urindx[beg] = cperm[dense_base+j];
            beg++;
        }
    }
    ur_inf[r].nzcnt = beg - ur_inf[r].rbeg;
 CLEANUP:
    ILL_RETURN (rval, "dense_replace_row");
}

static int dense_create_col (factor_work *f, int i)
{
    int      dcols      = f->dcols;
    int      drows      = f->drows;
    int      dense_base = f->dense_base;
    double  *dmat       = f->dmat;
    double   fzero_tol  = f->fzero_tol;
    lc_info *lc_inf     = f->lc_inf;
    int     *rperm      = f->rperm;
    double  *lccoef;
    int     *lcindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j=i+1; j<drows; j++) {
        if (dmat[j*dcols+i] > fzero_tol || dmat[j*dcols+i] < -fzero_tol) {
            nzcnt++;
        }
    }

    if (f->lc_freebeg + nzcnt >= f->lc_space) {
        rval = make_lc_space (f, nzcnt);
        ILL_CLEANUP_IF (rval);
    }
    beg = f->lc_freebeg;
    lc_inf[dense_base+i].cbeg = beg;
    lc_inf[dense_base+i].c = rperm[dense_base+i];
    lcindx = f->lcindx;
    lccoef = f->lccoef;
    
    for (j=i+1; j<drows; j++) {
        if (dmat[j*dcols+i] > fzero_tol || dmat[j*dcols+i] < -fzero_tol) {
            lccoef[beg] = dmat[j*dcols+i];
            lcindx[beg] = rperm[dense_base+j];
            beg++;
        }
    }
    lc_inf[dense_base+i].nzcnt = beg - lc_inf[dense_base+i].cbeg;
    f->lc_freebeg = beg;
 CLEANUP:
    ILL_RETURN (rval, "create_col");
}

static int dense_replace (factor_work *f)
{
    int drows = f->drows;
    int rval = 0;
    int i;

    for (i=0; i<drows; i++) {
        rval = dense_replace_row (f, i);
        ILL_CLEANUP_IF (rval);
        rval = dense_create_col (f, i);
        ILL_CLEANUP_IF (rval);
    }
    ILL_IFFREE (f->dmat, double);
    f->drows = 0;
    f->dcols = 0;
 CLEANUP:
    ILL_RETURN (rval, "dense_replace");
}

static int dense_factor (factor_work *f)
{
    int r;
    int c;
    int rval = 0;

/*
    printf ("dense kernel, %d rows, %d  cols...\n", f->nstages - f->stage,
            f->dim - f->stage);
    fflush (stdout);
*/
    
    rval = dense_build_matrix (f);
    ILL_CLEANUP_IF (rval);

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>1)
    printf ("before Dense ILLfactor\n");
    dump_matrix (f, 1);
    fflush (stdout);
#endif
#endif

    while (f->stage < f->nstages) {
        r = f->stage - f->dense_base;
        rval = dense_find_pivot (f, &r, &c);
        if (rval == E_NO_PIVOT) {
            rval = handle_singularity (f);
            ILL_CLEANUP_IF (rval);
            return E_SINGULAR_INTERNAL;
        } else {
            ILL_CLEANUP_IF (rval);
        }
#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>2)
        printf ("dense pivot elem: %d %d\n", r, c);
        fflush (stdout);
#endif
#endif /* FACTOR_DEBUG */
        dense_elim (f, r, c);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
        if (f->maxelem_factor > f->maxmult * f->maxelem_orig &&
            f->partial_cur < 1.0) {
            return E_FACTOR_BLOWUP;
        }
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>1)
        printf ("After dense pivot stage %d (%d) of %d (%d)\n",
                f->stage - f->dense_base, f->stage,
                f->nstages - f->dense_base, f->nstages);
        fflush (stdout);
#endif
#if (FACTOR_DEBUG+0>2)
        dump_matrix (f, 1);
        fflush (stdout);
#endif
#endif /* FACTOR_DEBUG */
    }

#ifdef FACTOR_DEBUG
    printf ("After dense ILLfactor:\n");
    dump_matrix (f, 0);
    fflush (stdout);
#endif /* FACTOR_DEBUG */

    rval = dense_replace (f);
    ILL_CLEANUP_IF (rval);

#ifdef FACTOR_DEBUG
    printf ("After replacement:\n");
    dump_matrix (f, 0);
    fflush (stdout);
#endif /* FACTOR_DEBUG */

 CLEANUP:
    ILL_RETURN (rval, "dense_factor");
}

#ifdef RECORD
FILE *fsave = (FILE *) NULL;
int fsavecnt = 0;
#endif /* RECORD */

static int ILLfactor_try (factor_work *f, int *basis, int *cbeg, int *clen,
        int *cindx, double *ccoef)
{
    int rval = 0;
    int r;
    int c;

#ifdef RECORD
    {
        int ncol = 0;
        int nzcnt = 0;
        int dim = f->dim;
        int i;
        int j;
        char fnambuf[40];
        
        for (i=0; i<dim; i++) {
            if (basis[i] > ncol) ncol = basis[i];
        }
        ncol++;
        for (i=0; i<ncol; i++) {
            nzcnt += clen[i];
        }
        if (fsave) fclose (fsave);
        sprintf (fnambuf, "prob.mat.%d", fsavecnt);
        fsavecnt++;
        fsave = fopen (fnambuf, "w");
        fprintf (fsave, "%d %d %d\n", f->dim, ncol, nzcnt);
        for (i=0; i<dim; i++) {
            fprintf (fsave, "%d ", basis[i]);
        }
        fprintf (fsave, "\n");
        for (i=0; i<ncol; i++) {
            fprintf (fsave, "%d", clen[i]);
            for (j=0; j<clen[i]; j++) {
                fprintf (fsave, " %d %.16e", cindx[cbeg[i]+j],
                         ccoef[cbeg[i]+j]);
            }
            fprintf (fsave, "\n");
        }
        fprintf (fsave, "\n");
        fflush (fsave);
    }
#endif /* RECORD */

    rval = init_matrix (f, basis, cbeg, clen, cindx, ccoef);
    ILL_CLEANUP_IF (rval);

    f->stage = 0;
    f->nstages = f->dim;

#ifdef FACTOR_DEBUG
    printf ("Initial matrix:\n");
#if (FACTOR_DEBUG+0>1)
    dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif /* FACTOR_DEBUG */
#ifdef FACTOR_STATS
    printf ("Initial matrix: ");
    dump_factor_stats (f);
#endif /* FACTOR_STATS */

    while (f->stage < f->nstages) {
        rval = find_pivot (f, &r, &c);
        if (rval == E_NO_PIVOT) {
            rval = handle_singularity (f);
            ILL_CLEANUP_IF (rval);
            return 0;
        } else {
            ILL_CLEANUP_IF (rval);
        }
        if (f->ur_inf[r].pivcnt > f->dense_fract * (f->nstages - f->stage) &&
            f->uc_inf[c].nzcnt  > f->dense_fract * (f->nstages - f->stage) &&
            f->nstages - f->stage > f->dense_min) {
            rval = dense_factor (f);
            if (rval == E_SINGULAR_INTERNAL) return 0;
            if (rval) return rval;
            break;
        }
#ifdef FACTOR_DEBUG
        printf ("pivot elem: %d %d\n", r, c);
        fflush (stdout);
#endif /* FACTOR_DEBUG */
        rval = elim (f, r, c);
        ILL_CLEANUP_IF (rval);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
        if (f->maxelem_factor > f->maxmult * f->maxelem_orig &&
            f->partial_cur < 1.0) {
            return E_FACTOR_BLOWUP;
        }
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>3)
        printf ("After pivot stage %d of %d\n", f->stage, f->nstages);
        dump_matrix (f, 0);
        fflush (stdout);
#endif
#endif /* FACTOR_DEBUG */
    }

    rval = build_iteration_u_data (f);
    ILL_CLEANUP_IF (rval);

    rval = build_iteration_l_data (f);
    ILL_CLEANUP_IF (rval);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
    if (f->maxelem_factor <= f->minmult * f->maxelem_orig &&
        f->partial_cur > f->partial_tol) {
        if (f->partial_cur > 0.5) {
            f->partial_cur = 0.5;
        } else if (f->partial_cur > 0.25) {
            f->partial_cur = 0.25;
        } else if (f->partial_cur > 0.1) {
            f->partial_cur = 0.1;
        } else {
            f->partial_cur /= 10.0;
        }
        if (f->partial_cur < f->partial_tol) {
            f->partial_cur = f->partial_tol;
        }
/*  Bico - comment out for dist 
        fprintf (stderr, "factor good, lowering partial tolerance to %.2f\n",
                 f->partial_cur);
*/
    }
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
    printf ("Factored matrix:\n");
#if (FACTOR_DEBUG+0>1)
    dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif /* FACTOR_DEBUG */

#ifdef FACTOR_STATS
    printf ("Factored matrix: ");
    dump_factor_stats (f);
#endif /* FACTOR_STATS */
 CLEANUP:
    ILL_RETURN (rval, "factor_try");
}

int ILLfactor (factor_work *f, int *basis, int *cbeg, int *clen, int *cindx,
        double *ccoef, int *p_nsing, int **p_singr, int **p_singc)
{
    int rval;

    f->p_nsing = p_nsing;
    f->p_singr = p_singr;
    f->p_singc = p_singc;

    *p_nsing = 0;
    
 AGAIN:    
    rval = ILLfactor_try (f, basis, cbeg, clen, cindx, ccoef);
    if (rval == E_FACTOR_BLOWUP) {
        if (f->partial_cur < 0.1) {
            f->partial_cur *= 10.0;
        } else if (f->partial_cur < 0.25) {
            f->partial_cur = 0.25;
        } else if (f->partial_cur < 0.5) {
            f->partial_cur = 0.5;
        } else if (f->partial_cur < 1.0) {
            f->partial_cur = 1.0;
        } else {
            return rval;
        }
/* Bico - comment out for dist
        fprintf (stderr, "factor blowup, changing partial tolerance to %.2f\n",
                 f->partial_cur);
*/
        goto AGAIN;
    }
    ILL_RETURN (rval, "ILLfactor");
}

static void ILLfactor_ftranl (factor_work *f, double *a)
{
    int      *lcindx = f->lcindx;
    double   *lccoef = f->lccoef;
    lc_info  *lc_inf = f->lc_inf;
    int       dim    = f->dim;
    int       beg;
    int       nzcnt;
    double v;
    int i;
    int j;
    
    for (i=0; i<dim; i++) {
        v = a[lc_inf[i].c];
        if (v != 0.0) {
            nzcnt = lc_inf[i].nzcnt;
            beg   = lc_inf[i].cbeg;
            for (j=0; j<nzcnt; j++) {
                a[lcindx[beg+j]] -= v * lccoef[beg+j];
            }
        }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_ftran a after l %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", a[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_ftran a after l:");
    for (j=0; j<f->dim; j++) {
      printf (" %.3f", a[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftranl2 (factor_work *f, double *a)
{
    int      *lrindx = f->lrindx;
    double   *lrcoef = f->lrcoef;
    lr_info  *lr_inf = f->lr_inf;
    int       dim    = f->dim;
    int       beg;
    int       nzcnt;
    double v;
    int i;
    int j;
    
    for (i=0; i<dim; i++) {
        nzcnt = lr_inf[i].nzcnt;
        beg   = lr_inf[i].rbeg;
        v = a[lr_inf[i].r];
        for (j=0; j<nzcnt; j++) {
            v -= lrcoef[beg+j] * a[lrindx[beg+j]];
        }
        a[lr_inf[i].r] = v;
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_ftran a after l2 %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", a[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_ftran a after l2:");
    for (j=0; j<f->dim; j++) {
        printf (" %.3f", a[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ftranl3_delay (factor_work *f, int c)
{
    lc_info  *lc_inf = f->lc_inf;
    int       nzcnt;
    int      *indx;
    int       i;

    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx  = f->lcindx + lc_inf[c].cbeg;
    for (i=0; i<nzcnt; i++) {
        c = indx[i];
        if (lc_inf[c].delay++ == 0) {
            ftranl3_delay (f, c);
        }
    }      
}

static void ftranl3_delay2 (factor_work *f, int c)
{
    lc_info  *lc_inf = f->lc_inf;
    int       nzcnt;
    int      *indx;
    int       i;
    int       last;

    do {
        c = lc_inf[c].crank;
        nzcnt = lc_inf[c].nzcnt;
        indx  = f->lcindx + lc_inf[c].cbeg;
        last = -1;
        for (i=0; i<nzcnt; i++) {
            c = indx[i];
            if (lc_inf[c].delay++ == 0) {
                if (last >= 0) {
                    ftranl3_delay2 (f, last);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}

static void ftranl3_process (factor_work *f, int c, svector *x)
{
    lc_info  *lc_inf = f->lc_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;

    v = work[c];
    work[c] = 0.0;
    if (v != 0.0) {
        x->indx[x->nzcnt] = c;
        x->coef[x->nzcnt] = v;
        x->nzcnt++;
    }
    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx  = f->lcindx + lc_inf[c].cbeg;
    coef  = f->lccoef + lc_inf[c].cbeg;
    for (i=0; i<nzcnt; i++) {
        c = indx[i];
        work[c] -= v * coef[i];
        if (--lc_inf[c].delay == 0) {
            ftranl3_process (f, c, x);
        }
    }
}

static void ftranl3_process2 (factor_work *f, int c, svector *x)
{
    lc_info  *lc_inf = f->lc_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;
    int       last;

    do {
        v = work[c];
        work[c] = 0.0;
        if (v != 0.0) {
            x->indx[x->nzcnt] = c;
            x->coef[x->nzcnt] = v;
            x->nzcnt++;
        }
        c = lc_inf[c].crank;
        nzcnt = lc_inf[c].nzcnt;
        indx  = f->lcindx + lc_inf[c].cbeg;
        coef  = f->lccoef + lc_inf[c].cbeg;
        last  = -1;
        for (i=0; i<nzcnt; i++) {
            c = indx[i];
            work[c] -= v * coef[i];
            if (--lc_inf[c].delay == 0) {
                if (last >= 0) {
                    ftranl3_process2 (f, last, x);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}
    
static void ILLfactor_ftranl3 (factor_work *f, svector *a, svector *x)
{
    double   *work   = f->work_coef;
    int       anzcnt = a->nzcnt;
    int      *aindx  = a->indx;
    double   *acoef  = a->coef;
    lc_info  *lc_inf = f->lc_inf;
    int i;

    for (i=0; i<anzcnt; i++) {
        if (lc_inf[aindx[i]].delay++ == 0) {
            ftranl3_delay2 (f, aindx[i]);
        }
        work[aindx[i]] = acoef[i];
    }
    x->nzcnt = 0;
    for (i=0; i<anzcnt; i++) {
        if (--lc_inf[aindx[i]].delay == 0) {
            ftranl3_process2 (f, aindx[i], x);
        }
    }
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_ftran x after l3:");
    for (i=0; i<x->nzcnt; i++) {
        printf (" %.3f*%d", x->coef[i], x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftrane (factor_work *f, double *a)
{
    int      *erindx = f->erindx;
    double   *ercoef = f->ercoef;
    er_info  *er_inf = f->er_inf;
    int       etacnt = f->etacnt;
    int       beg;
    int       nzcnt;
    double v;
    int i;
    int j;
    
    for (i=0; i<etacnt; i++) {
        v = a[er_inf[i].r];
        nzcnt = er_inf[i].nzcnt;
        beg   = er_inf[i].rbeg;
        for (j=0; j<nzcnt; j++) {
            v -= ercoef[beg+j] * a[erindx[beg+j]];
        }
        a[er_inf[i].r] = v;
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_ftran a after eta %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", a[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_ftran a after eta:");
    for (j=0; j<f->dim; j++) {
        printf (" %.3f", a[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftrane2 (factor_work *f, svector *a)
{
    int      *erindx = f->erindx;
    double   *ercoef = f->ercoef;
    er_info  *er_inf = f->er_inf;
    int       etacnt = f->etacnt;
    int       beg;
    int       nzcnt;
    int       anzcnt = a->nzcnt;
    int      *aindx  = a->indx;
    double   *acoef  = a->coef;
    double   *work_coef = f->work_coef;
    int      *work_indx = f->work_indx;
    double v;
    double fzero_tol;
    int i;
    int j;
    int r;

    for (i=0; i<anzcnt; i++) {
        work_coef[aindx[i]] = acoef[i];
        work_indx[aindx[i]] = i + 1;
    }
    for (i=0; i<etacnt; i++) {
        r = er_inf[i].r;
        v = work_coef[r];
        nzcnt = er_inf[i].nzcnt;
        beg   = er_inf[i].rbeg;
        for (j=0; j<nzcnt; j++) {
            v -= ercoef[beg+j] * work_coef[erindx[beg+j]];
        }
        if (v != 0.0) {
            work_coef[r] = v;
            if (work_indx[r] == 0) {
                acoef[anzcnt] = v;
                aindx[anzcnt] = r;
                work_indx[r] = anzcnt+1;
                anzcnt++;
            } else {
                acoef[work_indx[r]-1] = v;
            }
        } else {
            work_coef[r] = 0.0;
            if (work_indx[r]) {
                acoef[work_indx[r]-1] = 0.0;
            }
        }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_ftran a after eta2 %d:", i);
        for (j=0; j<anzcnt; j++) {
            printf (" %.3f*%d", acoef[j], aindx[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
    }
    i = 0;
    fzero_tol = f->fzero_tol;
    while (i < anzcnt) {
        work_coef[aindx[i]] = 0.0;
        work_indx[aindx[i]] = 0;
        if (acoef[i] > fzero_tol || acoef[i] < -fzero_tol) {
            i++;
        } else {
            --anzcnt;
            acoef[i] = acoef[anzcnt];
            aindx[i] = aindx[anzcnt];
        }
    }
    a->nzcnt = anzcnt;

#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_ftran a after eta2:");
    for (j=0; j<anzcnt; j++) {
        printf (" %.3f*%d", acoef[j], aindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftranu (factor_work *f, double *a, svector *x)
{
    int      *ucindx = f->ucindx;
    double   *uccoef = f->uccoef;
    uc_info  *uc_inf = f->uc_inf;
    int      *cperm  = f->cperm;
    int      *rperm  = f->rperm;
    int       dim    = f->dim;
    double    szero_tol = f->szero_tol;
    int       xnzcnt = 0;
    int      *xindx  = x->indx;
    double   *xcoef  = x->coef;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;

    for (i=dim-1; i>=0; i--) {
        if ((v = a[rperm[i]]) != 0.0) {
            j = cperm[i];
            beg   = uc_inf[j].cbeg;
            v     /= uccoef[beg];
            if (v > szero_tol || v < -szero_tol) {
                xindx[xnzcnt] = j;
                xcoef[xnzcnt] = v;
                xnzcnt++;
            }
            nzcnt = uc_inf[j].nzcnt;
            for (j=1; j<nzcnt; j++) {
              a[ucindx[beg+j]] -= v * uccoef[beg+j];
            }
            a[rperm[i]] = 0.0;
        }
    }
    x->nzcnt = xnzcnt;
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_ftran x after u:");
    for (j=0; j<x->nzcnt; j++) {
        printf (" %.3f*%d", x->coef[j], x->indx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftranu2 (factor_work *f, double *a, double *x)
{
    int      *urindx = f->urindx;
    double   *urcoef = f->urcoef;
    ur_info  *ur_inf = f->ur_inf;
    int      *rperm  = f->rperm;
    int       dim    = f->dim;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;

    for (i=dim-1; i>=0; i--) {
        nzcnt = ur_inf[rperm[i]].nzcnt;
        beg   = ur_inf[rperm[i]].rbeg;
        v     = a[rperm[i]];
        for (j=1; j<nzcnt; j++) {
            v -= urcoef[beg+j] * x[urindx[beg+j]];
        }
        x[urindx[beg]] = v / urcoef[beg];
    }
}

static void ftranu3_delay (factor_work *f, int c)
{
    uc_info  *uc_inf = f->uc_inf;
    int       nzcnt;
    int      *indx;
    int       i;

    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx  = f->ucindx + uc_inf[c].cbeg;
    for (i=1; i<nzcnt; i++) {
        c = indx[i];
        if (uc_inf[c].delay++ == 0) {
            ftranu3_delay (f, c);
        }
    }      
}

static void ftranu3_delay2 (factor_work *f, int c)
{
    uc_info  *uc_inf = f->uc_inf;
    int       nzcnt;
    int      *indx;
    int       i;
    int       last;

    do {
        c = f->cperm[f->rrank[c]];
        nzcnt = uc_inf[c].nzcnt;
        indx  = f->ucindx + uc_inf[c].cbeg;
        last  = -1;
        for (i=1; i<nzcnt; i++) {
            c = indx[i];
            if (uc_inf[c].delay++ == 0) {
                if (last >= 0) {
                    ftranu3_delay2 (f, last);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}

static void ftranu3_process (factor_work *f, int c, svector *x)
{
    uc_info  *uc_inf = f->uc_inf;
    double    szero_tol = f->szero_tol;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;

    v = work[c];
    work[c] = 0.0;
    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx  = f->ucindx + uc_inf[c].cbeg;
    coef  = f->uccoef + uc_inf[c].cbeg;
    v /= coef[0];
    if (v > szero_tol || v < -szero_tol) {
        x->indx[x->nzcnt] = c;
        x->coef[x->nzcnt] = v;
        x->nzcnt++;
    }
    for (i=1; i<nzcnt; i++) {
        c = indx[i];
        work[c] -= v * coef[i];
        if (--uc_inf[c].delay == 0) {
            ftranu3_process (f, c, x);
        }
    }
}

static void ftranu3_process2 (factor_work *f, int c, svector *x)
{
    uc_info  *uc_inf = f->uc_inf;
    double    szero_tol = f->szero_tol;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;
    int       last;

    do {
        v = work[c];
        work[c] = 0.0;
        c = f->cperm[f->rrank[c]];
        nzcnt = uc_inf[c].nzcnt;
        indx  = f->ucindx + uc_inf[c].cbeg;
        coef  = f->uccoef + uc_inf[c].cbeg;
        v /= coef[0];
        if (v > szero_tol || v < -szero_tol) {
            x->indx[x->nzcnt] = c;
            x->coef[x->nzcnt] = v;
            x->nzcnt++;
        }
        last  = -1;
        for (i=1; i<nzcnt; i++) {
            c = indx[i];
            work[c] -= v * coef[i];
            if (--uc_inf[c].delay == 0) {
                if (last >= 0) {
                    ftranu3_process2 (f, last, x);
                }
                last = c;
            }
        }
        c = last;
    } while (c >= 0);
}

static void ILLfactor_ftranu3 (factor_work *f, svector *a, svector *x)
{
    double   *work   = f->work_coef;
    int       anzcnt = a->nzcnt;
    int      *aindx  = a->indx;
    double   *acoef  = a->coef;
    uc_info  *uc_inf = f->uc_inf;
    int i;

    for (i=0; i<anzcnt; i++) {
        if (uc_inf[aindx[i]].delay++ == 0) {
            ftranu3_delay2 (f, aindx[i]);
        }
        work[aindx[i]] = acoef[i];
    }
    x->nzcnt = 0;
    for (i=0; i<anzcnt; i++) {
        if (--uc_inf[aindx[i]].delay == 0) {
            ftranu3_process2 (f, aindx[i], x);
        }
    }
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_ftran x after u3:");
    for (i=0; i<x->nzcnt; i++) {
        printf (" %.3f*%d", x->coef[i], x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

/* ILLfactor_ftran solves Bx=a for x */
void ILLfactor_ftran (factor_work *f, svector *a, svector *x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;

#ifdef RECORD
    {
        fprintf (fsave, "f %d", a->nzcnt);
        for (i=0; i<a->nzcnt; i++) {
            fprintf (fsave, " %d %.16e", a->indx[i], a->coef[i]);
        }
        fprintf (fsave, "\n");
        fflush (fsave);
    }
#endif /* RECORD */
#ifdef DEBUG
    {
        printf ("ILLfactor_ftran a:");
        for (i=0; i<a->nzcnt; i++) {
            printf (" %d %.3f", a->indx[i], a->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
        nzcnt = a->nzcnt;
        aindx = a->indx;
        acoef = a->coef;
        for (i=0; i<nzcnt; i++) {
            work_coef[aindx[i]] = acoef[i];
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        ILLfactor_ftranl3 (f, a, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
            nzcnt = f->xtmp.nzcnt;
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;
        
            for (i=0; i<nzcnt; i++) {
                work_coef[aindx[i]] = acoef[i];
            }
            sparse = 0;
        }
    } else {
        ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
        ILLfactor_ftrane2 (f, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
            nzcnt = f->xtmp.nzcnt;
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;
        
            for (i=0; i<nzcnt; i++) {
                work_coef[aindx[i]] = acoef[i];
            }
            sparse = 0;
        }
    } else {
        ILLfactor_ftrane (f, work_coef);
    }

    if (sparse) {
        ILLfactor_ftranu3 (f, &f->xtmp, x);
    } else {
        ILLfactor_ftranu (f, work_coef, x);
    }
    
#ifdef SORT_RESULTS
    sort_vector (x);
#endif
    
#ifdef DEBUG
    {
        printf ("ILLfactor_ftran x:");
        for (i=0; i<x->nzcnt; i++) {
            printf (" %d %.3f", x->indx[i], x->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */
}

/* ILLfactor_ftran_update solves Bx=a for x, and also returns upd, where Ux=upd */
void ILLfactor_ftran_update (factor_work *f, svector *a, svector *upd, svector *x)
{
    int i;
    int nzcnt;
    int dim;
    int sparse;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;
    double szero_tol;

#ifdef RECORD
    {
        fprintf (fsave, "F %d", a->nzcnt);
        for (i=0; i<a->nzcnt; i++) {
            fprintf (fsave, " %d %.16e", a->indx[i], a->coef[i]);
        }
        fprintf (fsave, "\n");
        fflush (fsave);
    }
#endif /* RECORD */
#ifdef DEBUG
    {
        printf ("ILLfactor_ftran_update a:");
        for (i=0; i<a->nzcnt; i++) {
            printf (" %d %.3f", a->indx[i], a->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
        aindx = a->indx;
        acoef = a->coef;
        nzcnt = a->nzcnt;

        for (i=0; i<nzcnt; i++) {
            work_coef[aindx[i]] = acoef[i];
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        ILLfactor_ftranl3 (f, a, upd);
        if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
            nzcnt = upd->nzcnt;
            aindx = upd->indx;
            acoef = upd->coef;
        
            for (i=0; i<nzcnt; i++) {
                work_coef[aindx[i]] = acoef[i];
            }
            sparse = 0;
        }
    } else {
        ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
        ILLfactor_ftrane2 (f, upd);
        if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
            nzcnt = upd->nzcnt;
            aindx = upd->indx;
            acoef = upd->coef;
        
            for (i=0; i<nzcnt; i++) {
                work_coef[aindx[i]] = acoef[i];
            }
            sparse = 0;
        }
    } else {
        ILLfactor_ftrane (f, work_coef);
        nzcnt = 0;
        dim = f->dim;
        szero_tol = f->szero_tol;
        aindx = upd->indx;
        acoef = upd->coef;
        for (i=0; i<dim; i++) {
            if (work_coef[i] != 0.0) {
                if (work_coef[i] > szero_tol || work_coef[i] < -szero_tol) {
                    aindx[nzcnt] = i;
                    acoef[nzcnt] = work_coef[i];
                    nzcnt++;
                }
            }
        }
        upd->nzcnt = nzcnt;
    }        

    if (sparse) {
        ILLfactor_ftranu3 (f, upd, x);
    } else {
        ILLfactor_ftranu (f, work_coef, x);
    }
    
#ifdef SORT_RESULTS
    sort_vector (upd);
    sort_vector (x);
#endif

#ifdef DEBUG
    {
        printf ("ILLfactor_ftran x:");
        for (i=0; i<x->nzcnt; i++) {
            printf (" %d %.3f", x->indx[i], x->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */
}

static void ILLfactor_btranl (factor_work *f, double *x)
{
    int      *lcindx = f->lcindx;
    double   *lccoef = f->lccoef;
    lc_info  *lc_inf = f->lc_inf;
    int       dim    = f->dim;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;

    for (i=dim-1; i>=0; i--) {
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_btran x before l %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", x[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
        v = x[lc_inf[i].c];
        nzcnt = lc_inf[i].nzcnt;
        beg   = lc_inf[i].cbeg;
        for (j=0; j<nzcnt; j++) {
            v -= x[lcindx[beg+j]] * lccoef[beg+j];
        }
        x[lc_inf[i].c] = v;
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_btran x after l:");
    for (j=0; j<f->dim; j++) {
        printf (" %.3f", x[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btranl2 (factor_work *f, double *x)
{
    int      *lrindx = f->lrindx;
    double   *lrcoef = f->lrcoef;
    lr_info  *lr_inf = f->lr_inf;
    int       dim    = f->dim;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;

    for (i=dim-1; i>=0; i--) {
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_btran x before l2 %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", x[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
        v = x[lr_inf[i].r];
        if (v != 0.0) {
            nzcnt = lr_inf[i].nzcnt;
            beg   = lr_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                x[lrindx[beg+j]] -= v * lrcoef[beg+j];
            }
        }
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_btran x after l2:");
    for (j=0; j<f->dim; j++) {
        printf (" %.3f", x[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void btranl3_delay (factor_work *f, int r)
{
    lr_info  *lr_inf = f->lr_inf;
    int       nzcnt;
    int      *indx;
    int       i;

    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx  = f->lrindx + lr_inf[r].rbeg;
    for (i=0; i<nzcnt; i++) {
        r = indx[i];
        if (lr_inf[r].delay++ == 0) {
            btranl3_delay (f, r);
        }
    }      
}

static void btranl3_delay2 (factor_work *f, int r)
{
    lr_info  *lr_inf = f->lr_inf;
    int       nzcnt;
    int      *indx;
    int       i;
    int       last;

    do {
        r = lr_inf[r].rrank;
        nzcnt = lr_inf[r].nzcnt;
        indx  = f->lrindx + lr_inf[r].rbeg;
        last = -1;
        for (i=0; i<nzcnt; i++) {
            r = indx[i];
            if (lr_inf[r].delay++ == 0) {
                if (last >= 0) {
                    btranl3_delay2 (f, last);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static void btranl3_process (factor_work *f, int r, svector *x)
{
    lr_info  *lr_inf = f->lr_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    double    szero_tol = f->szero_tol;
    int       i;
    double    v;

    v = work[r];
    work[r] = 0.0;
    if (v > szero_tol || v < -szero_tol) {
        x->indx[x->nzcnt] = r;
        x->coef[x->nzcnt] = v;
        x->nzcnt++;
    }
    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx  = f->lrindx + lr_inf[r].rbeg;
    coef  = f->lrcoef + lr_inf[r].rbeg;
    for (i=0; i<nzcnt; i++) {
        r = indx[i];
        work[r] -= v * coef[i];
        if (--lr_inf[r].delay == 0) {
            btranl3_process (f, r, x);
        }
    }
}
    
static void btranl3_process2 (factor_work *f, int r, svector *x)
{
    lr_info  *lr_inf = f->lr_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    double    szero_tol = f->szero_tol;
    int       i;
    double    v;
    int       last;

    do {
        v = work[r];
        work[r] = 0.0;
        if (v > szero_tol || v < -szero_tol) {
            x->indx[x->nzcnt] = r;
            x->coef[x->nzcnt] = v;
            x->nzcnt++;
        }
        r = lr_inf[r].rrank;
        nzcnt = lr_inf[r].nzcnt;
        indx  = f->lrindx + lr_inf[r].rbeg;
        coef  = f->lrcoef + lr_inf[r].rbeg;
        last  = -1;
        for (i=0; i<nzcnt; i++) {
            r = indx[i];
            work[r] -= v * coef[i];
            if (--lr_inf[r].delay == 0) {
                if (last >= 0) {
                    btranl3_process2 (f, last, x);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}
    
static void ILLfactor_btranl3 (factor_work *f, svector *a, svector *x)
{
    double   *work   = f->work_coef;
    int       anzcnt = a->nzcnt;
    int      *aindx  = a->indx;
    double   *acoef  = a->coef;
    lr_info  *lr_inf = f->lr_inf;
    int i;

    for (i=0; i<anzcnt; i++) {
        if (lr_inf[aindx[i]].delay++ == 0) {
            btranl3_delay2 (f, aindx[i]);
        }
        work[aindx[i]] = acoef[i];
    }
    x->nzcnt = 0;
    for (i=0; i<anzcnt; i++) {
        if (--lr_inf[aindx[i]].delay == 0) {
            btranl3_process2 (f, aindx[i], x);
        }
    }
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_btran x after l3:");
    for (i=0; i<x->nzcnt; i++) {
        printf (" %.3f*%d", x->coef[i], x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btrane (factor_work *f, double *x)
{
    int      *erindx = f->erindx;
    double   *ercoef = f->ercoef;
    er_info  *er_inf = f->er_inf;
    int       etacnt = f->etacnt;
    int       beg;
    int       nzcnt;
    double v;
    int i;
    int j;
    
    for (i=etacnt-1; i>=0; i--) {
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_btran x before eta %d:", i);
        for (j=0; j<f->dim; j++) {
            printf (" %.3f", x[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
        v = x[er_inf[i].r];
        if (v != 0.0) {
            nzcnt = er_inf[i].nzcnt;
            beg   = er_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                x[erindx[beg+j]] -= v * ercoef[beg+j];
            }
        }
    }
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_btran x after eta:");
    for (j=0; j<f->dim; j++) {
        printf (" %.3f", x[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btrane2 (factor_work *f, svector *x)
{
    int      *erindx = f->erindx;
    double   *ercoef = f->ercoef;
    er_info  *er_inf = f->er_inf;
    int       etacnt = f->etacnt;
    int       beg;
    int       nzcnt;
    int       xnzcnt = x->nzcnt;
    int      *xindx  = x->indx;
    double   *xcoef  = x->coef;
    double   *work_coef = f->work_coef;
    int      *work_indx = f->work_indx;
    double v;
    int i;
    int j;
    
    for (i=0; i<xnzcnt; i++) {
        work_coef[xindx[i]] = xcoef[i];
        work_indx[xindx[i]] = i + 1;
    }
    for (i=etacnt-1; i>=0; i--) {
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
        printf ("ILLfactor_btran x before eta2 %d:", i);
        for (j=0; j<xnzcnt; j++) {
            printf (" %.3f*%d", work_coef[xindx[j]], xindx[j]);
        }
        printf ("\n");
        fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
        v = work_coef[er_inf[i].r];
        if (v != 0.0) {
            nzcnt = er_inf[i].nzcnt;
            beg   = er_inf[i].rbeg;
            for (j=0; j<nzcnt; j++) {
                if (work_indx[erindx[beg+j]] == 0) {
                    work_indx[erindx[beg+j]] = xnzcnt;
                    xindx[xnzcnt++] = erindx[beg+j];
                }
                work_coef[erindx[beg+j]] -= v * ercoef[beg+j];
            }
        }
    }
    
    j = 0;
    while (j < xnzcnt) {
        xcoef[j] = work_coef[xindx[j]];
        work_coef[xindx[j]] = 0.0;
        work_indx[xindx[j]] = 0;
        if (xcoef[j] == 0.0) {
            --xnzcnt;
            xindx[j] = xindx[xnzcnt];
        } else {
            j++;
        }
    }
    x->nzcnt = xnzcnt;
    
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
    printf ("ILLfactor_btran x after eta2:");
    for (j=0; j<xnzcnt; j++) {
        printf (" %.3f*%d", xcoef[j], xindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btranu (factor_work *f, double *a, svector *x)
{
    int      *urindx = f->urindx;
    double   *urcoef = f->urcoef;
    ur_info  *ur_inf = f->ur_inf;
    int      *rperm  = f->rperm;
    int      *cperm  = f->cperm;
    int       dim    = f->dim;
    double    szero_tol = f->szero_tol;
    int       xnzcnt = 0;
    int      *xindx = x->indx;
    double   *xcoef = x->coef;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;
    
    for (i=0; i<dim; i++) {
        if ((v = a[cperm[i]]) != 0.0) {
            j = rperm[i];
            beg   = ur_inf[j].rbeg;
            v     /= urcoef[beg];
            if (v > szero_tol || v < -szero_tol) {
                xindx[xnzcnt] = j;
                xcoef[xnzcnt] = v;
                xnzcnt++;
            }
            nzcnt = ur_inf[j].nzcnt;
            for (j=1; j<nzcnt; j++) {
              a[urindx[beg+j]] -= v * urcoef[beg+j];
            }
            a[cperm[i]] = 0.0;
        }
    }
    x->nzcnt = xnzcnt;
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_btran x after u:");
    for (i=0; i<x->nzcnt; i++) {
        printf (" %.3f*%d", x->coef[i], x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btranu2 (factor_work *f, double *a, double *x)
{
    int      *ucindx = f->ucindx;
    double   *uccoef = f->uccoef;
    uc_info  *uc_inf = f->uc_inf;
    int      *cperm  = f->cperm;
    int       dim    = f->dim;
    int       nzcnt;
    int       beg;
    double    v;
    int i;
    int j;
    
    for (i=0; i<dim; i++) {
        nzcnt = uc_inf[cperm[i]].nzcnt;
        beg   = uc_inf[cperm[i]].cbeg;
        v     = a[cperm[i]];
        for (j=1; j<nzcnt; j++) {
            v -= uccoef[beg+j] * x[ucindx[beg+j]];
        }
        x[ucindx[beg]] = v / uccoef[beg];
    }
}

static void btranu3_delay (factor_work *f, int r)
{
    ur_info  *ur_inf = f->ur_inf;
    int       nzcnt;
    int      *indx;
    int       i;

    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx  = f->urindx + ur_inf[r].rbeg;
    for (i=1; i<nzcnt; i++) {
        r = indx[i];
        if (ur_inf[r].delay++ == 0) {
            btranu3_delay (f, r);
        }
    }      
}

static void btranu3_delay2 (factor_work *f, int r)
{
    ur_info  *ur_inf = f->ur_inf;
    int       nzcnt;
    int      *indx;
    int       i;
    int       last;

    do {
        r = f->rperm[f->crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx  = f->urindx + ur_inf[r].rbeg;
        last  = -1;
        for (i=1; i<nzcnt; i++) {
            r = indx[i];
            if (ur_inf[r].delay++ == 0) {
                if (last >= 0) {
                    btranu3_delay2 (f, last);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static void btranu3_process (factor_work *f, int r, svector *x)
{
    ur_info  *ur_inf = f->ur_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;

    v = work[r];
    work[r] = 0.0;
    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx  = f->urindx + ur_inf[r].rbeg;
    coef  = f->urcoef + ur_inf[r].rbeg;
    v /= coef[0];
    if (v != 0.0) {
        x->indx[x->nzcnt] = r;
        x->coef[x->nzcnt] = v;
        x->nzcnt++;
    }
    for (i=1; i<nzcnt; i++) {
        r = indx[i];
        work[r] -= v * coef[i];
        if (--ur_inf[r].delay == 0) {
            btranu3_process (f, r, x);
        }
    }
}

static void btranu3_process2 (factor_work *f, int r, svector *x)
{
    ur_info  *ur_inf = f->ur_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;
    int       last;

    do {
        v = work[r];
        work[r] = 0.0;
        r = f->rperm[f->crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx  = f->urindx + ur_inf[r].rbeg;
        coef  = f->urcoef + ur_inf[r].rbeg;
        v /= coef[0];
        if (v != 0.0) {
            x->indx[x->nzcnt] = r;
            x->coef[x->nzcnt] = v;
            x->nzcnt++;
        }
        last = -1;
        for (i=1; i<nzcnt; i++) {
            r = indx[i];
            work[r] -= v * coef[i];
            if (--ur_inf[r].delay == 0) {
                if (last >= 0) {
                    btranu3_process2 (f, last, x);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static void ILLfactor_btranu3 (factor_work *f, svector *a, svector *x)
{
    double   *work   = f->work_coef;
    int       anzcnt = a->nzcnt;
    int      *aindx  = a->indx;
    double   *acoef  = a->coef;
    ur_info  *ur_inf = f->ur_inf;
    int i;

    for (i=0; i<anzcnt; i++) {
        if (ur_inf[aindx[i]].delay++ == 0) {
            btranu3_delay2 (f, aindx[i]);
        }
        work[aindx[i]] = acoef[i];
    }
    x->nzcnt = 0;
    for (i=0; i<anzcnt; i++) {
        if (--ur_inf[aindx[i]].delay == 0) {
            btranu3_process2 (f, aindx[i], x);
        }
    }
#ifdef SOLVE_DEBUG
    printf ("ILLfactor_btran x after u3:");
    for (i=0; i<x->nzcnt; i++) {
        printf (" %.3f*%d", x->coef[i], x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif /* SOLVE_DEBUG */
}

/* ILLfactor_btran solves x^tB=a^t (or, B^t x = a) for x */
void ILLfactor_btran (factor_work *f, svector *a, svector *x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx = a->indx;
    double *acoef = a->coef;
    double *work_coef = f->work_coef;
    double szero_tol = f->szero_tol;
    int dim = f->dim;
    
#ifdef RECORD
    {
        fprintf (fsave, "b %d", a->nzcnt);
        for (i=0; i<a->nzcnt; i++) {
            fprintf (fsave, " %d %.16e", a->indx[i], a->coef[i]);
        }
        fprintf (fsave, "\n");
        fflush (fsave);
    }
#endif /* RECORD */
#ifdef DEBUG
    {
        printf ("ILLfactor_btran a:");
        for (i=0; i<a->nzcnt; i++) {
            printf (" %d %.3f", a->indx[i], a->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
        aindx = a->indx;
        acoef = a->coef;
        work_coef = f->work_coef;
        nzcnt = a->nzcnt;
        for (i=0; i<nzcnt; i++) {
            work_coef[aindx[i]] = acoef[i];
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        ILLfactor_btranu3 (f, a, &f->xtmp);
    } else {
        ILLfactor_btranu (f, work_coef, &f->xtmp);
    }

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
        aindx = f->xtmp.indx;
        acoef = f->xtmp.coef;
        work_coef = f->work_coef;
        nzcnt = f->xtmp.nzcnt;
        for (i=0; i<nzcnt; i++) {
            work_coef[aindx[i]] = acoef[i];
        }
        sparse = 0;
    } else {
        sparse = 1;
    }

    if (sparse) {
        ILLfactor_btrane2 (f, &f->xtmp);
        if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
            aindx = f->xtmp.indx;
            acoef = f->xtmp.coef;
            work_coef = f->work_coef;
            nzcnt = f->xtmp.nzcnt;
            for (i=0; i<nzcnt; i++) {
                work_coef[aindx[i]] = acoef[i];
            }
            sparse = 0;
        }
    } else {
        ILLfactor_btrane (f, work_coef);
    }

    if (sparse) {
        ILLfactor_btranl3 (f, &f->xtmp, x);
    } else {
        ILLfactor_btranl2 (f, work_coef);
        dim = f->dim;
        szero_tol = f->szero_tol;
        nzcnt = 0;
        aindx = x->indx;
        acoef = x->coef;
        for (i=0; i<dim; i++) {
            if (work_coef[i] != 0.0) {
                if (work_coef[i] > szero_tol || work_coef[i] < -szero_tol) {
                    aindx[nzcnt] = i;
                    acoef[nzcnt] = work_coef[i];
                    nzcnt++;
                }
                work_coef[i] = 0.0;
            }
        }
        x->nzcnt = nzcnt;
    }
    
#ifdef SORT_RESULTS
    sort_vector (x);
#endif

#ifdef DEBUG
    {
        printf ("ILLfactor_btran x:");
        for (i=0; i<x->nzcnt; i++) {
            printf (" %d %.3f", x->indx[i], x->coef[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */
}

static int expand_col (factor_work *f, int col)
{
    uc_info  *uc_inf     = f->uc_inf + col;
    int       uc_freebeg = f->uc_freebeg;
    int       nzcnt      = uc_inf->nzcnt;
    int       cbeg;
    double   *uccoef;
    int      *ucindx;
    int      *ucrind;
    int i;
    int rval = 0;

    if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
        rval = make_uc_space (f, nzcnt+1);
        ILL_CLEANUP_IF (rval);
        uc_freebeg = f->uc_freebeg;
    }
    cbeg = uc_inf->cbeg;
    uccoef = f->uccoef;
    ucindx = f->ucindx;
    ucrind = f->ucrind;

    for (i=0; i<nzcnt; i++) {
        uccoef[uc_freebeg+i] = uccoef[cbeg+i];
        ucindx[uc_freebeg+i] = ucindx[cbeg+i];
        ucrind[uc_freebeg+i] = ucrind[cbeg+i];
        ucindx[cbeg+i] = -1;
    }
    
    uc_inf->cbeg = uc_freebeg;
    f->uc_freebeg = uc_freebeg + nzcnt;
 CLEANUP:
    ILL_RETURN (rval, "expand_col");
}

static int expand_row (factor_work *f, int row)
{
    ur_info  *ur_inf     = f->ur_inf + row;
    int       ur_freebeg = f->ur_freebeg;
    int       nzcnt      = ur_inf->nzcnt;
    int       rbeg;
    double   *urcoef;
    int      *urindx;
    int      *urcind;
    int i;
    int rval = 0;

    if (ur_freebeg + nzcnt + 1 >= f->ur_space) {
        rval = make_ur_space (f, nzcnt+1);
        ILL_CLEANUP_IF (rval);
        ur_freebeg = f->ur_freebeg;
    }
    rbeg = ur_inf->rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    urcind = f->urcind;

    for (i=0; i<nzcnt; i++) {
        urcoef[ur_freebeg+i] = urcoef[rbeg+i];
        urindx[ur_freebeg+i] = urindx[rbeg+i];
        urcind[ur_freebeg+i] = urcind[rbeg+i];
        urindx[rbeg+i] = -1;
    }
    
    ur_inf->rbeg = ur_freebeg;
    f->ur_freebeg = ur_freebeg + nzcnt;
 CLEANUP:
    ILL_RETURN (rval, "expand_row");
}

static int add_nonzero (factor_work *f, int row, int col, double val)
{
    ur_info  *ur_inf = f->ur_inf + row;
    uc_info  *uc_inf = f->uc_inf + col;
    int       cnzcnt = uc_inf->nzcnt;
    int       rnzcnt = ur_inf->nzcnt;
    int       cloc   = uc_inf->cbeg + cnzcnt;
    int       rloc   = ur_inf->rbeg + rnzcnt;
    int rval = 0;

    if (f->ucindx[cloc] != -1) {
        rval = expand_col (f, col);
        ILL_CLEANUP_IF (rval);
        cloc = uc_inf->cbeg + cnzcnt;
    }
    if (f->urindx[rloc] != -1) {
        rval = expand_row (f, row);
        ILL_CLEANUP_IF (rval);
        rloc = ur_inf->rbeg + rnzcnt;
    }
    f->ucindx[cloc] = row;
    f->uccoef[cloc] = val;
    f->ucrind[cloc] = rnzcnt;
    f->urindx[rloc] = col;
    f->urcoef[rloc] = val;
    f->urcind[rloc] = cnzcnt;

    if (cloc == f->uc_freebeg) f->uc_freebeg++;
    if (rloc == f->ur_freebeg) f->ur_freebeg++;
    
    uc_inf->nzcnt = cnzcnt+1;
    ur_inf->nzcnt = rnzcnt+1;
 CLEANUP:
    ILL_RETURN (rval, "add_nonzero");
}

static void delete_nonzero_row (factor_work *f, int row, int ind)
{
    ur_info  *ur_inf = f->ur_inf;
    double   *urcoef = f->urcoef;
    int      *urindx = f->urindx;
    int      *urcind = f->urcind;
    int      *ucrind = f->ucrind;
    int       rbeg   = ur_inf[row].rbeg;
    int       nzcnt  = ur_inf[row].nzcnt - 1;
    int       cbeg;

    if (ind != nzcnt) {
        urcoef[rbeg+ind] = urcoef[rbeg+nzcnt];
        urindx[rbeg+ind] = urindx[rbeg+nzcnt];
        urcind[rbeg+ind] = urcind[rbeg+nzcnt];
        cbeg = f->uc_inf[urindx[rbeg+nzcnt]].cbeg;
        ucrind[cbeg + urcind[rbeg+nzcnt]] = ind;
        urindx[rbeg+nzcnt] = -1;
    }
    ur_inf[row].nzcnt = nzcnt;
}

static void delete_nonzero_col (factor_work *f, int col, int ind)
{
    uc_info  *uc_inf = f->uc_inf;
    double   *uccoef = f->uccoef;
    int      *ucindx = f->ucindx;
    int      *ucrind = f->ucrind;
    int      *urcind = f->urcind;
    int       cbeg   = uc_inf[col].cbeg;
    int       nzcnt  = uc_inf[col].nzcnt - 1;
    int       rbeg;

    if (ind != nzcnt) {
        uccoef[cbeg+ind] = uccoef[cbeg+nzcnt];
        ucindx[cbeg+ind] = ucindx[cbeg+nzcnt];
        ucrind[cbeg+ind] = ucrind[cbeg+nzcnt];
        rbeg = f->ur_inf[ucindx[cbeg+nzcnt]].rbeg;
        urcind[rbeg + ucrind[cbeg+nzcnt]] = ind;
        ucindx[cbeg+nzcnt] = -1;
    }
    uc_inf[col].nzcnt = nzcnt;
}

static void delete_column (factor_work *f, int col)
{
    uc_info  *uc_inf = f->uc_inf;
    int       beg    = uc_inf[col].cbeg;
    int       nzcnt  = uc_inf[col].nzcnt;
    int      *ucindx = f->ucindx + beg;
    int      *ucrind = f->ucrind + beg;
    int i;

    for (i=0; i<nzcnt; i++) {
        delete_nonzero_row (f, ucindx[i], ucrind[i]);
        ucindx[i] = -1;
    }
    uc_inf[col].nzcnt = 0;

#ifdef TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif
    
#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("delete_column corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */
}

static void delete_row (factor_work *f, int row, svector *x)
{
    ur_info  *ur_inf = f->ur_inf;
    int       beg    = ur_inf[row].rbeg;
    int       nzcnt  = ur_inf[row].nzcnt;
    int      *urindx = f->urindx + beg;
    double   *urcoef = f->urcoef + beg;
    int      *urcind = f->urcind + beg;
    int i;

    for (i=0; i<nzcnt; i++) {
        x->indx[i] = urindx[i];
        x->coef[i] = urcoef[i];
        delete_nonzero_col (f, urindx[i], urcind[i]);
        urindx[i] = -1;
    }
    x->nzcnt = nzcnt;
    ur_inf[row].nzcnt = 0;

#ifdef TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif
    
#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("delete_row corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */
}

static void load_row_noflag (factor_work *f, int row)
{
    double *prow_urcoef = f->urcoef + f->ur_inf[row].rbeg;
    int    *prow_urindx = f->urindx + f->ur_inf[row].rbeg;
    int     prow_nzcnt  = f->ur_inf[row].nzcnt;
    double *work_coef   = f->work_coef;
    int i;
    int j;
    
    for (i=0; i<prow_nzcnt; i++) {
        j = prow_urindx[i];
        work_coef[j] = prow_urcoef[i];
    }
}

static int create_column (factor_work *f, svector *a, int col,
        int *p_last_rank)
{
    int   *rrank    = f->rrank;
    int    nzcnt = a->nzcnt;
    int   *aindx = a->indx;
    double *acoef = a->coef;
#ifdef TRACK_FACTOR
    double max = f->maxelem_cur;
#endif /* TRACK_FACTOR */
    int i;
    int j;
    int rval = 0;
    int last_rank = -1;
    
    last_rank = 0;

    for (i=0; i<nzcnt; i++) {
        rval = add_nonzero (f, aindx[i], col, acoef[i]);
        ILL_CLEANUP_IF (rval);
#ifdef TRACK_FACTOR
        if (fabs(acoef[i]) > max) max = fabs(acoef[i]);
#endif /* TRACK_FACTOR */
        j = rrank[aindx[i]];
        if (j > last_rank) last_rank = j;
    }
    *p_last_rank = last_rank;

#ifdef TRACK_FACTOR
    f->nzcnt_cur += nzcnt;
    f->maxelem_cur = max;
#endif /* TRACK_FACTOR */

#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("create_column corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */
    
 CLEANUP:
    ILL_RETURN (rval, "create_column");
}

static int column_rank (factor_work *f, int col)
{
    int *cperm = f->cperm;
    int  dim   = f->dim;
    int i;

    for (i=0; i<dim; i++) {
        if (cperm[i] == col) {
            return i;
        }
    }
    return 0;
}

static void shift_permutations (factor_work *f, int rank_p, int rank_r)
{
    int *cperm = f->cperm;
    int *crank = f->crank;
    int *rperm = f->rperm;
    int *rrank = f->rrank;
    int  col_p = cperm[rank_p];
    int  row_p = rperm[rank_p];
    int i;
    
    for (i=rank_p; i<rank_r; i++) {
        cperm[i] = cperm[i+1];
        crank[cperm[i]] = i;
        rperm[i] = rperm[i+1];
        rrank[rperm[i]] = i;
    }
    cperm[rank_r] = col_p;
    crank[col_p] = rank_r;
    rperm[rank_r] = row_p;
    rrank[row_p] = rank_r;
}

static int eliminate_row (factor_work *f, int rank_p, int rank_r)
{
    ur_info  *ur_inf     = f->ur_inf;
    int      *rperm      = f->rperm;
    int      *cperm      = f->cperm;
    int      *urindx     = f->urindx;
    double   *urcoef     = f->urcoef;
    int      *erindx     = f->erindx;
    double   *ercoef     = f->ercoef;
    double   *work_coef  = f->work_coef;
    double    fzero_tol  = f->fzero_tol;
    double    pivot_mul;
    int       er_freebeg = f->er_freebeg;
    int       er_space   = f->er_space;
#ifdef TRACK_FACTOR
    double    max        = f->maxelem_cur;
#endif /* TRACK_FACTOR */
    int  beg;
    int  nzcnt;
    int i;
    int j;
    int c;
    int r;

    for (i=rank_p; i<rank_r; i++) {
        c = cperm[i];
        if (work_coef[c] > fzero_tol || work_coef[c] < -fzero_tol) {
            r = rperm[i];
            beg = ur_inf[r].rbeg;
            nzcnt = ur_inf[r].nzcnt;
            pivot_mul = work_coef[c] / urcoef[beg];
            work_coef[c] = 0.0;
            for (j=1; j<nzcnt; j++) {
                work_coef[urindx[beg+j]] -= pivot_mul * urcoef[beg+j];/* 0.85 */
            }
            if (er_freebeg >= er_space) {
                /* fprintf (stderr, "no space in eliminate_row\n"); */
                return E_UPDATE_NOSPACE;
            }
            erindx[er_freebeg] = r;
            ercoef[er_freebeg] = pivot_mul;
#ifdef TRACK_FACTOR
            if (fabs(pivot_mul) > max) max = fabs(pivot_mul);
#endif /* TRACK_FACTOR */
            er_freebeg++;
        } else {
            work_coef[c] = 0.0;
        }
    }
    f->er_freebeg = er_freebeg;
#ifdef TRACK_FACTOR
    f->maxelem_cur = max;
#endif /* TRACK_FACTOR */

    return 0;
}

static int create_row (factor_work *f, double *a, int row, int minrank)
{
    int   *cperm    = f->cperm;
    double fzero_tol = f->fzero_tol;
    int    dim      = f->dim;
#ifdef TRACK_FACTOR
    double max      = f->maxelem_cur;
#endif /* TRACK_FACTOR */
    int i;
    int j;
    int rval = 0;
    
    for (i=minrank; i<dim; i++) {
        if (a[cperm[i]] != 0.0) {
            j = cperm[i];
            if (a[j] > fzero_tol || a[j] < -fzero_tol) {
                rval = add_nonzero (f, row, j, a[j]);
                ILL_CLEANUP_IF (rval);
#ifdef TRACK_FACTOR
                if (fabs(a[j]) > max) max = fabs(a[j]);
#endif /* TRACK_FACTOR */
            }
            a[j] = 0.0;
        }
    }

#ifdef TRACK_FACTOR
    f->nzcnt_cur += f->ur_inf[row].nzcnt;
    f->maxelem_cur = max;
#endif /* TRACK_FACTOR */
    
#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("create_row corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */
 CLEANUP:
    ILL_RETURN (rval, "create_row");
}

static void serow_delay (factor_work *f, int r, int rank_r)
{
    ur_info  *ur_inf = f->ur_inf;
    int      *crank  = f->crank;
    int       nzcnt;
    int      *indx;
    int       i;
    int       last;

    do {
        r = f->rperm[crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx  = f->urindx + ur_inf[r].rbeg;
        last  = -1;
        for (i=1; i<nzcnt; i++) {
            r = indx[i];
            if (ur_inf[r].delay++ == 0 && crank[r] < rank_r) {
                if (last >= 0) {
                    serow_delay (f, last, rank_r);
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
}

static int serow_process (factor_work *f, int r, svector *newr, int rank_r)
{
    ur_info  *ur_inf = f->ur_inf;
    double   *work = f->work_coef;
    int       nzcnt;
    int      *indx;
    double   *coef;
    int       i;
    double    v;
    int       last;
    double    fzero_tol = f->fzero_tol;
    int       rval;

    do {
        v = work[r];
        work[r] = 0.0;
        if (f->crank[r] >= rank_r) {
            if (v > fzero_tol || v < -fzero_tol) {
            /* stash this nonzero in the resulting row */
#ifdef TRACK_FACTOR
                if (fabs(v) > f->maxelem_cur) f->maxelem_cur = fabs(v);
#endif /* TRACK_FACTOR */
                newr->indx[newr->nzcnt] = r;
                newr->coef[newr->nzcnt] = v;
                newr->nzcnt++;
                return 0;
            } else {
                return 0;
            }
        }
        r = f->rperm[f->crank[r]];
        nzcnt = ur_inf[r].nzcnt;
        indx  = f->urindx + ur_inf[r].rbeg;
        coef  = f->urcoef + ur_inf[r].rbeg;
        v /= coef[0];
        if (v > fzero_tol || v < -fzero_tol) {
            /* stash v in eta */
            if (f->er_freebeg >= f->er_space) {
                /* fprintf (stderr, "no space in eliminate_row\n"); */
                return E_UPDATE_NOSPACE;
            }
            f->erindx[f->er_freebeg] = r;
            f->ercoef[f->er_freebeg] = v;
#ifdef TRACK_FACTOR
            if (fabs(v) > f->maxelem_cur) f->maxelem_cur = fabs(v);
#endif /* TRACK_FACTOR */
            f->er_freebeg++;
        }
        last = -1;
        for (i=1; i<nzcnt; i++) {
            r = indx[i];
            work[r] -= v * coef[i];
            if (--ur_inf[r].delay == 0) {
                if (last >= 0) {
                    rval = serow_process (f, last, newr, rank_r);
                    if (rval) return rval;
                }
                last = r;
            }
        }
        r = last;
    } while (r >= 0);
    return 0;
}

static int sparse_eliminate_row (factor_work *f, svector *x, int row_p,
                                 int rank_r)
{
    double   *work   = f->work_coef;
    int       xnzcnt = x->nzcnt;
    int      *xindx  = x->indx;
    double   *xcoef  = x->coef;
    ur_info  *ur_inf = f->ur_inf;
    int      *crank  = f->crank;
    int i;
    int j;
    int rval = 0;
    svector newr;

    newr.indx = (int *) NULL;
    newr.coef = (double *) NULL;
    
    for (i=0; i<xnzcnt; i++) {
        j = xindx[i];
        if (ur_inf[j].delay++ == 0 && crank[j] < rank_r) {
            serow_delay (f, j, rank_r);
        }
        work[j] = xcoef[i];
    }
    
    newr.nzcnt = 0;
    ILL_SAFE_MALLOC (newr.indx, f->dim, int);
    ILL_SAFE_MALLOC (newr.coef, f->dim, double);

    for (i=0; i<xnzcnt; i++) {
        j = xindx[i];
        if (--ur_inf[j].delay == 0) {
            rval = serow_process (f, j, &newr, rank_r);
            ILL_CLEANUP_IF (rval);
        }
    }

    for (i=0; i<newr.nzcnt; i++) {
        rval = add_nonzero (f, row_p, newr.indx[i], newr.coef[i]);
        ILL_CLEANUP_IF (rval);
    }
        
#ifdef TRACK_FACTOR
    f->nzcnt_cur += newr.nzcnt;
#endif /* TRACK_FACTOR */

 CLEANUP:
    ILL_IFFREE (newr.coef, double);
    ILL_IFFREE (newr.indx, int);
    
    /* Bico 031210 - chg from ILL_RETURN */
    ILL_RESULT (rval, "sparse_eliminate_row");
}

static int move_pivot_row (factor_work *f, int r, int c)
{
    ur_info  *ur_inf = f->ur_inf + r;
    uc_info  *uc_inf = f->uc_inf;
    int       beg    = ur_inf->rbeg;
    int       nzcnt  = ur_inf->nzcnt;
    int      *urindx = f->urindx;
    int      *urcind = f->urcind;
    int      *ucrind = f->ucrind;
    double   *urcoef = f->urcoef;
    double dt;
    int i;

    if (urindx[beg] == c) return 0;
    
    for (i=1; i<nzcnt; i++) {
        if (urindx[beg+i] == c) {
            ILL_SWAP (urcoef[beg], urcoef[beg+i], dt);
            ILL_SWAP (urcind[beg], urcind[beg+i], dt);
            urindx[beg+i] = urindx[beg];
            urindx[beg] = c;
            ucrind[uc_inf[c].cbeg + urcind[beg]] = 0;
            ucrind[uc_inf[urindx[beg+i]].cbeg + urcind[beg+i]] = i;
            return 0;
        }
    }
    fprintf (stderr, "pivot row nonzero not found\n");
    return E_UPDATE_SINGULAR_ROW;
}

static int move_pivot_col (factor_work *f, int c, int r)
{
    uc_info  *uc_inf = f->uc_inf + c;
    ur_info  *ur_inf = f->ur_inf;
    int       beg    = uc_inf->cbeg;
    int       nzcnt  = uc_inf->nzcnt;
    int      *ucindx = f->ucindx;
    int      *ucrind = f->ucrind;
    int      *urcind = f->urcind;
    double   *uccoef = f->uccoef;
    double dt;
    int i;

    if (ucindx[beg] == r) return 0;
    
    for (i=1; i<nzcnt; i++) {
        if (ucindx[beg+i] == r) {
            ILL_SWAP (uccoef[beg], uccoef[beg+i], dt);
            ILL_SWAP (ucrind[beg], ucrind[beg+i], dt);
            ucindx[beg+i] = ucindx[beg];
            ucindx[beg] = r;
            urcind[ur_inf[r].rbeg + ucrind[beg]] = 0;
            urcind[ur_inf[ucindx[beg+i]].rbeg + ucrind[beg+i]] = i;
            return 0;
        }
    }
    fprintf (stderr, "pivot col nonzero not found\n");
    return E_UPDATE_SINGULAR_COL;
}

static int move_pivot (factor_work *f, int rank_r)
{
    int r = f->rperm[rank_r];
    int c = f->cperm[rank_r];
    int rval = 0;

    rval = move_pivot_row (f, r, c);
    ILL_CLEANUP_IF (rval);

#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("move_pivot_row corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */

    rval = move_pivot_col (f, c, r);
    ILL_CLEANUP_IF (rval);

#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("move_pivot_col corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */
    
 CLEANUP:
    ILL_RESULT (rval, "move_pivot");  /* Bico 031209 - chg from RETURN */
}

int ILLfactor_update (factor_work *f, svector *a, int col_p, int *p_refact)
{
    int row_p;
    int rank_r;
    int rank_p;
    int rval = 0;
    int nzcnt;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;
    int i;

#ifdef RECORD
    {
        int i;
        fprintf (fsave, "u %d %d", col_p, a->nzcnt);
        for (i=0; i<a->nzcnt; i++) {
            fprintf (fsave, " %d %.16e", a->indx[i], a->coef[i]);
        }
        fprintf (fsave, "\n");
        fflush (fsave);
    }
#endif /* RECORD */

#ifdef DEBUG
    {
        int i;
        printf ("ILLfactor_update col %d:", col_p);
        for (i=0; i<a->nzcnt; i++) {
            printf (" %.3f*%d", a->coef[i], a->indx[i]);
        }
        printf ("\n");
        fflush (stdout);
    }
#endif /* DEBUG */

#ifdef DEBUG
    if (check_matrix (f)) {
        printf ("ILLfactor_update received corrupted matrix\n");
        fflush (stdout);
    }
#endif /* DEBUG */

    if (f->etacnt >= f->etamax) {
        *p_refact = 1;
        return 0;
    }

#ifdef UPDATE_STUDY
    nupdate++;
#endif

    row_p = f->ucindx[f->uc_inf[col_p].cbeg];

    delete_column (f, col_p);

    rval = create_column (f, a, col_p, &rank_r);
    /* if (rval) fprintf (stderr, "create_column failed\n"); */
    ILL_CLEANUP_IF (rval);

    rank_p = f->crank[col_p];
#ifdef UPDATE_STUDY
    if (rank_p != f->rrank[row_p] || rank_p != column_rank (f, col_p)) {
      printf ("rank_p %d rrank[row_p] %d column_rank(f,col_p) %d\n",
              rank_p, f->rrank[row_p], column_rank(f,col_p));
    }
    if (rank_r > rank_p) {
        permshifttot += rank_r - rank_p;
    }
    for (i=0; i<a->nzcnt; i++) {
        if (f->rrank[a->indx[i]] > rank_p) colspiketot++;
    }
    for (i=0; i<f->ur_inf[row_p].nzcnt; i++) {
        if (f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] <= rank_r &&
            f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] != rank_p) {
            rowspiketot++;
        }
    }
#endif

    shift_permutations (f, rank_p, rank_r);

    delete_row (f, row_p, &f->xtmp);

    f->er_inf[f->etacnt].rbeg = f->er_freebeg;
    f->er_inf[f->etacnt].r    = row_p;

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
        nzcnt = f->xtmp.nzcnt;
        aindx = f->xtmp.indx;
        acoef = f->xtmp.coef;

        for (i=0; i<nzcnt; i++) {
            work_coef[aindx[i]] = acoef[i];
        }

        rval = eliminate_row (f, rank_p, rank_r);
        /* if (rval) fprintf (stderr, "eliminate_row failed\n"); */
        ILL_CLEANUP_IF (rval);

        rval = create_row (f, f->work_coef, row_p, rank_r);
        /* if (rval) fprintf (stderr, "create_row failed\n"); */
        ILL_CLEANUP_IF (rval);
    } else {
        rval = sparse_eliminate_row (f, &f->xtmp, row_p, rank_r);
        /* if (rval) fprintf (stderr, "sparse_eliminate_row failed\n"); */
        ILL_CLEANUP_IF (rval);
    }

    if (f->er_freebeg - f->er_inf[f->etacnt].rbeg > 0) {
        f->er_inf[f->etacnt].nzcnt = f->er_freebeg - f->er_inf[f->etacnt].rbeg;
#ifdef TRACK_FACTOR
        f->nzcnt_cur += f->er_inf[f->etacnt].nzcnt;
#endif /* TRACK_FACTOR */
#ifdef UPDATE_STUDY
        leftetatot += f->er_inf[f->etacnt].nzcnt;
#endif

#ifdef SORT_RESULTS
        sort_vector2 (f->er_inf[f->etacnt].nzcnt, f->erindx + f->er_inf[f->etacnt].rbeg, f->ercoef + f->er_inf[f->etacnt].rbeg);
#endif

        f->etacnt++;
    }

    rval = move_pivot (f, rank_r);
    /* if (rval) fprintf (stderr, "move_pivot failed\n"); */
    ILL_CLEANUP_IF (rval);

#ifdef UPDATE_DEBUG
    printf ("Updated factorization:\n");
#if (UPDATE_DEBUG+0>1)
    dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif /* UPDATE_DEBUG */

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
    if (f->maxelem_cur > f->updmaxmult * f->maxelem_orig) {
/* Bico - comment out for dist 
        fprintf (stderr, "factor_update blowup max cur %e max orig %e\n",
                 f->maxelem_cur, f->maxelem_orig);
*/
        return E_FACTOR_BLOWUP;
    }
#endif /* NOTICE_BLOWUP */
#endif
#ifdef UPDATE_STATS
    dump_factor_stats (f);
#endif
 CLEANUP:
    ILL_RESULT (rval, "ILLfactor_update");  /* Bico 031208 - chg from RETURN */
}

