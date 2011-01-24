/* RCS_INFO = "$RCSfile: basis.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iqsutil.h"
#include "lpdefs.h"
#include "qstruct.h"
#include "qsopt.h"
#include "basis.h"
#include "fct.h"
#include "lp.h"
#include "lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define BASIS_STATS 0
#define BASIS_DEBUG 0

static void
    get_var_info (lpinfo *lp, var_data *v);

static int
    init_slack_basis (lpinfo *lp, int *vstat, int *irow, int *rrow,
		      int *unitcol, int *icol, int *rcol),
    get_initial_basis1 (lpinfo *lp, int *vstat),
    get_initial_basis2 (lpinfo *lp, int *vstat),
    set_basis_indices (lpinfo *lp, int *vstat),
    choose_basis (int algorithm, double pinf1, double dinf1,
                  double pinf2, double dinf2);

void ILLbasis_init_basisinfo (lpinfo *lp)
{
   lp->baz    = (int *) NULL;
   lp->nbaz   = (int *) NULL;
   lp->vstat  = (int *) NULL;
   lp->vindex = (int *) NULL;
   lp->f      = (factor_work *) NULL;
}

void ILLbasis_free_basisinfo (lpinfo *lp)
{
   ILL_IFFREE (lp->baz, int);
   ILL_IFFREE (lp->nbaz, int);
   ILL_IFFREE (lp->vstat, int);
   ILL_IFFREE (lp->vindex, int);
   if (lp->f){
      ILLfactor_free_factor_work (lp->f);
      ILL_IFFREE (lp->f, factor_work);
   }
}

int ILLbasis_build_basisinfo (lpinfo *lp)
{
   int rval = 0;

   ILL_SAFE_MALLOC (lp->baz, lp->nrows, int);
   ILL_SAFE_MALLOC (lp->nbaz, lp->nnbasic, int);
   ILL_SAFE_MALLOC (lp->vstat, lp->ncols, int);
   ILL_SAFE_MALLOC (lp->vindex, lp->ncols, int);

   lp->fbasisid = -1;

 CLEANUP:
   if (rval) ILLbasis_free_basisinfo (lp);
   ILL_RETURN (rval, "ILLbasis_build_basisinfo");
}

int ILLbasis_load (lpinfo *lp, ILLlp_basis *B)
{
    int rval = 0;
    char *cstat = B->cstat;
    char *rstat = B->rstat;
    int *structmap = lp->O->structmap;
    int *rowmap = lp->O->rowmap;
    double *rng = lp->O->rangeval;
    int i, j, ncols = lp->ncols, nrows = lp->nrows, nstruct = lp->O->nstruct;
    int basic = 0, nonbasic = 0;

    ILLbasis_free_basisinfo (lp);
    ILLbasis_init_basisinfo (lp);
    rval = ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < nstruct; i++) {
        j = structmap[i];
        if (cstat[i] == QS_COL_BSTAT_BASIC) {
            lp->vstat[j] = STAT_BASIC;
            lp->baz[basic] = j;
            lp->vindex[j] = basic;
            basic++;
        } else {
            lp->nbaz[nonbasic] = j;
            lp->vindex[j] = nonbasic;
            nonbasic++;
            switch (cstat[i]) {
            case QS_COL_BSTAT_LOWER:
                lp->vstat[j] = STAT_LOWER;
                break;
            case QS_COL_BSTAT_UPPER:
                lp->vstat[j] = STAT_UPPER;
                break;
            case QS_COL_BSTAT_FREE:
                lp->vstat[j] = STAT_ZERO;
                break;
            default:
                fprintf (stderr, "unknown col basis stat 1: %c\n", cstat[i]);
                rval = 1; goto CLEANUP;
            }
        }
    }

    for (i = 0; i < nrows; i++) {
        j = rowmap[i];
        if (rng && rng[i] != 0.0) {
            if (rstat[i] == QS_ROW_BSTAT_BASIC) {
                lp->vstat[j] = STAT_BASIC;
                lp->baz[basic] = j;
                lp->vindex[j] = basic;
                basic++;
            } else {
                lp->nbaz[nonbasic] = j;
                lp->vindex[j] = nonbasic;
                nonbasic++;
                switch (rstat[i]) {
                case QS_ROW_BSTAT_LOWER:
                    lp->vstat[j] = STAT_LOWER;
                    break;
                case QS_ROW_BSTAT_UPPER:
                    lp->vstat[j] = STAT_UPPER;
                    break;
                default:
                    fprintf (stderr, "unknown range basis stat 2\n");
                    rval = 1; goto CLEANUP;
                }
            }
        } else {
            switch (rstat[i]) {
            case QS_ROW_BSTAT_BASIC:
                lp->vstat[j] = STAT_BASIC;
                lp->baz[basic] = j;
                lp->vindex[j] = basic;
                basic++;
                break;
            case QS_ROW_BSTAT_LOWER:
                lp->vstat[j] = STAT_LOWER;
                lp->nbaz[nonbasic] = j;
                lp->vindex[j] = nonbasic;
                nonbasic++;
                break;
            default:
                fprintf (stderr, "unknown row basis stat 3\n");
                rval = 1; goto CLEANUP;
            }
        }
    }

    if (basic + nonbasic != ncols) {
        fprintf (stderr, "error in counts in ILLopt_load_basis\n");
        rval = 1; goto CLEANUP;
    }

    if (lp->fbasisid != 0) lp->basisid = 0;
    else                   lp->basisid = 1;

CLEANUP:

    ILL_RETURN (rval, "ILLbasis_load");
}

static void get_var_info (lpinfo *lp, var_data *v)
{
   int  i=0;

   v->nartif   = 0;
   v->nslacks  = 0;
   v->nfree    = 0;
   v->nbndone  = 0;
   v->nbounded = 0;
   v->nfixed   = 0;
   v->cmax     = NINFTY;

   for (i=0; i<lp->ncols; i++){
      switch (lp->vtype[i]){
      case VARTIFICIAL: v->nartif ++; break;
      case VFREE:       v->nfree ++;  break;
      case VLOWER:
      case VUPPER:   if (lp->matcnt[i] == 1 && 
                         lp->O->rowmap[lp->matind[lp->matbeg[i]]] == i)
                        v->nslacks ++;
                     else
                        v->nbndone ++; 
                     break;

      case VFIXED:   v->nfixed ++;
      case VBOUNDED: if (lp->matcnt[i] == 1 && 
                         lp->O->rowmap[lp->matind[lp->matbeg[i]]] == i)
                        v->nslacks ++;
                     else
                        v->nbounded ++; 
                     break;
      }
      if (v->cmax < fabs (lp->cz[i]))  v->cmax = fabs (lp->cz[i]);
   }

#if BASIS_STATS > 0
   printf ("cols = %d, acols = %d, total  = %d, nrows = %d, nlog = %d\n",
           lp->ncols, lp->ncols-lp->nrows,
           v->nartif + v->nfree + v->nslacks + v->nbndone + v->nbounded,
           lp->nrows, v->nartif + v->nslacks);
#endif
}

static int init_slack_basis (lpinfo *lp, int *vstat, int *irow, int *rrow,
			     int *unitcol, int *icol, int *rcol)
{
   int  j, r, vt;
   int  nslacks = 0;

   for (j=0; j<lp->ncols; j++){
      r  = lp->matind[lp->matbeg[j]];
      vt = lp->vtype[j];

      if (( vt == VUPPER || vt == VLOWER || vt == VBOUNDED || vt == VFIXED) &&
          lp->matcnt[j] == 1 && 
          lp->O->rowmap[r] == j){

         vstat[j] = STAT_BASIC;
         irow[r]  = 1;
         rrow[r]  = 1;
         unitcol[r] = j;
         if (icol != NULL){ icol[j] = 1; rcol[j] = 1;}
         nslacks ++;
      }
      else if (vt == VARTIFICIAL){ unitcol[r] = j; vstat[j] = STAT_UPPER;}
      else if (vt == VFREE)                        vstat[j] = STAT_ZERO;
      else if (vt == VFIXED || vt == VUPPER)       vstat[j] = STAT_UPPER;
      else if (vt == VLOWER)                       vstat[j] = STAT_LOWER;
      else if (vt == VBOUNDED){
         if (fabs(lp->lz[j]) < fabs(lp->uz[j]))    vstat[j] = STAT_LOWER;
         else                                      vstat[j] = STAT_UPPER;
      }
   }
   return nslacks;
}

static int primal_col_select (lpinfo *lp, int *vstat, int *irow, int *rrow,
     int *unitcol, double *v, int *perm, int *porder, int nbelem, int pcols)
{
   int i, j, k, tr, r=0;
   int mcnt, mbeg;
   int *matbeg = lp->matbeg;
   int *matcnt = lp->matcnt;
   int *matind = lp->matind;
   double *matval = lp->matval;
   double alpha, val, maxelem;

   for (k=0; k<pcols; k++){
      j = porder[perm[k]];
      mcnt = matcnt[j];
      mbeg = matbeg[j];

      alpha = NINFTY;
      maxelem = NINFTY;

      for (i=0; i<mcnt; i++){
         val = fabs (matval[mbeg+i]);
         if (maxelem < val)  maxelem = val;
         if (rrow[matind[mbeg+i]] == 0 && alpha < val){
            alpha = val;
            r     = matind[mbeg+i];
         }
      }
      if (alpha > PARAM_IBASIS_RPIVOT * maxelem){
         vstat[j] = STAT_BASIC;
         nbelem++;
         irow[r] = 1;
         v[r]    = alpha;
         for (i=0; i<mcnt; i++)
            if (matval[mbeg+i] != 0.0) rrow[matind[mbeg+i]] ++;
      }
      else{
         alpha = NINFTY;
         for (i=0; i<mcnt; i++){
            tr = matind[mbeg+i];
            if (v[tr] != INFTY &&
                fabs (matval[mbeg+i]) > PARAM_IBASIS_RTRIANG * v[tr]){
                alpha = 0.0;
                break;
            }
            if (irow[tr] == 0 && alpha < fabs (matval[mbeg+i])){
               alpha = fabs (matval[mbeg+i]);
               r     = tr;
            }
         }
         if (alpha != 0.0 && alpha != NINFTY){
            vstat[j] = STAT_BASIC;
            nbelem++;
            irow[r] = 1;
            v[r]    = alpha;
            for (i=0; i<mcnt; i++)
               if (matval[mbeg+i] != 0.0) rrow[matind[mbeg+i]] ++;
         }
      }
   }
#if BASIS_STATS > 0
   printf ("nartifs = %d\n", lp->nrows - nbelem);
#endif

   if (nbelem < lp->nrows){
      for (i=0; i<lp->nrows; i++){
         if (irow[i] == 0){
            if (unitcol[i] != -1){
               vstat[unitcol[i]] = STAT_BASIC;
               nbelem++;
            }
            else{
               fprintf (stderr, "Error: Not enough artificials\n");
               return -1;
            }
         }
      }
   }
   return nbelem;
}   

/* This is an implementation of the initial basis procedure
   in: "Implementing the simplex method: the initial basis", by
   Bob Bixby.
   Goals: choose initial variables to go into basis which satisfy:
   1) vars are slacks, 2) vars have freedom to move
   3) initial submatrix is nonsingular, 4) low objective function
   contribution.
*/
static int get_initial_basis1 (lpinfo *lp, int *vstat)
{
   int  rval = 0;
   int  i, j, tot1, tot2;
   int  nbelem = 0, nslacks = 0;
   int  tfree  = 0, tbndone = 0;
   int  tbounded = 0;
   int  *irow = NULL, *rrow = NULL;
   int  *perm = NULL, *porder = NULL;
   int  *unitcol = NULL;
   double  cmax;
   double  *v = NULL;
   double  *qpenalty = NULL;
   var_data  vd;

   get_var_info (lp, &vd);
   cmax = (vd.cmax == 0.0) ? 1.0 : vd.cmax * 1000.0; 

   ILL_SAFE_MALLOC (irow, lp->nrows, int);
   ILL_SAFE_MALLOC (rrow, lp->nrows, int);
   ILL_SAFE_MALLOC (v, lp->nrows, double);
   ILL_SAFE_MALLOC (unitcol, lp->nrows, int);

   for (i=0; i<lp->nrows; i++){
      unitcol[i] = -1;
      v[i] = INFTY;
      irow[i] = 0;
      rrow[i] = 0;
   }

   nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, NULL, NULL);
   if (nslacks != vd.nslacks){
      printf ("complain: incorrect basis info(slacks)\n");
      rval = E_SIMPLEX_ERROR; ILL_CLEANUP;
   }
   if (nslacks == lp->nrows) ILL_CLEANUP;
   nbelem = nslacks;
   if (nbelem < lp->nrows){
      for (i=0; i<lp->nrows; i++){
         if (irow[i] == 0){
            if (unitcol[i] != -1){
               vstat[unitcol[i]] = STAT_BASIC;
               nbelem++;
            }
            else{
               fprintf (stderr, "Error: Not enough artificials\n");
               return -1;
            }
         }
      }
   }
   ILL_CLEANUP;

   tot1 = vd.nfree + vd.nbndone;
   tot2 = vd.nfree + vd.nbndone + vd.nbounded;
   ILL_SAFE_MALLOC (perm,     tot2, int);
   ILL_SAFE_MALLOC (porder,   tot2, int);
   ILL_SAFE_MALLOC (qpenalty, tot2, double);

   for (j=0; j<lp->ncols; j++){
      if (vstat[j] == STAT_BASIC) continue;

      switch (lp->vtype[j]){
      case VFREE:
         porder[tfree]   = j;
         perm[tfree]     = tfree;
         qpenalty[tfree] = lp->cz[j]/cmax;
         tfree++;
         break;

      case VLOWER:
      case VUPPER:
         porder[vd.nfree + tbndone]   = j;
         perm[vd.nfree + tbndone]     = tbndone;
         if (lp->vtype[j] == VLOWER)
            qpenalty[vd.nfree + tbndone] = lp->lz[j] + lp->cz[j]/cmax;
         else
            qpenalty[vd.nfree + tbndone] = -lp->uz[j] + lp->cz[j]/cmax;
         tbndone++;
         break;

      case VFIXED:
      case VBOUNDED:
         porder[tot1 + tbounded]   = j;
         perm[tot1 + tbounded]     = tbounded;
         qpenalty[tot1 + tbounded] = lp->lz[j] - lp->uz[j] + lp->cz[j]/cmax;
         tbounded++;
         break;
      }
   }
   if (tfree != vd.nfree || tbndone != vd.nbndone || tbounded != vd.nbounded){
      printf ("complain: incorrect basis info \n");
      rval = E_SIMPLEX_ERROR; ILL_CLEANUP;
   }

   ILLutil_double_perm_quicksort (perm, qpenalty, vd.nfree);
   ILLutil_double_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree, vd.nbndone);
   ILLutil_double_perm_quicksort (perm + tot1 , qpenalty + tot1, vd.nbounded);

   for (i=0; i<vd.nbndone; i++) perm[vd.nfree+i] += vd.nfree;
   for (i=0; i<vd.nbounded; i++) perm[tot1+i] += tot1;

   nbelem = primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem, tot2);
   if (nbelem != lp->nrows){
      printf ("complain: incorrect final basis size\n"); 
      rval = E_SIMPLEX_ERROR; ILL_CLEANUP;
   }

 CLEANUP:
   if (rval)
      ILLbasis_free_basisinfo (lp);
   ILL_IFFREE (irow, int);
   ILL_IFFREE (rrow, int);
   ILL_IFFREE (v, double);
   ILL_IFFREE (perm, int);
   ILL_IFFREE (porder, int);
   ILL_IFFREE (unitcol, int);
   ILL_IFFREE (qpenalty, double);
   ILL_RETURN (rval, "ILLbasis_get_initial");
}

static int get_initial_basis2 (lpinfo *lp, int *vstat)
{
   int  rval = 0;
   int  i, j, k, tot1, tot2;
   int  rbeg, rcnt, mcnt;
   int  nbelem = 0, nslacks = 0;
   int  tfree  = 0, tbndone = 0;
   int  tbounded = 0;
   int  *irow = NULL, *rrow = NULL;
   int  *perm = NULL, *porder = NULL;
   int  *unitcol = NULL;
   double  cmax;
   double  *v = NULL;
   double  *qpenalty = NULL;
   int  col = 0, s_i = 0, selc = 0;
   int  *icol = NULL, *rcol = NULL;
   int  *plen = NULL;
   double  selv = 0.0, seldj = 0.0;
   double  c_dj = 0.0, *dj = NULL;
   var_data  vd;

   get_var_info (lp, &vd);
   cmax = (vd.cmax == 0.0) ? 1.0 : vd.cmax * 1000.0; 

   ILL_SAFE_MALLOC (irow,   lp->nrows, int);
   ILL_SAFE_MALLOC (rrow,   lp->nrows, int);
   ILL_SAFE_MALLOC (v,      lp->nrows, double);
   ILL_SAFE_MALLOC (unitcol,lp->nrows, int);
   ILL_SAFE_MALLOC (icol,   lp->ncols, int);
   ILL_SAFE_MALLOC (rcol,   lp->ncols, int);
   ILL_SAFE_MALLOC (dj,     lp->ncols, double);

   for (i=0; i<lp->nrows; i++){
      unitcol[i] = -1;
      v[i] = INFTY;
      irow[i] = 0;
      rrow[i] = 0;
   }
   /* assign all d_j */
   for (i=0; i<lp->ncols; i++){
      icol[i] = 0;
      rcol[i] = 0;
      dj[i]   = lp->cz[i];
   }

   nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, icol, rcol);
   if (nslacks != vd.nslacks){
      printf ("complain: incorrect basis info\n");
      rval = E_SIMPLEX_ERROR; ILL_CLEANUP;
   }
   if (nslacks == lp->nrows) ILL_CLEANUP;
   nbelem = nslacks;

   /* allocate maximum required space for perm etc. */
   ILL_SAFE_MALLOC (perm,   lp->ncols, int);
   ILL_SAFE_MALLOC (porder, lp->ncols, int);
   ILL_SAFE_MALLOC (plen,   lp->nrows, int);
   ILL_SAFE_MALLOC (qpenalty, lp->ncols, double);

   /* find all unit rows and record lengths */
   for (i=0; i<lp->nrows; i++){
      if (irow[i] != 1){
         rbeg = lp->rowbeg[i];
         rcnt = lp->rowcnt[i];
         for (j=0; j<rcnt; j++)
            if (lp->rowval[rbeg+j] != 1.0 && lp->rowval[rbeg+j] != -1.0) break;
         if (j == rcnt) {
            perm[s_i]   = s_i;
            porder[s_i] = i;
            plen[s_i] = rcnt;
            s_i ++;
         }
      }
   }

   /*sort all unit rows */
   ILLutil_int_perm_quicksort (perm, plen, s_i);

   /* now go through the unit rows */
   for (k=0; k<s_i; k++){
      i    = porder[perm[k]];
      rbeg = lp->rowbeg[i];
      rcnt = lp->rowcnt[i];
      selc = -1;
      seldj= NINFTY;
      selv = 0.0;
 
      /* for every row s_i, compute min {d_j : d_j <0 , j is u or l or fr}*/
      for (j=0; j<rcnt; j++){
         col = lp->rowind[rbeg+j];
         if (rcol[col] == 1) break;
         if (dj[col] < 0.0){
            if (fabs (dj[col]) > seldj){
               selc  = col;
               seldj = fabs (dj[col]);
               selv  = lp->rowval[rbeg+j];
            }
         }
      }
      /* select pivot element and update all d_j's */
      if (selc != -1){
         nbelem ++;
         irow[i] = 1;
         rrow[i] = 1;
         icol[selc] = 1;
         c_dj = dj[selc]/selv;
         vstat[selc] = STAT_BASIC;

         for (j=0; j<rcnt; j++){
            col = lp->rowind[rbeg+j];
            dj[col] -= lp->rowval[rbeg+j]*c_dj;
            rcol[col] = 1;
         }
      }
   }
#if BASIS_STATS > 0
   printf ("unit rows = %d\n", s_i);
   printf ("nslacks %d, unit rows selected = %d\n", nslacks, nbelem - nslacks);
#endif
   /* now go through remaining cols with dj = 0 */
   tot1 = vd.nfree + vd.nbndone;

   for (j=0; j<lp->ncols; j++){
      if (vstat[j] == STAT_BASIC) continue;
      if (icol[j] == 1 || dj[j] > 1e-7 || dj[j] < -1e-7) continue;
      mcnt = lp->matcnt[j];

      switch (lp->vtype[j]){
      case VFREE:
         porder[tfree]   = j;
         perm[tfree]     = tfree;
         qpenalty[tfree] = lp->cz[j]/cmax + mcnt;
         tfree++;
         break;

      case VLOWER:
      case VUPPER:
         porder[vd.nfree + tbndone]   = j;
         perm[vd.nfree + tbndone]     = tbndone;
         if (lp->vtype[j] == VLOWER)
            qpenalty[vd.nfree + tbndone] = lp->lz[j] + lp->cz[j]/cmax + mcnt;
         else
            qpenalty[vd.nfree + tbndone] = -lp->uz[j] + lp->cz[j]/cmax + mcnt;
         tbndone++;
         break;

      case VFIXED:
      case VBOUNDED:
         porder[tot1 + tbounded]   = j;
         perm[tot1 + tbounded]     = tbounded;
         qpenalty[tot1 + tbounded] = lp->lz[j] - lp->uz[j] + lp->cz[j]/cmax + mcnt;
         tbounded++;
         break;
      }
   }
#if BASIS_STATS > 0
   printf ("bfree %d, bone %d, bbnd %d\n", tfree, tbndone, tbounded);
#endif

   ILLutil_double_perm_quicksort (perm, qpenalty, tfree);
   ILLutil_double_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree, tbndone);
   ILLutil_double_perm_quicksort (perm + tot1 , qpenalty + tot1, tbounded);

   tot2 = tfree+tbndone;
   for (i=0; i<tbndone; i++){
      perm[tfree+i] = perm[vd.nfree+i] + tfree;
      porder[tfree+i] = porder[vd.nfree+i];
   }
   for (i=0; i<tbounded; i++){
      perm[tot2+i] = perm[tot1+i] + tot2;
      porder[tot2+i] = porder[tot1+i];
   }
   tot2 += tbounded;

   nbelem = primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem, tot2);
   if (nbelem != lp->nrows){
      printf ("complain: incorrect final basis size\n"); 
      rval = E_SIMPLEX_ERROR; ILL_CLEANUP;
   }

 CLEANUP:
   if (rval)
      ILLbasis_free_basisinfo (lp);

   ILL_IFFREE (irow, int);
   ILL_IFFREE (rrow, int);
   ILL_IFFREE (v, double);
   ILL_IFFREE (unitcol, int);

   ILL_IFFREE (icol, int);
   ILL_IFFREE (rcol, int);
   ILL_IFFREE (dj, double);
   ILL_IFFREE (perm, int);
   ILL_IFFREE (porder, int);
   ILL_IFFREE (plen, int);
   ILL_IFFREE (qpenalty, double);
   ILL_RETURN (rval, "ILLbasis_get_initial");
}

static int set_basis_indices (lpinfo *lp, int *vstat)
{
   int  i, b=0, nb=0;
   int  vs;

   for (i=0; i<lp->ncols; i++){
      vs = vstat[i];
      lp->vstat[i] = vs;

      if (vs == STAT_BASIC){
         lp->baz[b] = i;
         lp->vindex[i] = b;
         b++;
      }
      else if (vs == STAT_UPPER || vs == STAT_LOWER || vs == STAT_ZERO){
         lp->nbaz[nb] = i;
         lp->vindex[i] = nb;
         nb ++;
      }
      else{
         fprintf (stderr, "Error in basis creation\n");
         return E_SIMPLEX_ERROR;
      }
   }
   if (b != lp->nrows){
      fprintf (stderr, "Error 2 in basis creation\n");
      return E_SIMPLEX_ERROR;
   }
   else if (nb != lp->nnbasic){
      fprintf (stderr, "Error 3 in basis creation\n");
      return E_SIMPLEX_ERROR;
   }
   return 0;
}

int ILLbasis_get_initial (lpinfo *lp, int algorithm)
{
   int  rval = 0;
   int  *vstat = NULL;

   ILLbasis_free_basisinfo (lp);
   ILLbasis_init_basisinfo (lp);
   rval = ILLbasis_build_basisinfo (lp);
   ILL_CLEANUP_IF (rval);

   ILL_SAFE_MALLOC (vstat, lp->ncols, int);

   if (algorithm == PRIMAL_SIMPLEX)
      rval = get_initial_basis1 (lp, vstat);
   else
      rval = get_initial_basis2 (lp, vstat);

   if (rval == E_SIMPLEX_ERROR){
      FILE *f = fopen ("bad.lp", "w");
      int tval = ILLwrite_lp_file(lp->O, f, NULL); 
	  if (tval) {
		  fprintf(stderr, "Error writing bad lp\n"); 
	  } 
      if (f != NULL) fclose (f);
   }
   ILL_CLEANUP_IF (rval);

   rval = set_basis_indices (lp, vstat);
   lp->basisid = 0;

 CLEANUP:
   ILL_IFFREE (vstat, int);
   ILL_RETURN (rval, "ILLbasis_get_initial");   
}

static int choose_basis (int algorithm, double pinf1, double dinf1,
                         double pinf2, double dinf2)
{
   int choice = 1;
   double  rp = 0.0, rd = 0.0;
   double  eps = 0.001;
#define PRI_RLIMIT  0.25
#define INF_RATIO 10.0

   if (algorithm == PRIMAL_SIMPLEX){
      if (pinf1 <= pinf2+eps && dinf1 <= dinf2+eps)
         choice = 1;
      else if (pinf2 <= pinf1+eps && dinf2 <= dinf1+eps)
         choice = 2;
      else if (pinf1 < pinf2 && dinf1 > dinf2){
         choice = 1;
         rp = pinf1 / pinf2;
         rd = dinf2 / dinf1;
         if (rp > PRI_RLIMIT && INF_RATIO*rd < rp) choice = 2;
      }
      else if (pinf1 > pinf2 && dinf1 < dinf2){
         choice = 2;
         rp = pinf2 / pinf1;
         rd = dinf1 / dinf2;
         if (rp > PRI_RLIMIT && INF_RATIO*rd < rp) choice = 1;
      }
      else
         choice = 1;
   }
   return choice;
}
   
int ILLbasis_get_cinitial (lpinfo *lp, int algorithm)
{
   int  rval = 0;
   int  *vstat1 = NULL;
   int  *vstat2 = NULL;
   int  singular;
   int  choice = 0;
#if BASIS_STATS > 0
   int  i, nz1 = 0, nz2 = 0;
#endif
   double  pinf1, pinf2;
   double  dinf1, dinf2;
   feas_info  fi;

   ILLbasis_free_basisinfo (lp);
   ILLbasis_init_basisinfo (lp);
   rval = ILLbasis_build_basisinfo (lp);
   ILL_CLEANUP_IF (rval);

   ILL_SAFE_MALLOC (vstat1, lp->ncols, int);
   ILL_SAFE_MALLOC (vstat2, lp->ncols, int);

   if (algorithm != PRIMAL_SIMPLEX){
      rval = get_initial_basis2 (lp, vstat2);
      ILL_CLEANUP_IF (rval);
      rval = set_basis_indices (lp, vstat2);
      lp->basisid = 0;
      ILL_CLEANUP;
   }

   rval = get_initial_basis1 (lp, vstat1);
   ILL_CLEANUP_IF (rval);
   rval = get_initial_basis2 (lp, vstat2);
   ILL_CLEANUP_IF (rval);
   lp->basisid = 0;

   /* handle first basis */
   rval = set_basis_indices (lp, vstat1);
   ILL_CLEANUP_IF (rval);
#if BASIS_STATS > 0
   for (i=0; i<lp->nrows; i++) nz1 += lp->matcnt[lp->baz[i]];
#endif
   rval = ILLbasis_factor (lp, &singular);
   ILL_CLEANUP_IF (rval);

   ILLfct_compute_piz (lp);
   ILLfct_compute_dz (lp);
   ILLfct_dual_adjust (lp, 0.0);
   ILLfct_compute_xbz (lp);

   ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
   ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
   pinf1 = lp->pinfeas; 
   dinf1 = lp->dinfeas;
   /*
   ILLfct_compute_pobj (lp);  obj1p = lp->objval;
   ILLfct_compute_dobj (lp);  obj1d = lp->objval;
   */

   /* handle second basis */
   rval = set_basis_indices (lp, vstat2);
   ILL_CLEANUP_IF (rval);
#if BASIS_STATS > 0
   for (i=0; i<lp->nrows; i++) nz2 += lp->matcnt[lp->baz[i]];
#endif
   rval = ILLbasis_factor (lp, &singular);
   ILL_CLEANUP_IF (rval);

   ILLfct_compute_piz (lp);
   ILLfct_compute_dz (lp);
   ILLfct_dual_adjust (lp, 0.0);
   ILLfct_compute_xbz (lp);

   ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
   ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
   pinf2 = lp->pinfeas;
   dinf2 = lp->dinfeas;

#if BASIS_STATS > 0
   printf ("b1: nz %d pinf %.2f dinf %.2f\n", nz1, pinf1, dinf1);
   printf ("b2: nz %d pinf %.2f dinf %.2f\n", nz2, pinf2, dinf2);
#endif
   choice = choose_basis (algorithm, pinf1, dinf1, pinf2, dinf2);
   if (choice == 1){
      lp->fbasisid = -1;
      rval = set_basis_indices (lp, vstat1);
      ILL_CLEANUP_IF (rval);
   }

 CLEANUP:
   if (rval == E_SIMPLEX_ERROR){
      FILE *fil = fopen ("bad.lp", "w");
      int tval = ILLwrite_lp_file(lp->O, fil, NULL);
	  if (tval) {
		  fprintf(stderr, "Error writing bad lp\n"); 
	  } 
      if (fil != NULL) fclose (fil);
   }
   ILL_IFFREE (vstat1, int);
   ILL_IFFREE (vstat2, int);
   ILL_RETURN (rval, "ILLbasis_get_initial");   
}

int ILLbasis_factor (lpinfo *lp, int *singular)
{
   int  rval = 0;
   int  i;
   int  eindex;
   int  lindex;
   int  ltype;
   int  lvstat;
   int  nsing  = 0;
   int  *singr = (int *) NULL;
   int  *singc = (int *) NULL;

   *singular = 0;
   do{
      if (lp->f){
         ILLfactor_free_factor_work (lp->f);
      }
      else {
         ILL_SAFE_MALLOC (lp->f, 1, factor_work);
         ILLfactor_init_factor_work (lp->f);
      }
      rval = ILLfactor_create_factor_work (lp->f, lp->nrows);
      ILL_CLEANUP_IF (rval);

      rval = ILLfactor (lp->f, lp->baz, lp->matbeg, lp->matcnt,
                     lp->matind, lp->matval, &nsing, &singr, &singc);
      ILL_CLEANUP_IF (rval);

      if (nsing != 0){
         *singular  = 1;
         for (i=0; i<nsing; i++){
            eindex = lp->vindex[lp->O->rowmap[singr[i]]];
            lindex = singc[i];
            ltype  = lp->vtype[lp->baz[lindex]];

            if (ltype == VBOUNDED || ltype == VLOWER || ltype == VARTIFICIAL)
               lvstat = STAT_LOWER;
            else if (ltype == VUPPER)
               lvstat = STAT_UPPER;
            else
               lvstat = STAT_ZERO;

            ILLfct_update_basis_info (lp, eindex, lindex, lvstat);
            lp->basisid ++;
         }
         ILL_IFFREE (singr, int);
         ILL_IFFREE (singc, int);
      }

   } while (nsing != 0);

   lp->fbasisid = lp->basisid;

 CLEANUP:
   ILL_IFFREE (singr, int);
   ILL_IFFREE (singc, int);
   ILL_RETURN (rval, "ILLbasis_factor");
}

int ILLbasis_refactor (lpinfo *lp)
{
   int  sing = 0;
   int  rval = 0;

   rval = ILLbasis_factor (lp, &sing);
   if (sing) {
      fprintf(stderr, "Singular basis in ILLbasis_refactor()\n");
      rval = -1;
   }

   ILL_RETURN (rval, "ILLbasis_refactor");
}

void ILLbasis_column_solve (lpinfo *lp, svector *rhs, svector *soln)
{
    ILLfactor_ftran (lp->f, rhs, soln);
}

void ILLbasis_column_solve_update (lpinfo *lp, svector *rhs, svector *upd, svector *soln)
{
    ILLfactor_ftran_update (lp->f, rhs, upd, soln);
}

void ILLbasis_row_solve (lpinfo *lp, svector *rhs, svector *soln)
{
    ILLfactor_btran (lp->f, rhs, soln);
}

int ILLbasis_update (lpinfo *lp, svector *y, int lindex, int *refactor, int *singular)
{
#if 0 /* To always refactor, change 0 to 1 */
   *refactor = 1;
   return ILLbasis_factor (lp, singular);
#else

   int  rval = 0;

   *refactor = 0;
   rval = ILLfactor_update (lp->f, y, lindex, refactor);
   if (rval == E_FACTOR_BLOWUP ||  rval == E_UPDATE_SINGULAR_ROW ||  rval == E_UPDATE_SINGULAR_COL) {
/* Bico - comment out for dist
       fprintf(stderr, "Warning: numerically bad basis in ILLfactor_update\n");
*/
       *refactor = 1;
       rval = 0;
   }
   if (rval == E_UPDATE_NOSPACE) {
     *refactor = 1;
     rval = 0;
   }

   if (*refactor) rval = ILLbasis_factor (lp, singular);

   if (rval) {
       FILE *eout = (FILE *) NULL;
       int tval;

       printf ("write bad lp to factor.lp\n"); fflush (stdout);
       eout = fopen ("factor.lp", "w");
       if (!eout) {
           fprintf (stderr, "could not open file to write bad factor lp\n");
       } else {
           tval = ILLwrite_lp_file(lp->O, eout, NULL); 
           if (tval) {
               fprintf (stderr, "error while writing bad factor lp\n");
           }
		   fclose (eout);
       }

       printf ("write bad basis to factor.bas\n"); fflush (stdout);
       tval = ILLlib_writebasis (lp,  (ILLlp_basis *) NULL, "factor.bas");
       if (tval) {
           fprintf (stderr, "error while writing factor basis\n");
       }
   }

   ILL_RETURN (rval, "ILLbasis_update");
#endif
}

