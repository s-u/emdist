/* RCS_INFO = "$RCSfile: fct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#define FCT_DEBUG 0

#include "iqsutil.h"
#include "lpdefs.h"
#include "stddefs.h"
#include "basis.h"
#include "fct.h"
#include "price.h"
#include "ratio.h"
#include "dstruct.h"

static int
   compute_zA1(lpinfo *lp, svector *z, svector *zA),
   compute_zA2(lpinfo *lp, svector *z, svector *zA),
   compute_zA3(lpinfo *lp, svector *z, svector *zA),
   expand_var_bounds (lpinfo *lp, double ftol, int *chgb),
   expand_var_coefs (lpinfo *lp, double ftol, int *chgc);

static void
   update_piv_values (count_struct *c, int phase, double piv),
   copy_vectors (svector *a, svector *b),
   add_vectors (lpinfo *lp, svector *a, svector *b, svector *c, double t);

static double
   my_rand (int bound, ILLrandstate *r);


void ILLfct_load_workvector (lpinfo *lp, svector *s)
{
   int  i;

   for (i=0; i<s->nzcnt; i++){
      lp->work.indx[i] = s->indx[i];
      lp->work.coef[s->indx[i]] = s->coef[i];
   }
   lp->work.nzcnt = s->nzcnt;
}

void ILLfct_zero_workvector (lpinfo *lp)
{
   int  i;

   for (i=0; i<lp->work.nzcnt; i++)
      lp->work.coef[lp->work.indx[i]] = 0.0;
   lp->work.nzcnt = 0;
}

void ILLfct_set_variable_type (lpinfo *lp)
{
   int  j;

   for (j=0; j<lp->ncols; j++){

      if (lp->uz[j] != INFTY && lp->lz[j] != NINFTY){

         if (lp->lz[j] < lp->uz[j])
            lp->vtype[j] = VBOUNDED;
         else if (lp->lz[j] == 0.0 &&
                  lp->matcnt[j] == 1 &&
                  lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
            lp->vtype[j] = VARTIFICIAL;
         else
            lp->vtype[j] = VFIXED;
      }
      else if (lp->uz[j] == INFTY && lp->lz[j] == NINFTY)
         lp->vtype[j] = VFREE;
      else if (lp->uz[j] == INFTY)
         lp->vtype[j] = VLOWER;
      else
         lp->vtype[j] = VUPPER;
   }
}

/* compute various vectors */

void ILLfct_compute_pobj (lpinfo *lp)
{
   int  i, j;
   int  col;
   double  sum = 0.0;

   for (i=0; i<lp->nrows; i++)
      sum += lp->cz[lp->baz[i]] * lp->xbz[i];

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];
      if (lp->vstat[col] == STAT_UPPER)
         sum += lp->cz[col] * lp->uz[col];
      else if (lp->vstat[col] == STAT_LOWER)
         sum += lp->cz[col] * lp->lz[col];
   }
   lp->pobjval = sum;
   lp->objval  = sum;
}

void ILLfct_compute_dobj (lpinfo *lp)
{
   int  i, j;
   int  col;
   double  sum = 0.0;

   for (i=0; i<lp->nrows; i++)
      sum += lp->piz[i] * lp->bz[i];

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];
      if (lp->vstat[col] == STAT_UPPER)
         sum += lp->dz[j] * lp->uz[col];
      else if (lp->vstat[col] == STAT_LOWER)
         sum += lp->dz[j] * lp->lz[col];
   }
   lp->dobjval = sum;
   lp->objval  = sum;
}

void ILLfct_compute_xbz (lpinfo *lp)
{
   int  i, j, r;
   int  col, mcnt, mbeg;
   double  xval;
   svector  *srhs  = &(lp->srhs);
   svector  *ssoln = &(lp->ssoln);

   for (i=0; i<lp->nrows; i++){
      lp->xbz[i]   = 0.0;
      srhs->coef[i] = lp->bz[i];
   }
   for (j=0; j<lp->nnbasic; j++){
      col  = lp->nbaz[j];
      xval = 0.0;
      if (lp->vstat[col] == STAT_UPPER && lp->uz[col] != 0.0)
         xval = lp->uz[col];
      else if (lp->vstat[col] == STAT_LOWER && lp->lz[col] != 0.0)
         xval = lp->lz[col];

      if (xval != 0.0){
         mcnt = lp->matcnt[col];
         mbeg = lp->matbeg[col];
         for (i=0; i<mcnt; i++)
            srhs->coef[lp->matind[mbeg+i]] -= xval*lp->matval[mbeg+i];
      }
   }
   for (i=0, r=0; i<lp->nrows; i++)
      if (srhs->coef[i] != 0.0){
         srhs->coef[r] = srhs->coef[i];
         srhs->indx[r] = i;
         r ++;
      }
   srhs->nzcnt = r;

   ILLbasis_column_solve (lp, srhs, ssoln);
   for (i=0; i<ssoln->nzcnt; i++)
      lp->xbz[ssoln->indx[i]] = ssoln->coef[i];
}

void ILLfct_compute_piz (lpinfo *lp)
{
   int  i, r;
   svector  *srhs  = &(lp->srhs);
   svector  *ssoln = &(lp->ssoln);

   for (i=0, r=0; i<lp->nrows; i++){
      lp->piz[i] = 0.0;
      if (lp->cz[lp->baz[i]] != 0.0){
         srhs->indx[r] = i;
         srhs->coef[r] = lp->cz[lp->baz[i]];
         r++;
      }
   }
   srhs->nzcnt = r;

   ILLbasis_row_solve (lp, srhs, ssoln);
   for (i=0; i<ssoln->nzcnt; i++)
      lp->piz[ssoln->indx[i]] = ssoln->coef[i];
}

void ILLfct_compute_dz (lpinfo *lp)
{
   int  i, j;
   int  col;
   int  mcnt, mbeg;
   double  sum;

   for (j=0; j<lp->nnbasic; j++){
      sum  = 0.0;
      col  = lp->nbaz[j];
      mcnt = lp->matcnt[col];
      mbeg = lp->matbeg[col];
      for (i=0; i<mcnt; i++)
         sum += lp->piz[lp->matind[mbeg+i]] * lp->matval[mbeg+i];
      lp->dz[j] = lp->cz[col] - sum;
   }
}

void ILLfct_compute_phaseI_xbz (lpinfo *lp)
{
   int  i, j, r;
   int  col, mcnt, mbeg;
   svector  *srhs  = &(lp->srhs);
   svector  *ssoln = &(lp->ssoln);

   for (i=0; i<lp->nrows; i++){
      lp->xbz[i]   = 0.0;
      srhs->coef[i] = 0.0;
   }
   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->dfeas[j]){
         mcnt = lp->matcnt[col];
         mbeg = lp->matbeg[col];
         if (lp->dfeas[j] == -1)
            for (i=0; i<mcnt; i++)
               srhs->coef[lp->matind[mbeg+i]] -= lp->matval[mbeg+i];
         else
            for (i=0; i<mcnt; i++)
               srhs->coef[lp->matind[mbeg+i]] += lp->matval[mbeg+i];
      }
   }
   for (i=0, r=0; i<lp->nrows; i++)
      if (srhs->coef[i] != 0.0){
         srhs->coef[r] = srhs->coef[i];
         srhs->indx[r] = i;
         r ++;
      }
   srhs->nzcnt = r;

   ILLbasis_column_solve (lp, srhs, ssoln);
   for (i=0; i<ssoln->nzcnt; i++)
      lp->xbz[ssoln->indx[i]] = ssoln->coef[i];
}

void ILLfct_compute_phaseI_piz (lpinfo *lp)
{
   int  i, r;
   svector  *srhs  = &(lp->srhs);
   svector  *ssoln = &(lp->ssoln);

   for (i=0, r=0; i<lp->nrows; i++){
      lp->pIpiz[i] = 0.0;
      if (lp->bfeas[i] != 0){
         srhs->indx[r] = i;
         srhs->coef[r] = (double) lp->bfeas[i];
         r++;
      }
   }
   srhs->nzcnt = r;

   ILLbasis_row_solve (lp, srhs, ssoln);
   for (i=0; i<ssoln->nzcnt; i++)
      lp->pIpiz[ssoln->indx[i]] = ssoln->coef[i];
   ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, 0.0);
}

void ILLfct_compute_phaseI_dz (lpinfo *lp)
{
   int  i, j;
   int  col;
   int  mcnt, mbeg;
   double  sum;

   for (j=0; j<lp->nnbasic; j++){
      sum  = 0.0;
      col  = lp->nbaz[j];
      mcnt = lp->matcnt[col];
      mbeg = lp->matbeg[col];
      for (i=0; i<mcnt; i++)
         sum += lp->pIpiz[lp->matind[mbeg+i]] * lp->matval[mbeg+i];
      lp->pIdz[j] = -sum;
   }
}

void ILLfct_compute_yz (lpinfo *lp, svector *yz, svector *updz, int col)
{
   svector  a;

   a.nzcnt = lp->matcnt[col];
   a.indx  = &(lp->matind[lp->matbeg[col]]);
   a.coef  = &(lp->matval[lp->matbeg[col]]);

   ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, PIVZ_TOLER);
   if (updz)
      ILLbasis_column_solve_update (lp, &a, updz, yz);
   else
      ILLbasis_column_solve (lp, &a, yz);
   ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, SZERO_TOLER);
}

void ILLfct_compute_zz (lpinfo *lp, svector *zz, int row)
{
   double  e = 1.0;
   svector a;

   a.nzcnt = 1;
   a.coef = &e;
   a.indx = &row;

   ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, PIVZ_TOLER);
   ILLbasis_row_solve (lp, &a, zz);
   ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, SZERO_TOLER);
}

void ILLfct_compute_psteep_upv (lpinfo *lp, svector *swz)
{
   ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void ILLfct_compute_dsteep_upv (lpinfo *lp, svector *swz)
{
   ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int compute_zA1(lpinfo *lp, svector *z, svector *zA)
{
   int  rval = 0;
   int  i, j, nz = 0;
   int  col, mcnt, mbeg;
   double  sum;
   double  *v = (double *) NULL;

   ILL_SAFE_MALLOC (v, lp->nrows, double);
   for (i=0; i<lp->nrows; i++) v[i] = 0.0;
   for (i=0; i<z->nzcnt; i++)  v[z->indx[i]] = z->coef[i];

   for (j=0; j<lp->nnbasic; j++){
      sum  = 0.0;
      col  = lp->nbaz[j];
      mcnt = lp->matcnt[col];
      mbeg = lp->matbeg[col];
      for (i=0; i<mcnt; i++)
         sum += v[lp->matind[mbeg+i]] * lp->matval[mbeg+i];

      if (sum < -PIVZ_TOLER || sum > PIVZ_TOLER) {
         zA->coef[nz] = sum;
         zA->indx[nz] = j;
         nz ++;
      }
   }
   zA->nzcnt = nz;

 CLEANUP:
   ILL_IFFREE (v, double);
   ILL_RETURN (rval, "compute_zA1");
}

static int compute_zA2(lpinfo *lp, svector *z, svector *zA)
{
   int  rval = 0;
   int  i, j, nz = 0;
   int  row, col;
   int  rcnt, rbeg;
   double  val;
   double  *t = (double *) NULL;

   ILL_SAFE_MALLOC (t, lp->ncols, double);
   for (j=0; j<lp->ncols; j++) t[j] = 0.0;

   for (i=0; i<z->nzcnt; i++){
      row  = z->indx[i];
      val  = z->coef[i];
      rcnt = lp->rowcnt[row];
      rbeg = lp->rowbeg[row];
      for (j=0; j<rcnt; j++){
         col = lp->rowind[rbeg+j];
         if (lp->vstat[col] != STAT_BASIC)
            t[col] += val * lp->rowval[rbeg+j];
      }
   }

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];
      if (t[col] < -PIVZ_TOLER || t[col] > PIVZ_TOLER) {
         zA->coef[nz] = t[col];
         zA->indx[nz] = j;
         nz ++;
      }
   }
   zA->nzcnt = nz;

 CLEANUP:
   ILL_IFFREE (t, double);
   ILL_RETURN (rval, "compute_zA2");
}

static int compute_zA3(lpinfo *lp, svector *z, svector *zA)
{
   int  rval = 0;
   int  i, j, k, ix;
   int  nz = 0;
   int  row, col;
   int  rcnt, rbeg;
   double  val;

   k = 0;
   for (i=0; i<z->nzcnt; i++){
      row  = z->indx[i];
      val  = z->coef[i];
      rcnt = lp->rowcnt[row];
      rbeg = lp->rowbeg[row];
      for (j=0; j<rcnt; j++){
         col = lp->rowind[rbeg+j];
         if (lp->vstat[col] != STAT_BASIC){
            ix = lp->vindex[col];
            if (lp->iwork[ix] == 0){
               lp->iwork[ix] = 1;
               lp->work.indx[k ++] = ix;
            }
            lp->work.coef[ix] += val * lp->rowval[rbeg+j];
         }
      }
   }

   for (j=0; j<k; j++){
      ix  = lp->work.indx[j];
      val = lp->work.coef[ix];
      lp->work.coef[ix] = 0.0;
      lp->iwork[ix] = 0;
      if (val < -PIVZ_TOLER || val > PIVZ_TOLER) {
         zA->coef[nz] = val;
         zA->indx[nz] = ix;
         nz ++;
      }
   }
   zA->nzcnt = nz;

   ILL_RETURN (rval, "compute_zA3");
}

int ILLfct_compute_zA (lpinfo *lp, svector *z, svector *zA)
{
   if (z->nzcnt < lp->nrows / 2)
      return compute_zA3(lp, z, zA);
   else
      return compute_zA1(lp, z, zA);
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void ILLfct_update_basis_info (lpinfo *lp, int eindex, int lindex, int lvstat)
{
   int  evar;
   int  lvar;

   evar = lp->nbaz[eindex];

   if (lindex >= 0){ /* variable leaves basis */
      lvar = lp->baz[lindex];
      lp->vstat[evar]  = STAT_BASIC;
      lp->vstat[lvar]  = lvstat;
      lp->vindex[evar] = lindex;
      lp->vindex[lvar] = eindex;
      lp->baz[lindex]  = evar;
      lp->nbaz[eindex] = lvar;
      (lp->basisid) ++;
   }
   else{
      lp->vstat[evar]  = (lp->vstat[evar] == STAT_LOWER) ? STAT_UPPER : STAT_LOWER;
   }
}

void ILLfct_update_xz (lpinfo *lp, double tz, int eindex, int lindex)
{
   int  i, evar, estat;

   if (tz != 0.0)
      for (i=0; i<lp->yjz.nzcnt; i++)
         lp->xbz[lp->yjz.indx[i]] -= tz * lp->yjz.coef[i];

   if (lindex >= 0){ /* variable leaves basis */
      evar  = lp->nbaz[eindex];
      estat = lp->vstat[evar];

      if (estat == STAT_LOWER)
         lp->xbz[lindex] = lp->lz[evar] + tz;
      else if (estat == STAT_UPPER)
         lp->xbz[lindex] = lp->uz[evar] + tz;
      else if (estat == STAT_ZERO)
         lp->xbz[lindex] = tz;
   }
}

void ILLfct_update_piz (lpinfo *lp, double alpha)
{
   int  i;

   for (i=0; i<lp->zz.nzcnt; i++)
      lp->piz[lp->zz.indx[i]] += alpha * lp->zz.coef[i];
}

void ILLfct_update_pIpiz (lpinfo *lp, svector *z, double alpha)
{
   int  i;

   if (alpha == 0.0)
      return;
   else if (alpha == 1.0){
      for (i=0; i<z->nzcnt; i++)
         lp->pIpiz[z->indx[i]] += z->coef[i];
   }
   else{
      for (i=0; i<z->nzcnt; i++)
         lp->pIpiz[z->indx[i]] += alpha * z->coef[i];
   }
}

void ILLfct_update_dz (lpinfo *lp, int eindex, double alpha)
{
   int  i;

   for (i=0; i<lp->zA.nzcnt; i++)
      lp->dz[lp->zA.indx[i]] -= alpha * lp->zA.coef[i];
   lp->dz[eindex] = -alpha;
}

void ILLfct_update_pIdz (lpinfo *lp, svector *zA, int eindex, double alpha)
{
   int  i;

   if (alpha == 0.0)
      return;
   else if (alpha == 1.0){
      for (i=0; i<zA->nzcnt; i++)
         lp->pIdz[zA->indx[i]] -= zA->coef[i];
   }
   else{
      for (i=0; i<zA->nzcnt; i++)
         lp->pIdz[zA->indx[i]] -= alpha * zA->coef[i];
   }
   if (eindex > -1)
      lp->pIdz[eindex] = -alpha;
}

/* bound and coef shift routines */

/* scale bound in my_rand to get more random digits, unless bound is large */
static double my_rand (int bound, ILLrandstate *r)
{
   int  k = bound, scale = 1;
   double v = 0.0;

   if (bound < 100000) {k = 20000 * bound; scale = 20000;}
   v = 1 + (ILLutil_lprand (r) % (k));
   return  v / (double) scale;
}

static int expand_var_bounds (lpinfo *lp, double ftol, int *chgb)
{
   int  rval = 0;
   int  i, col, nchg = 0;
   double  x, l, u;
   double  newb, ran;
   double  cftol = fabs (ftol / 10.0);
   ILLrandstate r;

   ILLutil_sprand (1, &r);

   for (i=0; i<lp->nrows; i++){
      col = lp->baz[i];
      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE) continue;
      x   = lp->xbz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      if (l != NINFTY && x < l + ftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newb = MIN (x, l) - (cftol * ran);
         rval = ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
      if (u != INFTY && x > u - ftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newb = MAX (x, u) + (cftol * ran);
         rval = ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
   }
   *chgb = nchg;

 CLEANUP:
   ILL_RETURN (rval, "expand_var_bounds");
}

static int expand_phaseI_bounds (lpinfo *lp, int *chgb)
{
   int  rval = 0;
   int  i, col, nchg = 0;
   double  x, l, u;
   double  newb, ran;
   double  cftol = fabs (lp->tol->ip_tol / 10.0);
   ILLrandstate r;

   ILLutil_sprand (1, &r);

   for (i=0; i<lp->nrows; i++){
      col = lp->baz[i];
      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE) continue;
      x   = lp->xbz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      if (l != NINFTY && x < l + cftol && x > l - cftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newb = l - (cftol * ran);
         rval = ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
      if (u != INFTY && x > u - cftol && x < u + cftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newb = u + (cftol * ran);
         rval = ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
   }
   *chgb = nchg;

 CLEANUP:
   ILL_RETURN (rval, "expand_var_bounds");
}

int ILLfct_adjust_viol_bounds (lpinfo *lp)
{
   int  rval = 0;
   int  chgb = 0;

   rval = expand_var_bounds (lp, -lp->tol->pfeas_tol, &chgb);
#if FCT_DEBUG > 0
   if (rval == 0) printf ("adjusting %d bounds\n", chgb);
#endif
   ILL_RETURN (rval, "ILLfct_adjust_viol_bounds");
}

int ILLfct_perturb_bounds (lpinfo *lp)
{
   int  rval = 0;
   int  chgb = 0;

   rval = expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if FCT_DEBUG > 0
   if (rval == 0) printf ("perturbing %d bounds\n", chgb);
#endif
   ILL_RETURN (rval, "ILLfct_perturb_bounds");
}

int ILLfct_perturb_phaseI_bounds (lpinfo *lp)
{
   int  rval = 0;
   int  chgb = 0;

   rval = expand_phaseI_bounds (lp, &chgb);
#if FCT_DEBUG > 0
   if (rval == 0) printf ("perturbing %d phase I bounds\n", chgb);
#endif
   ILL_RETURN (rval, "ILLfct_perturb_phaseI_bounds");
}

int ILLfct_bound_shift (lpinfo *lp, int col, int bndtype, double newbnd)
{
   int  rval = 0;
   bndinfo *nbnd = (bndinfo *) NULL;

   ILL_SAFE_MALLOC (nbnd, 1, bndinfo);

   nbnd->varnum = col;
   nbnd->btype  = bndtype;
   if (bndtype == BOUND_LOWER){
      nbnd->pbound = lp->lz[col];
      nbnd->cbound = newbnd;
      lp->lz[col]  = newbnd;
  }
   else{
      nbnd->pbound = lp->uz[col];
      nbnd->cbound = newbnd;
      lp->uz[col]  = newbnd;
   }

   if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL){
      /* printf ("changing f/a bound\n"); */
      if (lp->lz[col] < lp->uz[col])
         lp->vtype[col] = VBOUNDED;
   }

   nbnd->next   = lp->bchanges;
   lp->bchanges = nbnd;
   lp->nbchange ++;

 CLEANUP:
   if (rval) ILL_IFFREE (nbnd, bndinfo);
   ILL_RETURN (rval, "ILLfct_bound_shift");
}

void ILLfct_unroll_bound_change (lpinfo *lp)
{
   int  col;
   int  changex = 0;
   bndinfo  *bptr = lp->bchanges;
   bndinfo  *nptr = (bndinfo *) NULL;

   while (lp->nbchange != 0){
      col = bptr->varnum;

      if (bptr->btype == BOUND_UPPER)
         lp->uz[col] = bptr->pbound;
      else
         lp->lz[col] = bptr->pbound;

      if (lp->vtype[col] == VBOUNDED){
         if (lp->lz[col] == lp->uz[col])
            lp->vtype[col] = (lp->lz[col] == 0.0) ? VARTIFICIAL : VFIXED;
      }

      if (lp->vstat[col] != STAT_BASIC){
         if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
             (bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
            changex ++;
      }
      nptr = bptr->next;
      ILL_IFFREE (bptr, bndinfo);
      bptr = nptr;
      lp->nbchange--;
   }
   lp->bchanges = bptr;
   if (changex)
      ILLfct_compute_xbz (lp);
}

static int expand_var_coefs (lpinfo *lp, double ftol, int *chgc)
{
   int  rval = 0;
   int  i, col, vs, vt;
   int  nchg = 0;
   double  dj, c;
   double  ran, newc;
   double  cftol = fabs (ftol / 10.0);
   ILLrandstate r;

   ILLutil_sprand (1, &r);

   for (i=0; i<lp->nnbasic; i++){
      dj  = lp->dz[i];
      col = lp->nbaz[i];
      c   = lp->cz[col];
      vs  = lp->vstat[col];
      vt  = lp->vtype[col];
      
      if (vt == VARTIFICIAL || vt == VFIXED)
         continue;

      if (vs == STAT_ZERO){
         newc = c - dj;
         rval = ILLfct_coef_shift (lp, col, newc);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
      if (vs == STAT_LOWER && dj < ftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newc = c - MIN (0.0, dj) + (cftol * ran);
         rval = ILLfct_coef_shift (lp, col, newc);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
      if (vs == STAT_UPPER && dj > -ftol){
         ran  = my_rand (50, &(lp->rstate)) + 1.0;
         newc = c - MAX (0.0, dj) - (cftol * ran);
         rval = ILLfct_coef_shift (lp, col, newc);
         ILL_CLEANUP_IF (rval);
         nchg ++;
      }
   }
   *chgc = nchg;

 CLEANUP:
   ILL_RETURN (rval, "expand_var_coefs");
}

int ILLfct_adjust_viol_coefs (lpinfo *lp)
{
   int  rval = 0;
   int  chgc = 0;

   rval = expand_var_coefs (lp, -lp->tol->dfeas_tol, &chgc);
#if FCT_DEBUG > 0
   if (rval == 0) printf ("perturbing %d coefs\n", chgc);
#endif
   ILL_RETURN (rval, "ILLfct_adjust_viol_coefs");
}

int ILLfct_perturb_coefs (lpinfo *lp)
{
   int  rval = 0;
   int  chgc = 0;

   rval = expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if FCT_DEBUG > 0
   if (rval == 0) printf ("perturbing %d coefs\n", chgc);
#endif
   ILL_RETURN (rval, "ILLfct_perturb_coefs");
}

int ILLfct_coef_shift (lpinfo *lp, int col, double newcoef)
{
   int  rval = 0;
   coefinfo  *ncoef = (coefinfo *) NULL;

   ILL_SAFE_MALLOC (ncoef, 1, coefinfo);

   ncoef->varnum = col;
   ncoef->pcoef  = lp->cz[col];
   ncoef->ccoef  = newcoef;
   lp->cz[col]   = newcoef;
   ncoef->next   = lp->cchanges;
   lp->cchanges  = ncoef;
   lp->dz[lp->vindex[col]] += (ncoef->ccoef - ncoef->pcoef);
   lp->ncchange ++;

 CLEANUP:
   if (rval) ILL_IFFREE (ncoef, coefinfo);
   ILL_RETURN (rval, "ILLfct_coef_shift");
}

void ILLfct_unroll_coef_change (lpinfo *lp)
{
   int  bascoef = 0;
   coefinfo  *cptr = (coefinfo *) lp->cchanges;
   coefinfo  *nptr = (coefinfo *) NULL;

   while (lp->ncchange != 0){
      lp->cz[cptr->varnum] = cptr->pcoef;
      if (lp->vstat[cptr->varnum] != STAT_BASIC)
         lp->dz[lp->vindex[cptr->varnum]] += (cptr->pcoef - cptr->ccoef);
      else
         bascoef ++;

      nptr = cptr->next;
      ILL_IFFREE (cptr, coefinfo);
      cptr = nptr;
      lp->ncchange--;
   }
   lp->cchanges = cptr;
   if (bascoef){
      ILLfct_compute_piz (lp);
      ILLfct_compute_dz (lp);
   }
}

/* feasibility routines */
void ILLfct_check_pfeasible (lpinfo *lp, feas_info *fs, double ftol)
{
   int  i, col;
   double  infeas = 0.0;

   fs->pstatus = PRIMAL_FEASIBLE;
   fs->totinfeas = 0.0;

   for (i=0; i<lp->nrows; i++){
      col = lp->baz[i];

      if (lp->uz[col] != INFTY && lp->xbz[i] > lp->uz[col] + ftol){
         infeas += lp->xbz[i] - lp->uz[col];
         lp->bfeas[i] = 1;
      }
      else if (lp->lz[col] != NINFTY && lp->xbz[i] < lp->lz[col] - ftol){
         infeas -= lp->xbz[i] - lp->lz[col];
         lp->bfeas[i] = -1;
      }
      else
         lp->bfeas[i] = 0;
   }
   if (infeas != 0.0){
      fs->pstatus = PRIMAL_INFEASIBLE;
      fs->totinfeas = infeas;
   }
   /*HHH*/ lp->pinfeas = infeas;
}

/* feasibility routines */
void ILLfct_check_pIpfeasible (lpinfo *lp, feas_info *fs, double ftol)
{
   int  i, col;
   int  ninf = 0;

   fs->pstatus = PRIMAL_FEASIBLE;
   fs->totinfeas = 0.0;

   for (i=0; i<lp->nrows; i++){
      col = lp->baz[i];

      if (lp->uz[col] != INFTY && lp->xbz[i] > ftol){
         ninf ++;
      }
      else if (lp->lz[col] != NINFTY && lp->xbz[i] < -ftol){
         ninf ++;
      }
   }
   if (ninf != 0)
      fs->pstatus = PRIMAL_INFEASIBLE;
}

void ILLfct_check_dfeasible (lpinfo *lp, feas_info *fs, double ftol)
{
   int  j, col;
   double  infeas = 0.0;

   fs->dstatus   = DUAL_FEASIBLE;
   fs->totinfeas = 0.0;

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED){
         lp->dfeas[j] = 0;
         continue;
      }

      if (lp->dz[j] < -ftol &&
          (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)){
         infeas -= lp->dz[j];
         lp->dfeas[j] = -1;
      }
      else if (lp->dz[j] > ftol &&
               (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)){
         infeas += lp->dz[j];
         lp->dfeas[j] = 1;
      }
      else
         lp->dfeas[j] = 0;
   }

   if (infeas != 0.0){
      fs->totinfeas = infeas;
      fs->dstatus = DUAL_INFEASIBLE;
   }
   /*HHH*/ lp->dinfeas = infeas;
}

void ILLfct_check_pIdfeasible (lpinfo *lp, feas_info *fs, double ftol)
{
   int  j, col;
   int  ninf = 0;
   double  *dz = lp->pIdz;

   fs->dstatus = DUAL_FEASIBLE;

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;

      if (dz[j] < -ftol &&
          (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)){
         ninf ++;
      }
      else if (dz[j] > ftol &&
               (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)){
         ninf ++;
      }
   }

   if (ninf != 0)
      fs->dstatus = DUAL_INFEASIBLE;
}

void ILLfct_dual_adjust (lpinfo *lp, double ftol)
{
   int  j, col;

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->dz[j] < -ftol && lp->uz[col] != INFTY){
         lp->vstat[col] = STAT_UPPER;
      }
      else if (lp->dz[j] > ftol && lp->lz[col] != -INFTY){
         lp->vstat[col] = STAT_LOWER;
      }
   }
}

void ILLfct_dphaseI_simple_update (lpinfo *lp, double ftol)
{
   int  j, col;

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->dz[j] < -ftol && lp->vtype[col] == VBOUNDED)
         lp->vstat[col] = STAT_UPPER;
      else if (lp->dz[j] > ftol && lp->vtype[col] == VBOUNDED)
         lp->vstat[col] = STAT_LOWER;
   }
}

/* set status values */
void ILLfct_set_status_values (lpinfo *lp, int pstatus, int dstatus, int ptype, int dtype)
{
   if (dstatus == DUAL_FEASIBLE && dtype == PHASEII){
      if (!lp->ncchange){
         lp->probstat.dual_feasible = 1;
         lp->basisstat.dual_feasible = 1;
         lp->basisstat.dual_infeasible = 0;
      }
   }
   if (dstatus == DUAL_INFEASIBLE && dtype == PHASEII){
      if (!lp->ncchange){
         lp->basisstat.dual_feasible = 0;
         lp->basisstat.dual_infeasible = 1;
      }
      if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEI)
         if (!lp->ncchange) lp->probstat.dual_infeasible = 1;
   }
   if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEII){
      if (!lp->nbchange){
         lp->probstat.primal_feasible = 1;
         lp->basisstat.primal_feasible = 1;
         lp->basisstat.primal_infeasible = 0;
      }
   }
   if (pstatus == PRIMAL_INFEASIBLE && ptype == PHASEII){
      lp->basisstat.primal_feasible = 0;
      lp->basisstat.primal_infeasible = 1;

      if (dstatus == DUAL_FEASIBLE && dtype == PHASEI)
         lp->probstat.primal_infeasible = 1;
   }
   if (pstatus == PRIMAL_UNBOUNDED){
      if (!lp->nbchange){
         lp->probstat.primal_unbounded = 1;
         lp->basisstat.primal_unbounded = 1;
         lp->probstat.dual_infeasible = 1;
         lp->basisstat.dual_infeasible = 1;
         lp->basisstat.dual_feasible = 0;
      }
   }
   if (dstatus == DUAL_UNBOUNDED){
      if (!lp->ncchange){
         lp->probstat.dual_unbounded = 1;
         lp->basisstat.dual_unbounded = 1;
         lp->probstat.primal_infeasible = 1;
         lp->basisstat.primal_infeasible = 1;
         lp->basisstat.primal_feasible = 0;
      }
   }
   if (lp->probstat.primal_feasible && lp->probstat.dual_feasible)
      lp->probstat.optimal = 1;

   if (lp->basisstat.primal_feasible && lp->basisstat.dual_feasible)
      lp->basisstat.optimal = 1;
   else
      lp->basisstat.optimal = 0;
}

void ILLfct_init_counts (lpinfo *lp)
{
   int  i;
   count_struct  *c = lp->cnts;

   c->ynz_cnt   = 0;
   c->num_y     = 0;
   c->y_ravg    = 1.0 + (double)lp->nrows / (PARAM_HEAP_RATIO * ILLutil_our_log2 (lp->nrows));
   c->znz_cnt   = 0;
   c->num_z     = 0;
   c->zanz_cnt  = 0;
   c->num_za    = 0;
   c->za_ravg   = 1.0 + (double)lp->nnbasic / (PARAM_HEAP_RATIO * ILLutil_our_log2 (lp->nnbasic));
   c->pnorm_cnt = 0;
   c->dnorm_cnt = 0;
   c->pinz_cnt  = 0;
   c->num_pi    = 0;
   c->pi1nz_cnt = 0;
   c->num_pi1   = 0;
   c->upnz_cnt  = 0;
   c->num_up    = 0;
   c->pupv_cnt  = 0;
   c->dupv_cnt  = 0;

   c->pI_iter  = 0;
   c->pII_iter = 0;
   c->dI_iter  = 0;
   c->dII_iter = 0;
   c->tot_iter = 0;

   for (i=0; i<10; i++){
      c->pivpI[i] = 0; c->pivpII[i] = 0;
      c->pivdI[i] = 0; c->pivdII[i] = 0;
   }
}

static void update_piv_values (count_struct *c, int phase, double piv)
{
   int  i = 0;
   double  v = 1.0;

   if (piv == 0.0) return;
   if (piv < 0.0) piv = -piv;
   while (piv < v && i<9) {v /= 10.0; i++;}
   if (phase == PRIMAL_PHASEI) c->pivpI[i] ++;
   if (phase == PRIMAL_PHASEII) c->pivpII[i] ++;
   if (phase == DUAL_PHASEI) c->pivdI[i] ++;
   if (phase == DUAL_PHASEII) c->pivdII[i] ++;
}

void ILLfct_update_counts (lpinfo *lp, int f, int upi, double upd)
{
   count_struct  *c = lp->cnts;

   switch (f){
   case CNT_PPHASE1ITER: c->pI_iter ++; c->tot_iter ++; break;
   case CNT_PPHASE2ITER: c->pII_iter ++; c->tot_iter ++; break;
   case CNT_DPHASE1ITER: c->dI_iter ++; c->tot_iter ++; break;
   case CNT_DPHASE2ITER: c->dII_iter ++; c->tot_iter ++; break;
   case CNT_YNZ:         c->ynz_cnt += upi; c->num_y ++; break;
   case CNT_ZANZ:        c->zanz_cnt += upi; c->num_za ++; break;
   case CNT_PINZ:        c->pinz_cnt += upi; c->num_pi ++; break;
   case CNT_P1PINZ:      c->pi1nz_cnt += upi; c->num_pi1 ++; break;
   case CNT_UPNZ:        c->upnz_cnt += upi;  c->num_up ++; break;
   case CNT_PIPIV:       update_piv_values (c, PRIMAL_PHASEI, upd); break;
   case CNT_PIIPIV:      update_piv_values (c, PRIMAL_PHASEII, upd); break;
   case CNT_DIPIV:       update_piv_values (c, DUAL_PHASEI, upd); break;
   case CNT_DIIPIV:      update_piv_values (c, DUAL_PHASEII, upd); break;
   case CNT_YRAVG:       c->y_ravg = (c->tot_iter*c->y_ravg + upi)/(1+c->tot_iter); break;
   case CNT_ZARAVG:      c->za_ravg = (c->tot_iter*c->za_ravg + upi)/(1+c->tot_iter); break;
   }
}

void ILLfct_print_counts (lpinfo *lp)
{
   int  i, niter;
   count_struct  *c = lp->cnts;

   c->tot_iter = c->pI_iter + c->pII_iter + c->dI_iter + c->dII_iter;
   niter = (c->tot_iter == 0) ? 1 : c->tot_iter;
   printf ("Counts for problem %s\n", lp->O->probname);
   if (c->num_y != 0) printf ("avg ynz = %.2f\n", (double) c->ynz_cnt / c->num_y);
   if (c->num_z != 0) printf ("avg znz = %.2f\n", (double) c->znz_cnt / c->num_z);
   if (c->num_za != 0) printf ("avg zanz = %.2f\n", (double) c->zanz_cnt / c->num_za);
   printf ("avg pnorm = %.2f\n", (double) c->pnorm_cnt / lp->nnbasic);
   printf ("avg dnorm = %.2f\n", (double) c->dnorm_cnt / lp->nrows);
   if (c->num_pi != 0) printf ("avg pinz = %.2f\n", (double) c->pinz_cnt / c->num_pi);
   if (c->num_pi1 != 0) printf ("avg piInz = %.2f\n", (double) c->pi1nz_cnt / c->num_pi1);
   if (c->num_up != 0) printf ("avg upnz = %.2f\n", (double) c->upnz_cnt / c->num_up);

   for (i=0; i<10; i++)
      printf ("piv 1.0e-%d : %d %d %d %d\n",
              i, c->pivpI[i], c->pivpII[i], c->pivdI[i], c->pivdII[i]);
}

static void copy_vectors (svector *a, svector *b)
{
   int i;
   int nzcnt = a->nzcnt;

   for (i=0; i<nzcnt; i++) {
      b->indx[i] = a->indx[i];
      b->coef[i] = a->coef[i];
   }
   b->nzcnt = nzcnt;
}

/* c <- a + t*b */
static void add_vectors (lpinfo *lp, svector *a, svector *b, svector *c, double t)
{
   int  i, r, l;
   svector  *w = &(lp->work);

   for (i=0; i<b->nzcnt; i++){
      r = b->indx[i];
      w->indx[i]   = r;
      w->coef[r]   = t * b->coef[i];
      lp->iwork[r] = 1;
   }
   l = b->nzcnt;

   for (i=0; i<a->nzcnt; i++){
      r = a->indx[i];
      if (lp->iwork[r] == 0) w->indx[l++] = r;
      w->coef[r] += a->coef[i];
   }
   for (i=0; i<l; i++){
      r = w->indx[i];
      c->indx[i]   = r;
      c->coef[i]   = w->coef[r];
      w->coef[r]   = 0.0;
      lp->iwork[r] = 0;
   }
   w->nzcnt = 0;
   c->nzcnt = l;
}

void ILLfct_update_pfeas (lpinfo *lp, int lindex, svector *srhs)
{
   int  i, k, r;
   int  col, nz=0;
   int  cbnd, f;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   int  tctr  = lp->upd.tctr;
   double  tz = lp->upd.tz;
   double  x, l, u;
   double  pftol  = lp->tol->ip_tol;
   double  *t = lp->upd.t;
   double  dty = 0.0;

   for (i=0; i < tctr && t[perm[i]] <= tz+fabs(tz)*1e-2; i++){
      cbnd = ix[perm[i]] % 10;
      if (cbnd == BBOUND) continue;
      k    = ix[perm[i]] / 10;
      r    = lp->yjz.indx[k];

      if (lp->iwork[r] != 1){
         lp->iwork[r] = 1;
         x   = lp->xbz[r];
         col = lp->baz[r];
         l   = lp->lz[col];
         u   = lp->uz[col];

         if (r != lindex){
            if (l != NINFTY && x < l-pftol) f = -1;
            else if (u != INFTY && x > u+pftol) f = 1;
            else f = 0;

            if (f != lp->bfeas[r]){
               srhs->indx[nz] = r;
               srhs->coef[nz] = f - lp->bfeas[r];
               dty += (f - lp->bfeas[r]) * lp->yjz.coef[k];
               nz++;
               lp->bfeas[r] = f;
            }
         }
         else{
            lp->bfeas[r] = 0;
         }
      }
   }
   while (--i >= 0){
      cbnd = ix[perm[i]] % 10;
      if (cbnd == BBOUND) continue;
      k    = ix[perm[i]] / 10;
      r    = lp->yjz.indx[k];
      lp->iwork[r] = 0;
   }
   srhs->nzcnt = nz;
   lp->upd.dty = dty;
}

void ILLfct_compute_ppIzz (lpinfo *lp, svector *srhs, svector *ssoln)
{
   if (srhs->nzcnt != 0){
      ILLbasis_row_solve (lp, srhs, ssoln);
   }
}

void ILLfct_update_ppI_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                         svector *ssoln, int eindex, int lindex, double alpha)
{
   if (lindex == -1){
      if (srhs->nzcnt != 0){
         ILLfct_update_pIpiz (lp, ssoln, 1.0);
         if (pinf->p_strategy == COMPLETE_PRICING){
            ILLfct_compute_zA (lp, ssoln, &(lp->zA));
            ILLfct_update_pIdz (lp, &(lp->zA), -1, 1.0);
         }
      }
      else{
         if (pinf->p_strategy == COMPLETE_PRICING)
            ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
         else
            ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
         return;
      }
   }
   else{
      if (srhs->nzcnt == 0){
         ILLfct_update_pIpiz (lp, &(lp->zz), alpha);
         if (pinf->p_strategy == COMPLETE_PRICING)
            ILLfct_update_pIdz (lp, &(lp->zA), eindex, alpha);
      }
      else{
         alpha -= lp->upd.dty/lp->upd.piv;
         add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), alpha);
         ILLfct_update_pIpiz (lp, &(lp->zz), 1.0);
         if (pinf->p_strategy == COMPLETE_PRICING){
            ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
            ILLfct_update_pIdz (lp, &(lp->zA), eindex, 1.0);
         }
      }
      lp->pIdz[eindex] = -alpha-lp->upd.fs;
   }
   if (pinf->p_strategy == COMPLETE_PRICING){
      ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
                                 PRIMAL_PHASEI);
      if (eindex > -1)
         ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
      ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, 0.0);
   }
   else
      ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
}

void ILLfct_update_dfeas (lpinfo *lp, int eindex, svector *srhs)
{
   int  i, j, k, c;
   int  cbnd, col, nz=0;
   int  vs, vt, f;
   int  delta;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   int  tctr  = lp->upd.tctr;
   int  mcnt, mbeg;
   double  dj;
   double  tz = lp->upd.tz;
   double  dftol  = lp->tol->id_tol;
   double  *t = lp->upd.t;
   double  *w = lp->work.coef;
   double  dty = 0.0;

   for (j=0; j < tctr && t[perm[j]] <= tz+tz*1e-2; j++){
      k = ix[perm[j]] / 10;
      c = lp->zA.indx[k];

      if (lp->iwork[c] != 1){
         lp->iwork[c] = 1;
         cbnd = ix[perm[j]] % 10;
         col  = lp->nbaz[c];
         dj   = lp->dz[c];
         vs   = lp->vstat[col];
         vt   = lp->vtype[col];

         if (cbnd == BSKIP){
            if (dj < -dftol && vs == STAT_LOWER)
               lp->vstat[col] = STAT_UPPER;
            else if (dj > dftol && vs == STAT_UPPER)
               lp->vstat[col] = STAT_LOWER;
         }
         else if (c != eindex){
            if (dj < -dftol && (vs == STAT_LOWER || vs == STAT_ZERO)) f = -1;
            else if (dj > dftol && (vs == STAT_UPPER || vs == STAT_ZERO))f = 1;
            else f = 0;

            if (f != lp->dfeas[c]){
               delta = f - lp->dfeas[c];
               mcnt  = lp->matcnt[col];
               mbeg  = lp->matbeg[col];
               for (i=0; i<mcnt; i++)
                  w[lp->matind[mbeg+i]] += delta*lp->matval[mbeg+i];

               dty += delta * lp->zA.coef[k];
               nz   = 1;
               lp->dfeas[c] = f;
            }
         }
         else{
            lp->dfeas[c] = 0;
         }
      }
   }
   while (--j >= 0){
      k    = ix[perm[j]] / 10;
      c    = lp->zA.indx[k];
      lp->iwork[c] = 0;
   }

   if (nz){
      for (i=0, nz=0; i<lp->nrows; i++)
         if (w[i] != 0.0) {
            srhs->coef[nz] = w[i];
            srhs->indx[nz] = i;
            nz++;
            w[i] = 0.0;
         }
   }

   srhs->nzcnt = nz;
   lp->upd.dty = dty;
}

void ILLfct_compute_dpIy (lpinfo *lp, svector *srhs, svector *ssoln)
{
   if (srhs->nzcnt != 0){
      ILLbasis_column_solve (lp, srhs, ssoln);
   }
}

void ILLfct_update_dpI_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                         svector *ssoln, int lindex, double alpha)
{
   int  i;

   if (srhs->nzcnt == 0){
      ILLfct_update_xz (lp, alpha, -1, -1);
   }
   else{
      alpha += lp->upd.dty/lp->upd.piv;
      add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), -alpha);
      for (i=0; i<lp->yjz.nzcnt; i++)
         lp->xbz[lp->yjz.indx[i]] += lp->yjz.coef[i];
   }
   lp->xbz[lindex] = alpha - lp->upd.fs;

   if (pinf->d_strategy == COMPLETE_PRICING){
      ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt, DUAL_PHASEI);
      ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
      ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, 0.0);
   }
   else
      ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
}

void ILLfct_update_dIIfeas (lpinfo *lp, int eindex, svector *srhs)
{
   int  j, k;
   int  col, indx, vs;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   int  tctr  = lp->upd.tctr;
   double  l, u, zAj, dj;
   double  theta = 0.0;
   double  delta = 0.0;
   double  dty   = 0.0;
   double  t_max = lp->upd.tz;
   double  *t = lp->upd.t;
   svector  a;

   srhs->nzcnt = 0;
   for (j=0; j<tctr && t[perm[j]] <= t_max; j++){
      k    = ix[perm[j]];
      indx = lp->zA.indx[k];

      if (indx != eindex){
         zAj = lp->zA.coef[k];
         col = lp->nbaz[indx];
         l   = lp->lz[col];
         u   = lp->uz[col];
         vs  = lp->vstat[col];
         dj  = lp->dz[j];

         delta = (vs == STAT_UPPER) ? (l-u) : (u-l);
         theta = delta * zAj;
         lp->vstat[col]
               = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;
         dty  += theta;

         a.nzcnt = lp->matcnt[col];
         a.indx  = &(lp->matind[lp->matbeg[col]]);
         a.coef  = &(lp->matval[lp->matbeg[col]]);
         add_vectors (lp, srhs, &a, srhs, delta);
      }
   }
   lp->upd.dty = dty;
}

void ILLfct_compute_dpIIy (lpinfo *lp, svector *srhs, svector *ssoln)
{
   if (srhs->nzcnt != 0){
      ILLbasis_column_solve (lp, srhs, ssoln);
   }
}

void ILLfct_update_dpII_prices (lpinfo *lp, price_info *pinf, svector *srhs,
                         svector *ssoln, int eindex, int lindex, double eval, double alpha)
{
   int  i;
   svector  *u;

   if (srhs->nzcnt == 0){
      ILLfct_update_xz (lp, alpha, -1, -1);
      u = &(lp->yjz);
   }
   else{
      if (ssoln->nzcnt != 0)
         for (i=0; i<ssoln->nzcnt; i++)
            lp->xbz[ssoln->indx[i]] -= ssoln->coef[i];
      ILLfct_update_xz (lp, alpha, -1, -1);
      add_vectors (lp, ssoln, &(lp->yjz), ssoln, 1.0);
      u = ssoln;
   }
   lp->xbz[lindex] = eval + alpha;

   if (pinf->d_strategy == COMPLETE_PRICING){
      ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
      ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
      ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, 0.0);
   }
   else
      ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int ILLfct_test_pivot (lpinfo *lp, int indx, int indxtype, double piv_val)
{
   int  i;
   double  pval = 0.0;

   if (indxtype == ROW_PIVOT){
      for (i=0; i<lp->yjz.nzcnt; i++)
         if (lp->yjz.indx[i] == indx){
            pval = lp->yjz.coef[i];
            break;
         }
   }
   else{
      for (i=0; i<lp->zA.nzcnt; i++)
         if (lp->zA.indx[i] == indx){
            pval = lp->zA.coef[i];
            break;
         }
   }

   if (fabs(pval - piv_val) / fabs (piv_val) > ALTPIV_TOLER){
#if FCT_DEBUG > 1
      if (indxtype == ROW_PIVOT)
         printf ("y_i = %.8f, z_j = %.8f\n", pval, piv_val);
      else
         printf ("z_j = %.8f, y_i = %.8f\n", pval, piv_val);
#endif
      return 1;
   }
   return 0;
}

#if FCT_DEBUG > 0

void fct_test_workvector (lpinfo *lp)
{
   int i, err=0;
   for (i=0; i<lp->ncols; i++){
      if (lp->work.coef[i] != 0.0) {err++; lp->work.coef[i] = 0.0;}
      if (lp->iwork[i] != 0) {err++; lp->iwork[i] = 0;}
   }
   if (err) printf ("bad work vector, err=%d\n", err);
}

void fct_test_pfeasible (lpinfo *lp)
{
   int  i, col;
   int  err = 0;
   double  ftol = lp->tol->pfeas_tol;

   for (i=0; i<lp->nrows; i++){
      col = lp->baz[i];

      if (lp->uz[col] != INFTY && lp->xbz[i] > lp->uz[col] + ftol){
         if (lp->bfeas[i] != 1) {err++; lp->bfeas[i] = 1;}
      }
      else if (lp->lz[col] != NINFTY && lp->xbz[i] < lp->lz[col] - ftol){
         if (lp->bfeas[i] != -1) {err++; lp->bfeas[i] = -1;}
      }
      /* else if (lp->bfeas[i] != 0) {err++; lp->bfeas[i] = 0;} */
   }
   if (err != 0) printf ("test_pfeas err =%d\n", err);
}

void fct_test_dfeasible (lpinfo *lp)
{
   int  j, col;
   int  err = 0;
   double  ftol = lp->tol->dfeas_tol;

   for (j=0; j<lp->nnbasic; j++){
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;

      if (lp->dz[j] < -ftol &&
          (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)){
         if (lp->dfeas[j] != -1) {err++; lp->dfeas[j] = -1;}
      }
      else if (lp->dz[j] > ftol &&
               (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)){
         if (lp->dfeas[j] != 1) {err++; lp->dfeas[j] = 1;}
      }
      /* else if (lp->dfeas[j] != 0) {err++; lp->dfeas[j] = 0;}*/
   }
   if (err != 0) printf ("test_dfeas err =%d\n", err);
}

void fct_test_pI_x (lpinfo *lp, price_info *p)
{
   int  i;
   double  *x = calloc (lp->nrows, sizeof (double));
   double  err = 0.0;
   int  ern = 0;

   for (i=0; i<lp->nrows; i++) x[i] = lp->xbz[i];
   ILLfct_compute_phaseI_xbz (lp);
   for (i=0; i<lp->nrows; i++) if (fabs (x[i]-lp->xbz[i]) > 1e-6)
      {err += fabs (x[i]-lp->xbz[i]); ern ++; printf ("bad i = %d\n",i);}
   if (err != 0.0) printf ("dI x err = %.7f, ern = %d\n", err, ern);
   ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
   free (x);
}

void fct_test_pII_x (lpinfo *lp, price_info *p)
{
   int  i;
   double  *x = calloc (lp->nrows, sizeof (double));
   double  err = 0.0;
   int  ern = 0;

   for (i=0; i<lp->nrows; i++) x[i] = lp->xbz[i];
   ILLfct_compute_xbz (lp);
   for (i=0; i<lp->nrows; i++) if (fabs (x[i]-lp->xbz[i]) > 1e-6)
      {err += fabs (x[i]-lp->xbz[i]); ern ++; printf ("bad i = %d\n",i);}
   if (err != 0.0) printf ("dII x err = %.7f, ern = %d\n", err, ern);
   ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
   free (x);
}

void fct_test_pI_pi_dz (lpinfo *lp, price_info *p)
{
   int  i;
   double  *pidz = calloc (lp->ncols, sizeof (double));
   double err = 0.0;
   int  ern = 0;

   for (i=0; i<lp->nrows; i++) pidz[i] = lp->pIpiz[i];
   ILLfct_compute_phaseI_piz (lp);
   for (i=0; i<lp->nrows; i++) if (fabs (pidz[i]-lp->pIpiz[i]) > 1e-6)
      {err += fabs (pidz[i]-lp->pIpiz[i]); ern ++;}
   if (err != 0.0) printf ("pI pi err = %.7f, ern = %d\n", err, ern);

   err = 0.0; ern = 0;
   for (i=0; i<lp->nnbasic; i++) pidz[i] = lp->pIdz[i];
   ILLfct_compute_phaseI_dz (lp);
   for (i=0; i<lp->nnbasic; i++) if (fabs (pidz[i]-lp->pIdz[i]) > 1e-6)
      {err += fabs (pidz[i]-lp->pIdz[i]); ern ++;}
   if (err != 0.0) printf ("pI dz err = %.7f, ern = %d\n", err, ern);
   ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
   free (pidz);
}

void fct_test_pII_pi_dz (lpinfo *lp, price_info *p)
{
   int  i;
   double  *pidz = calloc (lp->ncols, sizeof (double));
   double err = 0.0;
   int  ern = 0;

   for (i=0; i<lp->nrows; i++) pidz[i] = lp->piz[i];
   ILLfct_compute_piz (lp);
   for (i=0; i<lp->nrows; i++) if (fabs (pidz[i]-lp->piz[i]) > 1e-6)
      {err += fabs (pidz[i]-lp->piz[i]); ern ++;}
   if (err != 0.0) printf ("pII pi err = %.7f, ern = %d\n", err, ern);

   err = 0.0; ern = 0;
   for (i=0; i<lp->nnbasic; i++) pidz[i] = lp->dz[i];
   ILLfct_compute_dz (lp);
   for (i=0; i<lp->nnbasic; i++) if (fabs (pidz[i]-lp->dz[i]) > 1e-6)
      {err += fabs (pidz[i]-lp->dz[i]); ern ++;}
   if (err != 0.0) printf ("pII dz err = %.7f, ern = %d\n", err, ern);
   /*
   ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
   */
   free (pidz);
}

#endif
