/* RCS_INFO = "$RCSfile: simplex.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"

#include "stddefs.h"
#include "fct.h"
#include "ratio.h"
#include "price.h"
#include "basis.h"
#include "simplex.h"
#include "dstruct.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lib.h"  /* for ILLlib_writebasis */
#include "lp.h"   /* for ILLwrite_lp */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define DENSE_PI 0
#define DENSE_PIIPI 0
#define DENSE_NORM 0
#define SIMPLEX_DEBUG 0

static void
    init_lp_status_info (lp_status_info *ls),
    init_simplex_tols (lpinfo *lp),
    init_internal_lpinfo (lpinfo *lp),
    free_internal_lpinfo (lpinfo *lp),
    monitor_iter (lpinfo *lp, price_info *p, iter_info *it, int cphase),
    get_current_stat (lp_status_info *p, int algorithm, int *bstat);

static int
    build_internal_lpinfo (lpinfo *lp),
    terminate_simplex (lpinfo *lp, int phase, iter_info *it),
    primal_phaseI_step (lpinfo *lp, price_info *pinf, svector *updz,
                        svector *wz, iter_info *it),
    primal_phaseII_step (lpinfo *lp, price_info *pinf, svector *updz,
                         svector *wz, iter_info *it),
                        
    dual_phaseI_step (lpinfo *lp, price_info *pinf, svector *updz,
                      svector *wz, iter_info *it),
    dual_phaseII_step (lpinfo *lp, price_info *pinf, svector *updz,
                       svector *wz, iter_info *it), 
	report_value(lpinfo *lp, iter_info *it, 
                 const char* value_name, double value);


void ILLsimplex_init_lpinfo (lpinfo *lp)
{
    ILLbasis_init_basisinfo (lp);
    init_internal_lpinfo (lp);
}

void ILLsimplex_free_lpinfo (lpinfo *lp)
{
   if (lp) {
      ILL_IFFREE (lp->lz, double);
      ILL_IFFREE (lp->uz, double);
      ILL_IFFREE (lp->cz, double);
      ILLbasis_free_basisinfo (lp);
      free_internal_lpinfo (lp);
   }
}

void ILLsimplex_load_lpinfo (ILLlpdata *qslp, lpinfo *lp)
{
   lp->basisid  = -1;
   lp->maxiter  = 500000;
   lp->maxtime  = 20000;
   lp->iterskip = 100;
   lp->objbound = INFTY;
   lp->O        = qslp;
}

void ILLsimplex_set_bound (lpinfo *lp, double objbound, int sense)
{
   lp->objbound = (sense == ILL_MAX) ? -objbound : objbound;
}

static void init_lp_status_info (lp_status_info *ls)
{
   ls->optimal           = 0;
   ls->primal_feasible   = 0;
   ls->primal_infeasible = 0;
   ls->primal_unbounded  = 0;
   ls->dual_feasible     = 0;
   ls->dual_infeasible   = 0;
   ls->dual_unbounded    = 0;
}

static void init_simplex_tols (lpinfo *lp)
{
   lp->tol->pfeas_tol = PFEAS_TOLER;
   lp->tol->dfeas_tol = DFEAS_TOLER;
   lp->tol->pivot_tol = PIVOT_TOLER;
   lp->tol->szero_tol = SZERO_TOLER;

   lp->tol->ip_tol    = lp->tol->pfeas_tol / 2.0;
   lp->tol->id_tol    = lp->tol->dfeas_tol / 2.0;
}

static void init_internal_lpinfo (lpinfo *lp)
{
   lp->localrows = 0;
   lp->rowcnt = (int *) NULL;
   lp->rowbeg = (int *) NULL;
   lp->rowind = (int *) NULL;
   lp->rowval = (double *) NULL;

   lp->cz = (double *) NULL;
   lp->lz = (double *) NULL;
   lp->uz = (double *) NULL;

   lp->xbz   = (double *) NULL;
   lp->piz   = (double *) NULL;
   lp->dz    = (double *) NULL;
   lp->pIxbz = (double *) NULL;
   lp->pIpiz = (double *) NULL;
   lp->pIdz  = (double *) NULL;

   lp->vtype = (int *) NULL;

   ILLsvector_init (&(lp->zz));
   ILLsvector_init (&(lp->yjz));
   ILLsvector_init (&(lp->zA));
   ILLsvector_init (&(lp->work));
   ILLsvector_init (&(lp->srhs));
   ILLsvector_init (&(lp->ssoln));
   lp->iwork    = (int *) NULL;
   lp->upd.perm = (int *) NULL;
   lp->upd.ix   = (int *) NULL;
   lp->upd.t    = (double *) NULL;

   lp->bfeas = (int *) NULL;
   lp->dfeas = (int *) NULL;
   lp->tol   = (tol_struct *) NULL;
   lp->cnts  = (count_struct *) NULL;
   lp->bchanges = (bndinfo *) NULL;
   lp->cchanges = (coefinfo *) NULL;
}

static void free_internal_lpinfo (lpinfo *lp)
{
   bndinfo   *binfo = (bndinfo *) NULL;
   coefinfo  *cinfo = (coefinfo *) NULL;

   if (lp->localrows) { 
       ILL_IFFREE (lp->rowcnt, int);
       ILL_IFFREE (lp->rowbeg, int);
       ILL_IFFREE (lp->rowind, int);
       ILL_IFFREE (lp->rowval, double);
       lp->localrows = 0;
   }

   ILL_IFFREE (lp->cz, double);
   ILL_IFFREE (lp->lz, double);
   ILL_IFFREE (lp->uz, double);

   ILL_IFFREE (lp->xbz, double);
   ILL_IFFREE (lp->piz, double);
   ILL_IFFREE (lp->dz, double);
   ILL_IFFREE (lp->pIxbz, double);
   ILL_IFFREE (lp->pIpiz, double);
   ILL_IFFREE (lp->pIdz, double);

   ILL_IFFREE (lp->vtype, int);

   ILLsvector_free (&(lp->zz));
   ILLsvector_free (&(lp->yjz));
   ILLsvector_free (&(lp->zA));
   ILLsvector_free (&(lp->work));
   ILLsvector_free (&(lp->srhs));
   ILLsvector_free (&(lp->ssoln));
   ILL_IFFREE (lp->iwork, int);
   ILL_IFFREE (lp->upd.perm, int);
   ILL_IFFREE (lp->upd.ix, int);
   ILL_IFFREE (lp->upd.t, double);

   ILL_IFFREE (lp->bfeas, int);
   ILL_IFFREE (lp->dfeas, int);
   ILL_IFFREE (lp->tol, tol_struct);
   ILL_IFFREE (lp->cnts, count_struct);

   while (lp->bchanges){
      binfo = lp->bchanges;
      lp->bchanges = binfo->next;
      ILL_IFFREE (binfo, bndinfo);
   }

   while (lp->cchanges){
      cinfo = lp->cchanges;
      lp->cchanges = cinfo->next;
      ILL_IFFREE (cinfo, coefinfo);
   }
}

static int build_internal_lpinfo (lpinfo *lp)
{
   int  rval = 0;
   int  i, n; 
   ILLlpdata   *qslp = lp->O;
   ILLlp_sinfo *S = lp->O->sinfo;
   double      *lower, *upper, *obj;
   ILLlp_rows  lprows; 
   ILLmatrix   *A;

   init_lp_status_info (&(lp->probstat));
   init_lp_status_info (&(lp->basisstat));

   if (S != (ILLlp_sinfo *) NULL) {
       lp->nrows = S->nrows;
       lp->ncols = S->ncols;
       lp->bz    = S->rhs;
       lower = S->lower;
       upper = S->upper;
       obj   = S->obj;
       A     = &(S->A);
   } else {
       lp->nrows = qslp->nrows;
       lp->ncols = qslp->ncols;
       lp->bz    = qslp->rhs;
       lower = qslp->lower;
       upper = qslp->upper;
       obj   = qslp->obj;
       A     = &(qslp->A);
   }

   lp->matbeg = A->matbeg;
   lp->matcnt = A->matcnt;
   lp->matind = A->matind;
   lp->matval = A->matval;

   lp->nnbasic = lp->ncols - lp->nrows;

   ILL_SAFE_MALLOC (lp->lz, lp->ncols, double);
   ILL_SAFE_MALLOC (lp->uz, lp->ncols, double);
   ILL_SAFE_MALLOC (lp->cz, lp->ncols, double);
   if (!lp->lz || !lp->uz || !lp->cz) {
       fprintf (stderr, "build_internal_lpinfo\n");
       rval = 1; goto CLEANUP;
   }
   for (i = 0; i < lp->ncols; i++) {
       lp->lz[i] = lower[i];
       lp->uz[i] = upper[i];
       if (qslp->objsense ==  ILL_MAX) {
           lp->cz[i] = -obj[i];
       } else {
           lp->cz[i] = obj[i];
       }
   }

   if (!lp->O->rA) {
       rval = ILLlp_rows_init(&lprows, lp->O, 1); 
       ILL_CLEANUP_IF(rval);
       lp->rowbeg = lprows.rowbeg; 
       lp->rowcnt = lprows.rowcnt; 
       lp->rowind = lprows.rowind; 
       lp->rowval = lprows.rowval; 
       lp->localrows = 1;
   } else {
       /* row format exists, just use pointers */
       lp->rowbeg = lp->O->rA->rowbeg; 
       lp->rowcnt = lp->O->rA->rowcnt; 
       lp->rowind = lp->O->rA->rowind; 
       lp->rowval = lp->O->rA->rowval; 
       lp->localrows = 0;
   }

   ILL_SAFE_MALLOC (lp->xbz, lp->nrows, double);
   ILL_SAFE_MALLOC (lp->piz, lp->nrows, double);
   ILL_SAFE_MALLOC (lp->dz, lp->nnbasic, double);
   lp->final_phase = -1;
   lp->infub_ix    = -1;

   ILL_SAFE_MALLOC (lp->vtype, lp->ncols, int);

   rval = ILLsvector_alloc (&(lp->zz), lp->nrows);
   ILL_CLEANUP_IF (rval);
   rval = ILLsvector_alloc (&(lp->yjz), lp->nrows);
   ILL_CLEANUP_IF (rval);
   rval = ILLsvector_alloc (&(lp->zA), lp->nnbasic);
   ILL_CLEANUP_IF (rval);
   rval = ILLsvector_alloc (&(lp->work), lp->ncols);
   ILL_CLEANUP_IF (rval);
   rval = ILLsvector_alloc (&(lp->srhs), lp->nrows);
   ILL_CLEANUP_IF (rval);
   rval = ILLsvector_alloc (&(lp->ssoln), lp->nrows);
   ILL_CLEANUP_IF (rval);
   ILL_SAFE_MALLOC (lp->iwork, lp->ncols, int);
   for (i=0; i<lp->ncols; i++){
      lp->work.indx[i] = 0;
      lp->work.coef[i] = 0.0;
      lp->iwork[i] = 0;
   }
   n = 2*(lp->nrows)+1 > lp->ncols ? 2*(lp->nrows)+1 : lp->ncols;
   ILL_SAFE_MALLOC (lp->upd.perm, n, int);
   ILL_SAFE_MALLOC (lp->upd.ix,   n, int);
   ILL_SAFE_MALLOC (lp->upd.t,    n, double);

   ILL_SAFE_MALLOC (lp->bfeas, lp->nrows, int);
   ILL_SAFE_MALLOC (lp->dfeas, lp->nnbasic, int);

   ILL_SAFE_MALLOC (lp->tol, 1, tol_struct);
   init_simplex_tols (lp);
   ILL_SAFE_MALLOC (lp->cnts, 1, count_struct);
   ILLfct_init_counts (lp);

   lp->nbchange = 0;
   lp->ncchange = 0;

   lp->pIratio = RATIOTEST_HARRIS;
   lp->pIIratio = RATIOTEST_HARRIS;
   lp->dIratio = RATIOTEST_HARRIS;
   lp->dIIratio = RATIOTEST_HARRIS;
   lp->starttime = ILLutil_zeit ();
   ILLutil_sprand (1, &(lp->rstate));

 CLEANUP:
   if (rval) free_internal_lpinfo (lp);
   ILL_RETURN (rval, "build_internal_lpinfo");
}

int ILLsimplex_retest_psolution (lpinfo *lp, price_info *p, int phase, feas_info *fi)
{
   int rval    = 0;
   int fbid    = lp->fbasisid;
   int bid     = lp->basisid;
   double ptol = lp->tol->pfeas_tol;
   double dtol = lp->tol->dfeas_tol;
   double iptol= lp->tol->ip_tol;
   double idtol= lp->tol->id_tol;

   fi->pstatus = -1;
   fi->dstatus = -1;
   if (fbid < bid - PARAM_PRIMAL_REFACTORGAP) {
      rval = ILLbasis_refactor (lp);
      ILL_CLEANUP_IF (rval);
   }
   if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
      ILLfct_compute_xbz (lp);

   if (phase == PRIMAL_PHASEII){
      if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP) {
         ILLfct_compute_piz (lp);
         ILLfct_compute_dz (lp);
         if (p != NULL && p->p_strategy == COMPLETE_PRICING)
            ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
      }
      ILLfct_compute_pobj (lp);
      ILLfct_check_pfeasible (lp, fi, ptol);
      ILLfct_check_dfeasible (lp, fi, dtol);
   }
   else if (phase == PRIMAL_PHASEI){
      ILLfct_check_pfeasible (lp, fi, iptol);
      if (fi->pstatus != PRIMAL_FEASIBLE){
         ILLfct_compute_phaseI_piz (lp);
         ILLfct_compute_phaseI_dz (lp);
         ILLfct_check_pIdfeasible (lp, fi, idtol);
         if (p != NULL && p->p_strategy == COMPLETE_PRICING)
            ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
      }
   }
 CLEANUP:
   ILL_RETURN (rval, "ILLsimplex_retest_psolution");
}

int ILLsimplex_retest_dsolution (lpinfo *lp, price_info *p, int phase, feas_info *fi)
{
   int rval    = 0;
   int fbid    = lp->fbasisid;
   int bid     = lp->basisid;
   double ptol = lp->tol->pfeas_tol;
   double dtol = lp->tol->dfeas_tol;
   double iptol= lp->tol->ip_tol;
   double idtol= lp->tol->id_tol;

   fi->pstatus = -1;
   fi->dstatus = -1;
   if (fbid < bid - PARAM_DUAL_REFACTORGAP) {
      rval = ILLbasis_refactor (lp);
      ILL_CLEANUP_IF (rval);
   }
   if (fbid < bid - PARAM_DUAL_RESOLVEGAP) {
      ILLfct_compute_piz (lp);
      ILLfct_compute_dz (lp);
   }

   if (phase == DUAL_PHASEII){
      if (fbid < bid - PARAM_DUAL_RESOLVEGAP) {
         ILLfct_compute_xbz (lp);
         ILL_CLEANUP_IF (rval);
         if (p != NULL){
            if (p->d_strategy == COMPLETE_PRICING)
               ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
            else
               ILLprice_update_mpartial_price (lp, p, DUAL_PHASEII, ROW_PRICING);
         }
      }
      ILLfct_compute_dobj (lp);
      ILLfct_check_dfeasible (lp, fi, dtol);
      ILLfct_check_pfeasible (lp, fi, ptol);
   }
   else if (phase == DUAL_PHASEI){
      ILLfct_check_dfeasible (lp, fi, idtol);
      if (fi->dstatus != DUAL_FEASIBLE){
         ILLfct_compute_phaseI_xbz (lp);
         ILLfct_check_pIpfeasible (lp, fi, iptol);
         if (p != NULL){
            if (p->d_strategy == COMPLETE_PRICING)
               ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
            else
               ILLprice_update_mpartial_price (lp, p, DUAL_PHASEI, ROW_PRICING);
         }
      }
   }
 CLEANUP:
   ILL_RETURN (rval, "ILLsimplex_retest_dsolution");
}

int ILLsimplex_solution (lpinfo *lp, double *xz, double *piz, double *dz, double *objval)
{
   int  i, j;
   int  col;

   if (xz != NULL) {
      if (lp->basisstat.optimal == 0) { ILL_RETURN (1, "ILLsimplex_solution"); }
      for (i=0; i<lp->nrows; i++)
         xz[lp->baz[i]] = lp->xbz[i];
      for (j=0; j<lp->nnbasic; j++){
         col = lp->nbaz[j];
         if (lp->vstat[col] == STAT_UPPER)
            xz[col] = lp->uz[col];
         else if (lp->vstat[col] == STAT_LOWER)
            xz[col] = lp->lz[col];
         else
            xz[col] = 0;
      }
   }
   if (piz != NULL){
      if (lp->basisstat.optimal == 0) { ILL_RETURN (1, "ILLsimplex_solution"); }
      for (i=0; i<lp->nrows; i++)
         piz[i] = lp->piz[i];
   }
   if (dz != NULL){
      if (lp->basisstat.optimal == 0) { ILL_RETURN (1, "ILLsimplex_solution"); }
      for (i=0; i<lp->nrows; i++)
         dz[lp->baz[i]] = 0;
      for (j=0; j<lp->nnbasic; j++)
         dz[lp->nbaz[j]] = lp->dz[j];
   }
   if (objval != NULL)
      *objval = lp->objval;
   return 0;
}

int ILLsimplex_infcertificate (lpinfo *lp, double *pi)
{
   int  i, col, nz;
   char  *sense;
   double  x, l, u;
   lp_status_info  *ls;

   if (pi == NULL) return 0;

   ls = &(lp->basisstat);
   if (ls->primal_infeasible == 0 && ls->dual_unbounded == 0) {
      ILL_RETURN (1, "ILLsimplex_infcertificate");
   }

   if (lp->final_phase == PRIMAL_PHASEI && lp->pIpiz != NULL){
      for (i=0; i<lp->nrows; i++) pi[i] = lp->pIpiz[i];
   }
   else if (lp->final_phase == DUAL_PHASEII && lp->infub_ix != -1){
      col = lp->baz[lp->infub_ix];
      x   = lp->xbz[lp->infub_ix];
      l   = lp->lz[col];
      u   = lp->uz[col];

      for (i=0; i<lp->nrows; i++) pi[i] = 0.0;

      if (l != NINFTY && x < l){
	 for (i=0, nz=lp->zz.nzcnt; i<nz; i++)
	    pi[lp->zz.indx[i]] = -lp->zz.coef[i];
      }
      else{
	 for (i=0, nz=lp->zz.nzcnt; i<nz; i++)
	    pi[lp->zz.indx[i]] = lp->zz.coef[i];
      }
   }
   else{
     fprintf (stderr, "Invalid call to inf. certificate routine\n");
     ILL_RETURN (1, "ILLsimplex_infcertificate");
   }

   sense = lp->O->sense;
   for (i=0; i<lp->nrows; i++){
      if (sense[i] == 'G' && pi[i] < 0.0) pi[i] = 0.0;
      if (sense[i] == 'L' && pi[i] > 0.0) pi[i] = 0.0;
   }
   return 0;
}

#if SIMPLEX_DEBUG > 1
void test_cert (lpinfo *lp, double *pi)
{
   int  i, j;
   int  mcnt, mbeg;
   double  fsum = 0.0, sum;

   for (i=0; i<lp->nrows; i++){
      if (lp->O->sense[i] == 'G' && pi[i] < 0.0) printf ("compl \n");
      if (lp->O->sense[i] == 'L' && pi[i] > 0.0) printf ("compll \n");
   }

   for (i=0; i<lp->nrows; i++)
      fsum += pi[i] * lp->bz[i];

   for (j=0; j<lp->nnbasic; j++){
      sum  = 0.0;
      mcnt = lp->matcnt[j];
      mbeg = lp->matbeg[j];
      for (i=0; i<mcnt; i++)
         sum += pi[lp->matind[mbeg+i]] * lp->matval[mbeg+i];

      if (sum < -1e-6 && (lp->vtype[j] == VUPPER || lp->vtype[j] == VFREE))
         printf ("compl1\n");
      if (sum > 1e-6 && (lp->vtype[j] == VLOWER || lp->vtype[j] == VFREE))
         printf ("compl2\n");

      if (sum < 0.0 && lp->lz[j] != NINFTY)
         fsum += (-sum)*lp->lz[j];
      else if (sum > 0.0 && lp->uz[j] != INFTY)
         fsum -= sum*lp->uz[j];
   }
   printf ("fsum = %.8f\n", fsum);
}
#endif

static void save_paraminfo (price_info *pinf, iter_info *it)
{
   param_info  *pr = &(it->oldinfo);

   pr->origalgo = it->algorithm;
   pr->pphaseI  = pinf->pI_price;
   pr->pphaseII = pinf->pII_price;
   pr->dphaseI  = pinf->dI_price;
   pr->dphaseII = pinf->dII_price;
   pr->p_strategy = pinf->p_strategy;
   pr->d_strategy = pinf->d_strategy;
}

static void restore_paraminfo (iter_info *it, price_info *pinf)
{
   param_info  *pr = &(it->oldinfo);

   it->algorithm    = pr->origalgo;
   pinf->pI_price   = pr->pphaseI;
   pinf->pII_price  = pr->pphaseII;
   pinf->dI_price   = pr->dphaseI;
   pinf->dII_price  = pr->dphaseII;
   pinf->p_strategy = pr->p_strategy;
   pinf->d_strategy = pr->d_strategy;
}

int ILLsimplex (lpinfo *lp, int algorithm, ILLlp_basis *B, price_info *pinf,
                int *status, int sdisplay)
{
   int  phase     = -1;
   int  singular  = -1;
   int  rval      = 0;
   int  new_price = -1;
   svector  wz;
   svector  updz;
   feas_info  fi;
   iter_info  it;

   it.newphase  = -1;
   it.nextphase = -1;
   it.nextstep  = -1;
   it.sdisplay  = sdisplay;
   it.itercnt   = 0;
   it.solstatus = ILL_LP_UNSOLVED;
   it.curtime   = 0;
   it.rounds    = 0;
   it.prevobj   = INFTY;
   it.nosolve   = 0;
   it.noprog    = 0;
   it.objtol    = OBJBND_TOLER;
   it.chkobj    = PARAM_MAX_NOPROG;
   it.inner     = 0;
   it.algorithm = algorithm;
   it.pricetype = -1;
   it.resumeid  = -1;
   save_paraminfo (pinf, &it);

#if SIMPLEX_DEBUG > 0
   if (lp->O->nrows > 1000) it.sdisplay = 1;
#endif
   if (status) *status = QS_LP_UNSOLVED;

   free_internal_lpinfo (lp);
   init_internal_lpinfo (lp);
   rval = build_internal_lpinfo (lp);
   ILL_CLEANUP_IF (rval);

   ILLsvector_init (&wz);
   rval = ILLsvector_alloc (&wz, lp->nrows);
   ILL_CLEANUP_IF (rval);
   ILLsvector_init (&updz);
   rval = ILLsvector_alloc (&updz, lp->nrows);
   ILL_CLEANUP_IF (rval);

   if (it.sdisplay){
	  char buffer[256];
      sprintf (buffer, "starting ILLsimplex on %s...\n", lp->O->probname);
      /* depending on LP's reporter 
	       string is printed to stdout 
	       or handed to GUI */
	  rval = rval || ILLstring_report(buffer, &lp->O->reporter); 
      printf ("Problem has %d rows and %d cols\n", lp->nrows, lp->ncols);
	  fflush (stdout);
   }
   ILLfct_set_variable_type (lp);

   if (B != (ILLlp_basis *) NULL) {
      rval = ILLbasis_load (lp, B);
      ILL_CLEANUP_IF (rval);
      if (it.algorithm == DUAL_SIMPLEX) {
         if (B->rownorms) {
            rval = ILLprice_load_rownorms (lp, B->rownorms, pinf);
            ILL_CLEANUP_IF (rval);
         }
         else
            ILL_IFFREE (pinf->dsinfo.norms, double);
      }
      else if (it.algorithm == PRIMAL_SIMPLEX) {
         if (B->colnorms) {
            rval = ILLprice_load_colnorms (lp, B->colnorms, pinf);
            ILL_CLEANUP_IF (rval);
         }
         else 
            ILL_IFFREE (pinf->psinfo.norms, double);
      }
      else if (it.algorithm != PRIMAL_OR_DUAL){
         fprintf (stderr, "Unknown algorithm %d in ILLsimplex\n", it.algorithm);
         rval = 1;  ILL_CLEANUP;
      }
   }
   else if (lp->basisid == -1){
      if (lp->nrows < 200 && lp->ncols < 400)
         rval = ILLbasis_get_initial (lp, it.algorithm);
      else
         rval = ILLbasis_get_cinitial (lp, it.algorithm);
      ILL_CLEANUP_IF (rval);
      ILLprice_free_pricing_info (pinf);
   }

   if (lp->fbasisid != lp->basisid){
      rval = ILLbasis_factor (lp, &singular);
      ILL_CLEANUP_IF (rval);
      if (singular) ILLprice_free_pricing_info (pinf);
   }

 START:
   it.solstatus = ILL_LP_UNSOLVED;
   init_lp_status_info (&(lp->basisstat));

   ILLfct_compute_piz (lp);
   ILLfct_compute_dz (lp);
   if (it.algorithm == DUAL_SIMPLEX){
      if (B != NULL || it.resumeid == SIMPLEX_RESUME_UNSHIFT)
         ILLfct_dual_adjust (lp, lp->tol->dfeas_tol);
      else
         ILLfct_dual_adjust (lp, 0.0);
   }
   ILLfct_compute_xbz (lp);

   ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
   ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
   ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

   if (fi.dstatus == DUAL_FEASIBLE && lp->objbound != INFTY){
      ILLfct_compute_dobj (lp);
      if (lp->dobjval > lp->objbound){
         it.solstatus = ILL_BND_REACHED;
         goto TERMINATE;
      }
   }
   if (fi.pstatus == PRIMAL_FEASIBLE && fi.dstatus == DUAL_FEASIBLE){
      it.solstatus = ILL_LP_SOLVED;
      ILLfct_compute_pobj (lp);
      goto TERMINATE;
   }

   if (it.algorithm == PRIMAL_OR_DUAL){
      if (fi.pstatus == PRIMAL_FEASIBLE) it.algorithm = PRIMAL_SIMPLEX;
      else if (fi.dstatus == DUAL_FEASIBLE) it.algorithm = DUAL_SIMPLEX;
      else if (lp->pinfeas < 10*lp->dinfeas) it.algorithm = PRIMAL_SIMPLEX;
      else it.algorithm = DUAL_SIMPLEX;
   }

   if (it.algorithm == PRIMAL_SIMPLEX){
      if (fi.pstatus == PRIMAL_FEASIBLE) phase = PRIMAL_PHASEII;
      else                               phase = PRIMAL_PHASEI;
   }
   else if (it.algorithm == DUAL_SIMPLEX){
      if (fi.dstatus == DUAL_FEASIBLE) phase = DUAL_PHASEII;
      else                             phase = DUAL_PHASEI;
   }

   rval = ILLprice_build_pricing_info (lp, pinf, phase);
   ILL_CLEANUP_IF (rval);

   it.newphase = SIMPLEX_PHASE_NEW;
   it.nextstep = SIMPLEX_CONTINUE;

   while (it.nextstep == SIMPLEX_CONTINUE){

      if (phase == PRIMAL_PHASEI){
         rval = primal_phaseI_step (lp, pinf, &updz, &wz, &it);
         ILL_CLEANUP_IF (rval);
      }

      else if (phase == PRIMAL_PHASEII){
         rval = primal_phaseII_step (lp, pinf, &updz, &wz, &it);
         ILL_CLEANUP_IF (rval);
      }

      else if (phase == DUAL_PHASEI){
         rval = dual_phaseI_step (lp, pinf, &updz, &wz, &it);
         ILL_CLEANUP_IF (rval);
      }

      else if (phase == DUAL_PHASEII){
         rval = dual_phaseII_step (lp, pinf, &updz, &wz, &it);
         ILL_CLEANUP_IF (rval);
      }
      if (it.nextstep == SIMPLEX_RESUME){
         ILLprice_free_pricing_info (pinf);

         if (it.resumeid == SIMPLEX_RESUME_UNSHIFT){
            if (it.pricetype == QS_PRICE_PDEVEX){
               pinf->pI_price   = QS_PRICE_PDEVEX;
               pinf->pII_price  = QS_PRICE_PDEVEX;
            }
            else if (it.pricetype == QS_PRICE_DDEVEX){
               pinf->dI_price   = QS_PRICE_DDEVEX;
               pinf->dII_price  = QS_PRICE_DDEVEX;
            }
         }
         else if (it.resumeid == SIMPLEX_RESUME_NUMER){
            ILLfct_unroll_bound_change (lp);
            ILLfct_unroll_coef_change (lp);
            rval = ILLbasis_get_initial (lp, it.algorithm);
            ILL_CLEANUP_IF (rval);
            rval = ILLbasis_factor (lp, &singular);
            ILL_CLEANUP_IF (rval);
         }
         it.pricetype = -1;
         goto START;
      }
      else if (it.nextstep == SIMPLEX_CONTINUE){
         it.itercnt++;

         if (it.nextphase != phase){
            it.newphase = SIMPLEX_PHASE_NEW;
            phase       = it.nextphase;
            new_price   = ILLprice_get_price (pinf, phase);

            if (pinf->cur_price != new_price){
               ILLprice_free_pricing_info (pinf);
               rval = ILLprice_build_pricing_info (lp, pinf, phase);
               ILL_CLEANUP_IF (rval);
            }
         }
      }
   }

#if SIMPLEX_DEBUG > 0
   ILLfct_print_counts (lp);
#endif
   rval = terminate_simplex (lp, phase, &it);
   ILL_CLEANUP_IF (rval);

 TERMINATE:
   restore_paraminfo (&it, pinf);

   if (it.sdisplay) {
       printf ("completed ILLsimplex\n");
       printf ("%s: ", lp->O->probname);
       fflush (stdout); 
   }

   if (status) {
      if (it.solstatus == ILL_MAX_ITER) {
         *status = QS_LP_ITER_LIMIT;
      } else if (it.solstatus == ILL_MAX_TIME) {
         *status = QS_LP_TIME_LIMIT;
      } else if (it.solstatus == ILL_LP_ABORTED) {
         *status = QS_LP_ABORTED;
      } else if (it.solstatus == ILL_PPHASEI_ERROR ||
                 it.solstatus == ILL_PPHASEII_ERROR ||
                 it.solstatus == ILL_DPHASEI_ERROR ||
                 it.solstatus == ILL_DPHASEII_ERROR ||
                 it.solstatus == ILL_LP_UNSOLVED) {
         *status = QS_LP_UNSOLVED;
      } else if (it.solstatus == ILL_LP_SOLVED) {
         if (lp->basisstat.optimal) {
            *status = QS_LP_OPTIMAL;
         } else if (lp->basisstat.primal_infeasible ||
                    lp->basisstat.dual_unbounded) {
            *status = QS_LP_INFEASIBLE;
         } else if (lp->basisstat.primal_unbounded) {
            *status = QS_LP_UNBOUNDED;
         }
      } else {
         fprintf (stderr, "unknown solution status in ILLsimplex\n");
         rval = 1;
         ILL_CLEANUP_IF (rval);
      }
   }

#if SIMPLEX_DEBUG > 1
   {
      int rva = 0;
      double *pi = NULL;
      ILL_SAFE_MALLOC (pi, lp->nrows, double);

      rva = ILLsimplex_infcertificate (lp, pi);
      printf ("rva = %d\n", rva);
      if (!rva){
         test_cert (lp, pi);
      }
      ILL_IFFREE (pi, double);
   }
#endif
   if (it.sdisplay){
      int bstat = 0;

      printf ("time = %.3f, pI = %d, pII = %d, dI = %d, dII = %d, ",
             ILLutil_zeit () - lp->starttime, lp->cnts->pI_iter, lp->cnts->pII_iter,
             lp->cnts->dI_iter, lp->cnts->dII_iter);
      fflush (stdout);
      get_current_stat (&(lp->basisstat), it.algorithm, &bstat);
      switch (bstat){
      case OPTIMAL:
         printf ("opt = %f\n", lp->objval);
         break;
      case PRIMAL_INFEASIBLE:
         printf ("no primal soln\n");
         break;
      case PRIMAL_UNBOUNDED:
         printf ("primal unbounded\n");
         break;
      case PRIMAL_FEASIBLE:
         printf ("primal obj = %f\n", lp->pobjval);
         break;
      case DUAL_INFEASIBLE:
         printf ("no dual soln\n");
         break;
      case DUAL_UNBOUNDED:
         printf ("dual unbounded\n");
         break;
      case DUAL_FEASIBLE:
         printf ("dual obj = %f\n", lp->dobjval);
         break;
      }
      fflush (stdout);

      if (it.sdisplay > 1){
         if (it.algorithm == PRIMAL_SIMPLEX && pinf->pI_price == QS_PRICE_PDEVEX)
            printf ("Devex norms initialised %d times\n", pinf->pdinfo.ninit);
         fflush (stdout);
      }
   }

 CLEANUP:
   ILLsvector_free (&wz);
   ILLsvector_free (&updz);
   ILL_RETURN (rval, "ILLsimplex");
}

static int terminate_simplex (lpinfo *lp, int phase, iter_info *it)
{
   int  rval  = 0;
   int  sphase;
   feas_info  fi;

   if (it->solstatus != ILL_MAX_TIME && it->solstatus != ILL_MAX_ITER)
      ILL_CLEANUP;

   if (it->algorithm == PRIMAL_SIMPLEX){
      if (lp->nbchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d bound shifts\n", lp->nbchange); fflush (stdout);
         }
         ILLfct_unroll_bound_change (lp);
      }
      rval = ILLsimplex_retest_psolution (lp, NULL, phase, &fi);
      ILL_CLEANUP_IF (rval);

      sphase = (phase == PRIMAL_PHASEI) ? PHASEI : PHASEII;
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, sphase);
   }
   else if (it->algorithm == DUAL_SIMPLEX){
      if (lp->ncchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d coef shifts\n", lp->ncchange); fflush (stdout);
         }
         ILLfct_unroll_coef_change (lp);
      }
      rval = ILLsimplex_retest_dsolution (lp, NULL, phase, &fi);
      ILL_CLEANUP_IF (rval);

      sphase = (phase == DUAL_PHASEI) ? PHASEI : PHASEII;
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, sphase, PHASEII);
   }

 CLEANUP:
   ILL_RETURN (rval, "terminate_simplex");
}

static int test_progress (double objval, double prevobj)
{
   double  denom = 0.0;

   denom = (fabs (objval) < 1e-7) ? 1.0 : objval;
   if (fabs ((objval - prevobj) / denom) < 1e-5)
      return 0;
   else
      return 1;
}

static void monitor_iter (lpinfo *lp, price_info *p, iter_info *it, int phase)
{
   double  print_val = 0.0;
   double  tottime  = ILLutil_zeit () - lp->starttime;
   int     curtime  = (int) ILLutil_our_floor (tottime); /* MONIKA */
   char    print_str[20];
   feas_info  fi;
   int aborted = 0; 

   /* one of the following two time display mechanisms */
   switch (phase){
   case PRIMAL_PHASEI:
      print_val = lp->pinfeas; strcpy (print_str,"primal infeas"); break;
   case PRIMAL_PHASEII:
      print_val = lp->pobjval; strcpy (print_str,"primal objval"); break;
   case DUAL_PHASEI:
      print_val = lp->dinfeas; strcpy (print_str,"dual infeas"); break;
   case DUAL_PHASEII:
      print_val = lp->dobjval; strcpy (print_str, "dual objval"); break;
   }

   aborted = report_value(lp, it, print_str, print_val); 
   /*if (it->sdisplay && it->itercnt % lp->iterskip == 0) {
	  // printf ("(%d): %s = %f\n", it->itercnt, print_str, print_val);
	  // fflush (stdout);
	  }*/
   if (curtime != it->curtime){
      it->curtime = curtime;
      /*
      if (it->sdisplay){
         printf ("time = %d.0, ", curtime);
         printf ("(%d): %s = %f\n", it->itercnt, print_str, print_val);
         fflush (stdout);
      }
      */
   }

   if (phase == DUAL_PHASEII && lp->objbound != INFTY){
      if (lp->dobjval > lp->objbound + it->objtol){
         ILLfct_unroll_coef_change (lp);
         ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
         ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

         if (fi.dstatus == DUAL_FEASIBLE){
            ILLfct_compute_dobj (lp);
            if (lp->dobjval > lp->objbound){
               it->solstatus = ILL_BND_REACHED;
               it->nextstep  = SIMPLEX_TERMINATE;
               if (it->sdisplay) {printf ("bound reached\n"); fflush (stdout);}
            }
            else
               it->objtol *= 10.0;
         }
         else{
            it->nextphase = DUAL_PHASEI;
            it->newphase  = SIMPLEX_PHASE_NEW;
            it->objtol   *= 5.0;
         }
      }
   }
   if (it->itercnt >= lp->maxiter){
      it->solstatus = ILL_MAX_ITER;
      it->nextstep = SIMPLEX_TERMINATE;
      if (it->sdisplay) {printf ("iter limit reached\n"); fflush (stdout);}
      ILL_CLEANUP;
   }
   else if (tottime >= lp->maxtime){
      it->solstatus = ILL_MAX_TIME;
      it->nextstep = SIMPLEX_TERMINATE;
      if (it->sdisplay) {printf ("time limit reached\n"); fflush (stdout);}
      ILL_CLEANUP;
   } else if (aborted) {
	  it->solstatus = ILL_LP_ABORTED;
      it->nextstep = SIMPLEX_TERMINATE;
      if (it->sdisplay) {printf ("aborted\n"); fflush (stdout);}
      ILL_CLEANUP;
   } 
   /*
   if (it->rounds && it->inner){
      it->inner --;
      if (it->inner == 0){
         printf ("restoring ..\n");
         restore_paraminfo (it, p);
         it->newphase   = SIMPLEX_PHASE_NEW;
         it->nextstep   = SIMPLEX_RESUME;
         it->resumeid   = SIMPLEX_RESUME_OUTER;
         ILL_CLEANUP;
      }
   }
   */
   if (phase == DUAL_PHASEII){
      if (it->noprog > it->chkobj){
         ILLfct_perturb_coefs (lp);
         it->noprog  = 0;
         it->prevobj = lp->dobjval;
      }
   }
   else if (phase == PRIMAL_PHASEII){
      if (it->noprog > it->chkobj){
         ILLfct_perturb_bounds (lp);
         it->noprog  = 0;
         it->prevobj = lp->pobjval;
      }
   }
   else if (phase == PRIMAL_PHASEI){
      if (it->noprog > it->chkobj){
         it->algorithm  = DUAL_SIMPLEX;
         it->nextstep   = SIMPLEX_RESUME;
         it->resumeid   = SIMPLEX_RESUME_NUMER;
      }
   }
   else if (phase == DUAL_PHASEI){
      if (it->noprog > it->chkobj){
         it->algorithm  = PRIMAL_SIMPLEX;
         it->nextstep   = SIMPLEX_RESUME;
         it->resumeid   = SIMPLEX_RESUME_NUMER;
      }
   }
 CLEANUP:
   return;
}

static int primal_phaseI_step (lpinfo *lp, price_info *pinf, svector *updz,
                               svector *wz, iter_info *it)
{
   int  rval     = 0;
   int  singular = 0;
   int  refactor = 0;
   int  cphase   = PRIMAL_PHASEI;
   double  alpha = 0.0;
   feas_info  fi;
   ratio_res  rs;
   price_res  pr;

   ILLfct_update_counts (lp, CNT_PPHASE1ITER, 0, 0.0);
   it->nextstep    = SIMPLEX_CONTINUE;
   it->nextphase   = PRIMAL_PHASEI;
   lp->final_phase = PRIMAL_PHASEI;
   it->nosolve ++;
    
   if (it->newphase != 0){
      ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
      if (it->newphase == SIMPLEX_PHASE_NEW){
         it->noprog  = 0;
         if (it->sdisplay){
            printf ("starting primal phase I\n"); fflush (stdout);
         }
      }
      it->newphase = 0;
      it->nosolve  = 0;
      it->prevobj  = lp->pinfeas;
      ILL_SAFE_MALLOC (lp->pIpiz, lp->nrows, double);
      ILL_SAFE_MALLOC (lp->pIdz, lp->nnbasic, double);

      ILLfct_compute_phaseI_piz (lp);
      if (pinf->p_strategy == COMPLETE_PRICING){
         ILLfct_compute_phaseI_dz (lp);
#if USEHEAP > 0
         ILLprice_free_heap (pinf);
#endif
         ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEI);
#if USEHEAP > 0
         rval = ILLprice_test_for_heap (lp, pinf, lp->nnbasic,pinf->d_scaleinf,
                                        PRIMAL_SIMPLEX, 0);
         ILL_CLEANUP_IF (rval);
#endif
      }
      else if (pinf->p_strategy == MULTI_PART_PRICING)
         ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
   }

   monitor_iter (lp, pinf, it, cphase);
   if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
       it->newphase != 0)
      ILL_CLEANUP;

   ILLprice_primal (lp, pinf, &pr, cphase);

   if (pr.price_stat == PRICE_OPTIMAL){
      if (it->sdisplay > 1) {
         printf ("primal phase I seemingly done\n");
         printf ("retesting soln\n");
         fflush (stdout);
      }
      rval = ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);

      ILL_CLEANUP_IF (rval);
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEI);

      if (fi.pstatus == PRIMAL_FEASIBLE){
         it->nextphase = PRIMAL_PHASEII;
      }
      else if (fi.dstatus == DUAL_FEASIBLE){
         it->solstatus = ILL_LP_SOLVED;
         it->nextstep  = SIMPLEX_TERMINATE;
      }
      ILL_CLEANUP;
   }

   ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
   ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, 0.0);
   ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, 0.0);

   ILLratio_pI_test (lp, pr.eindex, pr.dir, &rs);

   if (rs.ratio_stat == RATIO_FAILED){
      /*
      rval = E_SIMPLEX_ERROR;
      it->solstatus = ILL_PPHASEI_ERROR;
      */
      it->algorithm = DUAL_SIMPLEX;
      it->nextstep  = SIMPLEX_RESUME;
      it->resumeid  = SIMPLEX_RESUME_NUMER;
      ILL_CLEANUP;
   }
   else if (rs.ratio_stat == RATIO_NEGATIVE){
      double itol = lp->tol->ip_tol;

      lp->tol->ip_tol = 0.0;
      lp->pinfeas += lp->upd.c_obj;
      if (!test_progress (lp->pinfeas, it->prevobj))
         it->noprog ++;
      else{
         it->prevobj = lp->pinfeas;
         it->noprog  = 0;
      }
      ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
      lp->tol->ip_tol = itol;
      ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
      ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex, rs.lindex, 0.0);
   }
   else if (rs.ratio_stat == RATIO_NOBCHANGE){
      lp->pinfeas += lp->upd.c_obj;
      if (!test_progress (lp->pinfeas, it->prevobj))
         it->noprog ++;
      else{
         it->prevobj = lp->pinfeas;
         it->noprog  = 0;
      }

      ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
      ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
      ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
      ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
      fct_test_workvector (lp);
      fct_test_pfeasible (lp);
#endif
      ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex, rs.lindex, 0.0);
   }
   else if (rs.ratio_stat == RATIO_BCHANGE){
      alpha = lp->pIdz[pr.eindex] / rs.pivotval;
      lp->pinfeas += lp->upd.c_obj;

      if (!test_progress (lp->pinfeas, it->prevobj)){
         if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
             lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL){
            if (it->noprog > 0) it->noprog --;
         }
         else
            it->noprog ++;
      }
      else{
         it->prevobj = lp->pinfeas;
         it->noprog  = 0;
      }

      ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
      ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, 0.0);
      if (pinf->p_strategy == COMPLETE_PRICING){
         ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
         ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, 0.0);

         if (pinf->pI_price == QS_PRICE_PSTEEP)
            ILLfct_compute_psteep_upv (lp, wz);
      }

      rval = ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex, rs.pivotval);
      ILL_CLEANUP_IF (rval);

      ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
      ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
      ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
      ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
      fct_test_workvector (lp);
      fct_test_pfeasible (lp);
#endif
      rval = ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
      ILL_CLEANUP_IF (rval);

      if (singular){
         it->nextstep = SIMPLEX_RESUME;
         it->resumeid = SIMPLEX_RESUME_SING;
         ILL_CLEANUP;
      }
      if (!refactor){
         ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex, rs.lindex, alpha);
      }
      if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE){
         ILLfct_compute_xbz (lp);
         ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
         ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);
         if (fi.pstatus == PRIMAL_FEASIBLE) it->nextphase = PRIMAL_PHASEII;

         it->newphase = SIMPLEX_PHASE_RECOMP;
         ILL_CLEANUP;
      }
   }

#if DENSE_PI > 1
   fct_test_workvector (lp);
   fct_test_pi_dz (lp, pinf);
#endif

 CLEANUP:
   if (it->nextphase != PRIMAL_PHASEI || it->nextstep == SIMPLEX_RESUME ||
       it->newphase  != 0 || rval != 0){
      ILL_IFFREE (lp->pIpiz, double);
      ILL_IFFREE (lp->pIdz, double);
   }
   ILL_RETURN (rval, "primal_phaseI_step");
}

static int primal_phaseII_step (lpinfo *lp, price_info *pinf, svector *updz,
                                svector *wz, iter_info *it)
{
   int  boundch;
   int  rval     = 0;
   int  bndtype  = 0;
   int  singular = 0;
   int  refactor = 0;
   int  ratio_iter = 0;
   int  cphase   = PRIMAL_PHASEII;
   double  lbound;
   double  alpha;
   feas_info  fi;
   ratio_res  rs;
   price_res  pr;

   ILLfct_update_counts (lp, CNT_PPHASE2ITER, 0, 0.0);
   it->nextstep    = SIMPLEX_CONTINUE;
   it->nextphase   = PRIMAL_PHASEII;
   lp->final_phase = PRIMAL_PHASEII;
   it->nosolve ++;

   if (it->newphase != 0){
      ILLfct_compute_pobj (lp);
      if (it->newphase == SIMPLEX_PHASE_NEW){
         it->noprog  = 0;
         if (it->sdisplay) {
            printf ("starting primal phase II\n"); fflush (stdout);
         }
      }
      it->newphase = 0;
      it->nosolve  = 0;
      it->prevobj  = lp->pobjval;
      ILLfct_compute_piz (lp);
      if (pinf->p_strategy == COMPLETE_PRICING){
         ILLfct_compute_dz (lp);
#if USEHEAP > 0
         ILLprice_free_heap (pinf);
#endif
         ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEII);
#if USEHEAP > 0
         rval = ILLprice_test_for_heap (lp, pinf, lp->nnbasic,pinf->d_scaleinf,
                                        PRIMAL_SIMPLEX, 0);
         ILL_CLEANUP_IF (rval);
#endif
      }
      else if (pinf->p_strategy == MULTI_PART_PRICING){
         ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
      }
   }

   monitor_iter (lp, pinf, it, cphase);
   if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
       it->newphase != 0)
      ILL_CLEANUP;

   ILLprice_primal (lp, pinf, &pr, cphase);

   if (pr.price_stat == PRICE_OPTIMAL){
      if (lp->nbchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d bound shifts\n", lp->nbchange);
            fflush (stdout);
         }
         ILLfct_unroll_bound_change (lp);
         ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
         ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);

         /*HHH*/ ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol); 
         /*HHH* printf ("primal (opt) infeas %.6f\n", lp->pinfeas); fflush (stdout); 
         *HHH* printf ("dual (opt) infeas %.6f\n", lp->dinfeas); fflush (stdout);*/

         if (fi.pstatus != PRIMAL_FEASIBLE){
            it->algorithm = DUAL_SIMPLEX;
            it->nextstep  = SIMPLEX_RESUME;
            it->resumeid  = SIMPLEX_RESUME_UNSHIFT;
            it->pricetype = QS_PRICE_DDEVEX;
            ILL_CLEANUP;
            /*
            it->nextphase = PRIMAL_PHASEI;
            lp->tol->ip_tol /= 5.0;
            lp->tol->id_tol /= 5.0;
            ILL_CLEANUP;
            */
         }
      }

      if (it->sdisplay > 1) {
         printf ("problem seemingly solved\n");
         printf ("seemingly opt = %f\nretesting soln\n",lp->pobjval);
         fflush (stdout);
      }
      rval = ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);
      ILL_CLEANUP_IF (rval);
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

      if (fi.pstatus == PRIMAL_INFEASIBLE){
         it->nextphase = PRIMAL_PHASEI;
         lp->tol->ip_tol /= 5.0;
         lp->tol->id_tol /= 5.0;
      }
      else if (fi.dstatus == DUAL_FEASIBLE){
         it->solstatus = ILL_LP_SOLVED;
         lp->objval   = lp->pobjval;
         it->nextstep = SIMPLEX_TERMINATE;
      }
      ILL_CLEANUP;
   }

   ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
   ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, 0.0);
   ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, 0.0);
   ratio_iter = 0;
   do{
      ILLratio_pII_test (lp, pr.eindex, pr.dir, &rs);
      lbound  = rs.lbound;
      boundch = rs.boundch;
      ratio_iter ++;

      if (boundch){
         /*
         if (ratio_iter > PARAM_PRATIOTESTS){
            lbound = lp->xbz[rs.lindex];
            boundch = 0;
         }
         */
         boundch = 0;
         bndtype = (rs.lvstat == STAT_UPPER) ? BOUND_UPPER : BOUND_LOWER;
         rval = ILLfct_bound_shift (lp, lp->baz[rs.lindex], bndtype, lbound);
         ILL_CLEANUP_IF (rval);
      }
   } while (boundch);

   if (rs.ratio_stat == RATIO_FAILED){
      /*
      rval = E_SIMPLEX_ERROR;
      it->solstatus = ILL_PPHASEII_ERROR;
      */
      it->algorithm = DUAL_SIMPLEX;
      it->nextstep  = SIMPLEX_RESUME;
      it->resumeid  = SIMPLEX_RESUME_NUMER;
      ILL_CLEANUP;
   }
   else if (rs.ratio_stat == RATIO_UNBOUNDED){
      if (lp->nbchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d bound shifts\n", lp->nbchange);
            fflush (stdout);
         }
         ILLfct_unroll_bound_change (lp);
      }
      ILLfct_set_status_values (lp, PRIMAL_UNBOUNDED, -1, PHASEII, -1);
      it->solstatus = ILL_LP_SOLVED;
      it->nextstep = SIMPLEX_TERMINATE;
      ILL_CLEANUP;
   }
   else if (rs.ratio_stat == RATIO_NOBCHANGE){
      lp->pobjval += rs.tz * lp->dz[pr.eindex];
      lp->objval = lp->pobjval;
      if (!test_progress (lp->pobjval, it->prevobj))
         it->noprog ++;
      else{
         it->prevobj = lp->pobjval;
         it->noprog  = 0;
      }

      ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
      ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
      if (pinf->p_strategy == COMPLETE_PRICING)
         ILLprice_compute_dual_inf (lp, pinf, &pr.eindex, 1, PRIMAL_PHASEII);
      else if (pinf->p_strategy == MULTI_PART_PRICING)
         ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
   }
   else if (rs.ratio_stat == RATIO_BCHANGE){
      alpha = lp->dz[pr.eindex] / rs.pivotval;
      lp->pobjval += rs.tz * lp->dz[pr.eindex];
      lp->objval = lp->pobjval;

      if (!test_progress (lp->pobjval, it->prevobj)){
         if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
             lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL){
            if (it->noprog > 0) it->noprog --;
         }
         else
            it->noprog ++;
      }
      else{
         it->prevobj = lp->pobjval;
         it->noprog  = 0;
      }

      ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
      ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, 0.0);
      if (pinf->p_strategy == COMPLETE_PRICING){
         ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
         ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, 0.0);
         if (pinf->pII_price == QS_PRICE_PSTEEP)
            ILLfct_compute_psteep_upv (lp, wz);
      }
      rval = ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex, rs.pivotval);
      ILL_CLEANUP_IF (rval);

      ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
      ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
      rval = ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
      ILL_CLEANUP_IF (rval);

      if (singular){
         it->nextstep = SIMPLEX_RESUME;
         it->resumeid = SIMPLEX_RESUME_SING;
         ILL_CLEANUP;
      }
      if (!refactor){
         ILLfct_update_piz (lp, alpha);

         if (pinf->p_strategy == COMPLETE_PRICING){
            ILLfct_update_dz (lp, pr.eindex, alpha);
            ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt, PRIMAL_PHASEII);
            ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, 0.0);
         }
         else if (pinf->p_strategy == MULTI_PART_PRICING){
            ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
         }
      }
      if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE){ 
         ILLfct_compute_xbz (lp);
         it->newphase = SIMPLEX_PHASE_RECOMP;
      }
   }

 CLEANUP:
   ILL_RETURN (rval, "primal_phaseII_step");
}

static int dual_phaseI_step (lpinfo *lp, price_info *pinf, svector *updz,
                             svector *wz, iter_info *it)
{
   int  rval = 0;
   int  singular = 0;
   int  refactor = 0;
   int  cphase = DUAL_PHASEI;
   double  alpha;
   double  alpha1 = 0.0;
   feas_info  fi;
   ratio_res  rs;
   price_res  pr;

   ILLfct_update_counts (lp, CNT_DPHASE1ITER, 0, 0.0);
   it->nextstep    = SIMPLEX_CONTINUE;
   it->nextphase   = DUAL_PHASEI;
   lp->final_phase = DUAL_PHASEI;
   it->nosolve ++;

   if (it->newphase != 0){
      ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
      if (it->newphase == SIMPLEX_PHASE_NEW){
         it->noprog  = 0;
         if (it->sdisplay) {
            printf ("starting dual phase I\n"); fflush (stdout);
         }
      }
      it->newphase = 0;
      it->nosolve  = 0;
      it->prevobj  = lp->dinfeas;

      ILLfct_compute_phaseI_xbz (lp);
      if (pinf->d_strategy == COMPLETE_PRICING){
#if USEHEAP > 0
         ILLprice_free_heap (pinf);
#endif
         ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEI);
#if USEHEAP > 0
         rval = ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
                                        DUAL_SIMPLEX, 0);
         ILL_CLEANUP_IF (rval);
#endif
      }
      else if (pinf->d_strategy == MULTI_PART_PRICING){
         ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
      }
   }

   monitor_iter (lp, pinf, it, cphase);
   if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
       it->newphase != 0)
      ILL_CLEANUP;

   ILLprice_dual (lp, pinf, cphase, &pr);

   if (pr.price_stat == PRICE_OPTIMAL){
      if (it->sdisplay > 1) {
         printf ("dual phase I seemingly done\n");
         printf ("retesting soln\n"); fflush (stdout);
      }

      rval = ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
      ILL_CLEANUP_IF (rval);
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEI, PHASEII);

      if (fi.dstatus == DUAL_FEASIBLE){
         it->nextphase = DUAL_PHASEII;
      }
      else if (fi.pstatus == PRIMAL_FEASIBLE){
         it->solstatus = ILL_LP_SOLVED;
         it->nextstep  = SIMPLEX_TERMINATE;
      }
      it->newphase = SIMPLEX_PHASE_NEW;
      ILL_CLEANUP;
   }

   ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
   ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
   ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, 0.0);
   ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, 0.0);

   ILLratio_dI_test (lp, pr.lindex, pr.lvstat, &rs);

   if (rs.ratio_stat == RATIO_FAILED){
      /*
      rval = E_SIMPLEX_ERROR;
      it->solstatus = ILL_DPHASEI_ERROR;
      */
      it->algorithm = PRIMAL_SIMPLEX;
      it->nextstep  = SIMPLEX_RESUME;
      it->resumeid  = SIMPLEX_RESUME_NUMER;
      ILL_CLEANUP;
   }
   else if (rs.ratio_stat == RATIO_BCHANGE){
      ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[rs.eindex]);
      rval = ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
      if (rval){
         rval = ILLbasis_factor (lp, &singular);
         ILL_CLEANUP_IF (rval);
         if (singular == 0) refactor = 1;
         goto END;
      }
      ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, 0.0);
      ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, 0.0);

      if (pinf->dI_price == QS_PRICE_DSTEEP)
         ILLfct_compute_dsteep_upv (lp, wz);
      rval = ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex, rs.pivotval);
      ILL_CLEANUP_IF (rval);

      lp->dinfeas -= lp->upd.c_obj;

      if (!test_progress (lp->dinfeas, it->prevobj)){
         if (lp->vtype[lp->baz[pr.lindex]] == VARTIFICIAL ||
             lp->vtype[lp->nbaz[rs.eindex]] == VFREE){
            if (it->noprog > 0) it->noprog --;
         }
         else
            it->noprog ++;
      }
      else{
         it->prevobj = lp->dinfeas;
         it->noprog  = 0;
      }

      alpha  = lp->dz[rs.eindex] / rs.pivotval;
      alpha1 = lp->xbz[pr.lindex] / rs.pivotval;

      ILLfct_update_piz (lp, alpha);
      ILLfct_update_dz (lp, rs.eindex, alpha);
      ILLfct_update_dfeas (lp, rs.eindex, &(lp->srhs));
      ILLfct_compute_dpIy (lp, &(lp->srhs), &(lp->ssoln));
      ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);

#if DENSE_PI > 0
      fct_test_workvector (lp);
      fct_test_dfeasible (lp);
#endif
      rval = ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
      ILL_CLEANUP_IF (rval);

#if DENSE_NORM > 0
      test_dsteep_norms (lp, pinf);
#endif

      ILLfct_update_dpI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.lindex, alpha1);

   END:
      if (singular){
         it->nextstep = SIMPLEX_RESUME;
         it->resumeid = SIMPLEX_RESUME_SING;
         ILL_CLEANUP;
      }
      if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE){
         ILLfct_compute_piz (lp);
         ILLfct_compute_dz (lp);
         ILLfct_dual_adjust (lp, 0.0);
         ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
         ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);
         if (fi.dstatus == DUAL_FEASIBLE) it->nextphase = DUAL_PHASEII;

         it->newphase = SIMPLEX_PHASE_RECOMP;
         ILL_CLEANUP;
      }
   }

#if DENSE_PI > 1
   fct_test_workvector (lp);
   fct_test_pI_x (lp, pinf);
#endif

 CLEANUP:
   ILL_RETURN (rval, "dual_phaseI_step");
}

static int dual_phaseII_step (lpinfo *lp, price_info *pinf, svector *updz,
                              svector *wz, iter_info *it)
{
   int  coeffch;
   int  rval = 0;
   int  singular = 0;
   int  refactor = 0;
   int  ratio_iter = 0;
   int  cphase = DUAL_PHASEII;
   int  lcol, ecol;
   int  estat, newphase;
   double  x_bi, v_l, eval;
   double  alpha;
   double  alpha1 = 0.0;
   double  ecoeff;
   feas_info  fi;
   ratio_res  rs;
   price_res  pr;

   ILLfct_update_counts (lp, CNT_DPHASE2ITER, 0, 0.0);
   it->nextstep    = SIMPLEX_CONTINUE;
   it->nextphase   = DUAL_PHASEII;
   lp->final_phase = DUAL_PHASEII;
   newphase        = it->newphase;
   it->nosolve ++;

   if (it->newphase != 0){
      ILLfct_compute_dobj (lp);
      if (it->newphase == SIMPLEX_PHASE_NEW){
         it->noprog  = 0;
         if (it->sdisplay) {
            printf ("starting dual phase II\n"); fflush (stdout);
         }
      }
      it->newphase = 0;
      it->nosolve  = 0;
      it->prevobj  = lp->dobjval;
      ILLfct_compute_xbz (lp);

      if (pinf->d_strategy == COMPLETE_PRICING){
#if USEHEAP > 0
         ILLprice_free_heap (pinf);
#endif
         ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEII);
#if USEHEAP > 0
         rval = ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
                                        DUAL_SIMPLEX, 0);
         ILL_CLEANUP_IF (rval);
#endif
      }
      else if (pinf->d_strategy == MULTI_PART_PRICING){
         ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
      }
   }

   monitor_iter (lp, pinf, it, cphase);
   if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
       it->newphase != 0)
      ILL_CLEANUP;

   ILLprice_dual (lp, pinf, cphase, &pr);

   if (pr.price_stat == PRICE_OPTIMAL){
      if (lp->ncchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d coef shifts\n", lp->ncchange); fflush (stdout);
         }
         ILLfct_unroll_coef_change (lp);
         ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
         ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

         /*HHH*/ ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
         /*HHH* printf ("dual (opt) infeas %.6f\n", lp->dinfeas); fflush (stdout);
         *HHH* printf ("primal (opt) infeas %.6f\n", lp->pinfeas); fflush (stdout);*/

         if (fi.dstatus != DUAL_FEASIBLE){
            it->algorithm = PRIMAL_SIMPLEX;
            it->nextstep  = SIMPLEX_RESUME;
            it->resumeid  = SIMPLEX_RESUME_UNSHIFT;
            it->pricetype = QS_PRICE_PDEVEX;
            ILL_CLEANUP;
            /*
            it->nextphase = DUAL_PHASEI;
            lp->tol->ip_tol /= 5.0;
            lp->tol->id_tol /= 5.0;
            ILL_CLEANUP;
            */
         }
      }
      if (it->sdisplay > 1) {
         printf ("problem seemingly solved\n");
         printf ("seemingly dual opt = %f\n",lp->dobjval);
         printf ("retesting soln\n");
         fflush (stdout);
      }

      rval = ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
      ILL_CLEANUP_IF (rval);
      ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

      if (fi.dstatus == DUAL_INFEASIBLE){
         it->nextphase = DUAL_PHASEI;
         lp->tol->ip_tol /= 5.0;
         lp->tol->id_tol /= 5.0;
      }
      else if (fi.pstatus == PRIMAL_FEASIBLE){
         lp->objval = lp->dobjval;
         it->solstatus = ILL_LP_SOLVED;
         it->nextstep = SIMPLEX_TERMINATE;
      }
      ILL_CLEANUP;
   }

   ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
   ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
   ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, 0.0);
   ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, 0.0);

   ratio_iter = 0;
   do{
      ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
      if (rs.ratio_stat == RATIO_NEGATIVE){
         if (it->sdisplay > 1){
            printf ("adjust coefs to remove negative ratio tests\n"); fflush (stdout);
         }
         ILLfct_adjust_viol_coefs (lp);
         ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
         if (rs.ratio_stat == RATIO_NEGATIVE){
            printf ("internal error: bad ratio test\n"); fflush (stdout);
            rs.ratio_stat = RATIO_FAILED;
            break;
         }
      }

      coeffch = rs.coeffch;
      ecoeff  = rs.ecoeff;
      ratio_iter ++;

      if (coeffch){
         /*
         if (ratio_iter > PARAM_DRATIOTESTS){
            ecoeff = lp->cz[lp->nbaz[rs.eindex]] - lp->dz[rs.eindex];
            coeffch = 0;
         }
         */
         coeffch = 0;
         rval = ILLfct_coef_shift (lp, lp->nbaz[rs.eindex], ecoeff);
         ILL_CLEANUP_IF (rval);
      }
      if (rs.ratio_stat == RATIO_BCHANGE)
         if (lp->vstat[lp->nbaz[rs.eindex]] == STAT_ZERO)
            break;

   } while (coeffch);

   if (rs.ratio_stat == RATIO_FAILED){
      /*
      rval = E_SIMPLEX_ERROR;
      it->solstatus = ILL_DPHASEII_ERROR;
      */
      it->algorithm = PRIMAL_SIMPLEX;
      it->nextstep  = SIMPLEX_RESUME;
      it->resumeid  = SIMPLEX_RESUME_NUMER;
      ILL_CLEANUP;
   }
   else if (rs.ratio_stat == RATIO_UNBOUNDED){
      lp->infub_ix = pr.lindex;
      if (lp->ncchange != 0){
         if (it->sdisplay > 1) {
            printf ("unrolling %d coef shifts\n", lp->ncchange); fflush (stdout);
         }
         ILLfct_unroll_coef_change (lp);
      }
      ILLfct_set_status_values (lp, -1, DUAL_UNBOUNDED, -1, PHASEII);
      it->solstatus = ILL_LP_SOLVED;
      it->nextstep = SIMPLEX_TERMINATE;
   }
   else if (rs.ratio_stat == RATIO_BCHANGE){
      lcol = lp->baz[pr.lindex];
      ecol = lp->nbaz[rs.eindex];

      ILLfct_compute_yz (lp, &(lp->yjz), updz, ecol);
      ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, 0.0);
      ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, 0.0);
      rval = ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
      if (rval != 0){
         if (newphase == 0){
            rval = ILLbasis_factor (lp, &singular);
            ILL_CLEANUP_IF (rval);
            if (singular == 0) refactor = 1;
            goto END;
         }
         else{
            if (it->sdisplay > 1){
               printf ("warning: bad step\n"); fflush (stdout);
            }
         }
      }

      lp->dobjval += lp->upd.c_obj;
      lp->objval   = lp->dobjval;

      if (!test_progress (lp->dobjval, it->prevobj)){
         if (lp->vtype[lcol] == VARTIFICIAL || lp->vtype[ecol] == VFREE){
            if (it->noprog > 0) it->noprog --;
         }
         else
            it->noprog ++;
      }
      else{
         it->prevobj = lp->dobjval;
         it->noprog  = 0;
      }

      if (pinf->dII_price == QS_PRICE_DSTEEP)
         ILLfct_compute_dsteep_upv (lp, wz);
      rval = ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex, rs.pivotval);
      ILL_CLEANUP_IF (rval);

      x_bi  = lp->xbz[pr.lindex];
      v_l   = (pr.lvstat == STAT_LOWER) ? lp->lz[lcol] : lp->uz[lcol];
      alpha = (pr.lvstat == STAT_LOWER) ? -rs.tz : rs.tz;
      estat = lp->vstat[ecol];
      if (estat == STAT_LOWER)     eval = lp->lz[ecol];
      else if (estat == STAT_ZERO) eval = 0.0;
      else                         eval = lp->uz[ecol];

      ILLfct_update_piz (lp, alpha);
      ILLfct_update_dz (lp, rs.eindex, alpha);
      ILLfct_update_dIIfeas (lp, rs.eindex, &(lp->srhs));
      ILLfct_compute_dpIIy (lp, &(lp->srhs), &(lp->ssoln));
      alpha1 = (x_bi - v_l - lp->upd.dty) / rs.pivotval;
      ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);
      rval = ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
      ILL_CLEANUP_IF (rval);

      ILLfct_update_dpII_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), rs.eindex, pr.lindex, eval, alpha1);

#if DENSE_NORM > 0
      test_dsteep_norms (lp, pinf);
#endif

   END:
      if (singular){
         it->nextstep = SIMPLEX_RESUME;
         it->resumeid = SIMPLEX_RESUME_SING;
         ILL_CLEANUP;
      }
      if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE){
         ILLfct_compute_piz (lp);
         ILLfct_compute_dz (lp);
         ILLfct_dual_adjust (lp, 0.0);
         it->newphase = SIMPLEX_PHASE_RECOMP;
      }
   }

#if DENSE_PIIPI > 0
   fct_test_workvector (lp);
   if (!refactor){
      fct_test_pII_x (lp, pinf);
      fct_test_pII_pi_dz (lp, pinf);
   }
#endif

 CLEANUP:
   ILL_RETURN (rval, "dual_phaseII_step");
}

static void get_current_stat (lp_status_info *p, int algorithm, int *bstat)
{
   if (p->optimal) *bstat = OPTIMAL;
   else if (algorithm == PRIMAL_SIMPLEX){
      if (p->primal_feasible) *bstat = PRIMAL_FEASIBLE;
      else if (p->primal_infeasible) *bstat = PRIMAL_INFEASIBLE;
      else if (p->primal_unbounded) *bstat = PRIMAL_UNBOUNDED;
      else *bstat = NONOPTIMAL;
   }
   else if (algorithm == DUAL_SIMPLEX){
      if (p->dual_feasible) *bstat = DUAL_FEASIBLE;
      else if (p->dual_infeasible) *bstat = DUAL_INFEASIBLE;
      else if (p->dual_unbounded) *bstat = DUAL_UNBOUNDED;
      else *bstat = NONOPTIMAL;
   }
}

int ILLsimplex_pivotin (lpinfo *lp, price_info *pinf, int rcnt, int *rlist,
                        int pivot_opt, int *basis_mod)
{
   int  i, npiv  = 0;
   int  eindex;
   int  rval     = 0;
   int  singular = 0;
   int  refactor = 0;
   int  *rowmap  = lp->O->rowmap;
   int  *clist = NULL;
   double  alpha = 0.0;
   svector  wz;
   svector  updz;
   ratio_res  rs;
   feas_info  fi;

   *basis_mod = 0;
   if (rcnt <= 0) { ILL_RETURN (rval, "ILLsimplex_pivotin"); } 

   if (pivot_opt == SIMPLEX_PIVOTINROW){
      ILL_SAFE_MALLOC (clist, rcnt, int);
      for (i=0; i<rcnt; i++) clist[i] = rowmap[rlist[i]];
   }
   else
      clist = rlist;

   for (i=0; i<rcnt; i++){
      if (lp->vstat[clist[i]] != STAT_BASIC){
         *basis_mod = 1;
         break;
      }
   }
   if (*basis_mod == 0) {
       if (pivot_opt == SIMPLEX_PIVOTINROW) {
           ILL_IFFREE (clist, int);
       }
       ILL_RETURN (rval, "ILLsimplex_pivotin");
   }

   /* printf ("Forcing vars into basis in ILLsimplex_pivotin \n"); */
   ILLsvector_init (&wz);
   rval = ILLsvector_alloc (&wz, lp->nrows);
   ILL_CLEANUP_IF (rval);
   ILLsvector_init (&updz);
   rval = ILLsvector_alloc (&updz, lp->nrows);
   ILL_CLEANUP_IF (rval);

   lp->pobjval = lp->dobjval;
   for (i=0; i<rcnt; i++){
      if (lp->vstat[clist[i]] == STAT_BASIC)
         continue;
      npiv ++;

      eindex = lp->vindex[clist[i]];
      ILLfct_compute_yz (lp, &(lp->yjz), &updz, lp->nbaz[eindex]);
      ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, 0.0);
      ILLfct_update_counts (lp, CNT_UPNZ, updz.nzcnt, 0.0);

      ILLratio_pivotin_test (lp, clist, rcnt, &rs);

      if (rs.ratio_stat == RATIO_UNBOUNDED || rs.ratio_stat == RATIO_FAILED){
         fprintf (stderr, "Pivot_in failed\n");
         rval = E_SIMPLEX_ERROR;
         ILL_CLEANUP;
      }
      else if (rs. ratio_stat == RATIO_BCHANGE){
         alpha = lp->dz[eindex] / rs.pivotval;
         if (rs.lvstat == STAT_LOWER)
            lp->dobjval += rs.tz * (lp->lz[lp->baz[rs.lindex]] - lp->xbz[rs.lindex]);
         else
            lp->dobjval += rs.tz * (lp->xbz[rs.lindex] - lp->uz[lp->baz[rs.lindex]]);
         lp->objval = lp->dobjval;

         ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
         ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
         ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, 0.0);
         ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, 0.0);

         if (pinf->dsinfo.norms && pinf->dII_price == QS_PRICE_DSTEEP){
            ILLfct_compute_dsteep_upv (lp, &wz);
            rval = ILLprice_update_pricing_info (lp, pinf, DUAL_PHASEII, &wz,
                                               eindex, rs.lindex, rs.pivotval);
            ILL_CLEANUP_IF (rval);
         }
         else if (pinf->psinfo.norms && pinf->pII_price == QS_PRICE_PSTEEP){
            ILLfct_compute_psteep_upv (lp, &wz);
            rval = ILLprice_update_pricing_info (lp, pinf, PRIMAL_PHASEII, &wz,
                                               eindex, rs.lindex, rs.pivotval);
            ILL_CLEANUP_IF (rval);
         }

         ILLfct_update_xz (lp, rs.tz, eindex, rs.lindex);
         ILLfct_update_basis_info (lp, eindex, rs.lindex, rs.lvstat);
         rval = ILLbasis_update (lp, &updz, rs.lindex, &refactor, &singular);
         ILL_CLEANUP_IF (rval);

         if (singular){
            fprintf (stderr, "singular matrix in pivot_in\n");
            rval = E_SIMPLEX_ERROR;
            ILL_CLEANUP;
         }
         if (!refactor){
            ILLfct_update_piz (lp, alpha);
            ILLfct_update_dz (lp, eindex, alpha);
         }
         else{
            ILLfct_compute_xbz (lp);
            ILLfct_compute_piz (lp);
            ILLfct_compute_dz (lp);
            ILLfct_compute_dobj (lp);
         }
      }
   }
   /*
   ILLfct_dphaseI_simple_update (lp, lp->tol->dfeas_tol);
   ILLfct_compute_xbz (lp);
   ILLfct_compute_dobj (lp);
   */

   ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
   ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
   ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

 CLEANUP:
   if (pivot_opt == SIMPLEX_PIVOTINROW)
      ILL_IFFREE (clist, int);
   ILLsvector_free (&wz);
   ILLsvector_free (&updz);
   ILL_RETURN (rval, "ILLsimplex_pivotin");
}

static int report_value(lpinfo *lp, iter_info *it, 
                        const char* value_name, double value) 
{
	int rval = 0;
	char buffer[256];
	
	if (it->sdisplay && it->itercnt % lp->iterskip == 0) {
		sprintf(buffer, "(%d): %s = %f\n", it->itercnt, value_name, value); 
		rval = ILLstring_report(buffer, &lp->O->reporter); 
		fflush (stdout);
	} else { 
		/* make sure ILLstring_report is called at least every 10 iterations */
		if (it->itercnt % (lp->iterskip / 10)) {
			rval = ILLstring_report(NULL, &lp->O->reporter); 
		} 
	}
	if (rval != 0) { /* ILLstring_report was called and failed, which means we should abort */
		it->solstatus = QS_LP_ABORTED; 
	} 
	return rval;
}
