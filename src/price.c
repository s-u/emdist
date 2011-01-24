/* RCS_INFO = "$RCSfile: price.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#include "stddefs.h"
#include "qsopt.h"
#include "lpdefs.h"
#include "fct.h"
#include "price.h"
#include "basis.h"
#include "iqsutil.h"
#include "dstruct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define  MULTIP 1
#define  PRICE_DEBUG 0

static void
    update_d_scaleinf (price_info *p, heap *h, int j, double inf, int prule),
    update_p_scaleinf (price_info *p, heap *h, int i, double inf, int prule);

static double
    compute_dualI_inf (lpinfo *lp, int j),
    compute_dualII_inf (lpinfo *lp, int j),
    compute_primalI_inf (lpinfo *lp, int i),
    compute_primalII_inf (lpinfo *lp, int i);

void ILLprice_free_heap (price_info *pinf)
{
   ILLheap_free (&(pinf->h));
}

int ILLprice_build_heap (price_info *pinf, int nkeys, double *keylist)
{
   ILLheap_init (&(pinf->h));
   pinf->htrigger = 1.0 + (double)nkeys / (PARAM_HEAP_RATIO * ILLutil_our_log2 (nkeys));
   return ILLheap_build (&(pinf->h), nkeys, keylist);
}

int ILLprice_test_for_heap (lpinfo *lp, price_info *pinf, int nkeys,
                            double *keylist, int algo, int upd)
{
   int  rval = 0;
   heap   *h = &(pinf->h);
   double  ravg = 0.0;

   ravg = (algo == PRIMAL_SIMPLEX) ? lp->cnts->za_ravg : lp->cnts->y_ravg;
   if (upd != 0){
      if (ravg <= pinf->htrigger) pinf->hineff --;
      else if (ravg > 2.0*pinf->htrigger) pinf->hineff ++;
   }
   if (h->hexist == 0 && pinf->hineff <= 0){
      rval = ILLprice_build_heap (pinf, nkeys, keylist);
      ILL_CLEANUP_IF (rval);
   }
   else if (h->hexist != 0 && pinf->hineff >= PARAM_HEAP_UTRIGGER){
      ILLprice_free_heap (pinf);
      /*
      printf ("freeing heap ..\n");
      printf ("iter = %d, ravg = %.2f, trigger = %.2f\n",
              lp->cnts->tot_iter, ravg, pinf->htrigger);
      */
   }

 CLEANUP:
   if (rval) ILLprice_free_heap (pinf);
   return rval;
}

void ILLprice_init_pricing_info (price_info *pinf)
{
   pinf->p_strategy = -1;
   pinf->d_strategy = -1;
   pinf->pI_price   = -1;
   pinf->pII_price  = -1;
   pinf->dI_price   = -1;
   pinf->dII_price  = -1;
   pinf->cur_price  = -1;
   pinf->p_scaleinf = (double *) NULL;
   pinf->d_scaleinf = (double *) NULL;
   pinf->pdinfo.norms    = (double *) NULL;
   pinf->pdinfo.refframe = (int *) NULL;
   pinf->psinfo.norms    = (double *) NULL;
   pinf->ddinfo.norms    = (double *) NULL;
   pinf->ddinfo.refframe = (int *) NULL;
   pinf->dsinfo.norms    = (double *) NULL;
   pinf->dmpinfo.gstart  = pinf->pmpinfo.gstart  = (int *) NULL;
   pinf->dmpinfo.gshift  = pinf->pmpinfo.gshift  = (int *) NULL;
   pinf->dmpinfo.gsize   = pinf->pmpinfo.gsize   = (int *) NULL;
   pinf->dmpinfo.bucket  = pinf->pmpinfo.bucket  = (int *) NULL;
   pinf->dmpinfo.perm    = pinf->pmpinfo.perm    = (int *) NULL;
   pinf->dmpinfo.infeas  = pinf->pmpinfo.infeas  = (double *) NULL;
   ILLheap_init (&(pinf->h));
   pinf->htrigger = 0.0;
   pinf->hineff   = 0;
}

void ILLprice_free_pricing_info (price_info *pinf)
{
   ILL_IFFREE (pinf->p_scaleinf, double);
   ILL_IFFREE (pinf->d_scaleinf, double);
   ILL_IFFREE (pinf->pdinfo.norms, double);
   ILL_IFFREE (pinf->pdinfo.refframe, int);
   ILL_IFFREE (pinf->psinfo.norms, double);
   ILL_IFFREE (pinf->ddinfo.norms, double);
   ILL_IFFREE (pinf->ddinfo.refframe, int);
   ILL_IFFREE (pinf->dsinfo.norms, double);
   ILLprice_free_mpartial_info (&(pinf->pmpinfo));
   ILLprice_free_mpartial_info (&(pinf->dmpinfo));
   ILLprice_free_heap (pinf);
}

int ILLprice_build_pricing_info (lpinfo *lp, price_info *pinf, int phase)
{
   int  rval = 0;
   int  p_price = -1;
   int  d_price = -1;

   switch (phase){
   case PRIMAL_PHASEI:  p_price = pinf->pI_price; break;
   case PRIMAL_PHASEII: p_price = pinf->pII_price; break;
   case DUAL_PHASEI:    d_price = pinf->dI_price; break;
   case DUAL_PHASEII:   d_price = pinf->dII_price; break;
   }

   if (p_price != -1){
      pinf->cur_price = p_price;

      if (p_price == QS_PRICE_PDANTZIG || p_price == QS_PRICE_PDEVEX ||
          p_price == QS_PRICE_PSTEEP){
         pinf->p_strategy = COMPLETE_PRICING;
         ILL_IFFREE (pinf->d_scaleinf, double);
         ILL_SAFE_MALLOC (pinf->d_scaleinf, lp->nnbasic, double);
      }
      else if (p_price == QS_PRICE_PMULTPARTIAL)
         pinf->p_strategy = MULTI_PART_PRICING;

      switch(p_price){
      case QS_PRICE_PDEVEX:
         if (pinf->pdinfo.norms) return rval;
         rval = ILLprice_build_pdevex_norms (lp, &(pinf->pdinfo), 0);
         ILL_CLEANUP_IF (rval);
         break;
      case QS_PRICE_PSTEEP:
         if (pinf->psinfo.norms) return rval;
         rval = ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
         ILL_CLEANUP_IF (rval);
         break;
      case QS_PRICE_PMULTPARTIAL:
         rval = ILLprice_build_mpartial_info (lp, pinf, COL_PRICING);
         ILL_CLEANUP_IF (rval);
         break;
      }
   }
   else if (d_price != -1){
      pinf->cur_price = d_price;

      if (d_price == QS_PRICE_DDANTZIG || d_price == QS_PRICE_DSTEEP ||
          d_price == QS_PRICE_DDEVEX){
         pinf->d_strategy  = COMPLETE_PRICING;
         ILL_IFFREE (pinf->p_scaleinf, double);
         ILL_SAFE_MALLOC (pinf->p_scaleinf, lp->nrows, double);
      }
      else if (d_price == QS_PRICE_DMULTPARTIAL)
         pinf->d_strategy = MULTI_PART_PRICING;

      switch(d_price){
      case QS_PRICE_DSTEEP:
         if (pinf->dsinfo.norms) return rval;
         rval = ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
         ILL_CLEANUP_IF (rval);
         break;
      case QS_PRICE_DMULTPARTIAL:
         rval = ILLprice_build_mpartial_info (lp, pinf, ROW_PRICING);
         ILL_CLEANUP_IF (rval);
         break;
      case QS_PRICE_DDEVEX:
         if (pinf->ddinfo.norms) return rval;
         rval = ILLprice_build_ddevex_norms (lp, &(pinf->ddinfo), 0);
         ILL_CLEANUP_IF (rval);
         break;
      }
   }

 CLEANUP:
   if (rval) ILLprice_free_pricing_info (pinf);
   ILL_RETURN (rval, "ILLprice_build_pricing_info");
}

int ILLprice_update_pricing_info (lpinfo *lp, price_info *pinf, int phase,
                                 svector *wz, int eindex, int lindex, double y)
{
   int  rval = 0;
   int  p_price = -1;
   int  d_price = -1;

   switch (phase){
   case PRIMAL_PHASEI:  p_price = pinf->pI_price; break;
   case PRIMAL_PHASEII: p_price = pinf->pII_price; break;
   case DUAL_PHASEI:    d_price = pinf->dI_price; break;
   case DUAL_PHASEII:   d_price = pinf->dII_price; break;
   }

   if (p_price != -1){
      if (p_price == QS_PRICE_PDEVEX){
         rval = ILLprice_update_pdevex_norms (lp, &(pinf->pdinfo), eindex, y);
         ILL_CLEANUP_IF (rval);
      }
      else if (p_price == QS_PRICE_PSTEEP)
         ILLprice_update_psteep_norms (lp, &(pinf->psinfo), wz, eindex, y);
   }
   else if (d_price != -1){
      if (d_price == QS_PRICE_DSTEEP)
         ILLprice_update_dsteep_norms (lp, &(pinf->dsinfo), wz, lindex, y);
      else if (d_price == QS_PRICE_DDEVEX){
         rval = ILLprice_update_ddevex_norms (lp, &(pinf->ddinfo), lindex, y);
         ILL_CLEANUP_IF (rval);
      }
   }
 CLEANUP:
   ILL_RETURN (rval, "ILLprice_update_pricing_info");
}

int ILLprice_get_price (price_info *p, int phase)
{
   int  pri = -1;

   switch (phase){
   case PRIMAL_PHASEI:  return p->pI_price;
   case PRIMAL_PHASEII: return p->pII_price;
   case DUAL_PHASEI:    return p->dI_price;
   case DUAL_PHASEII:   return p->dII_price;
   }
   return pri;
}

void ILLprice_free_mpartial_info (mpart_info *p)
{
   ILL_IFFREE (p->gstart, int);
   ILL_IFFREE (p->gshift, int);
   ILL_IFFREE (p->gsize, int);
   ILL_IFFREE (p->bucket, int);
   ILL_IFFREE (p->infeas, double);
   ILL_IFFREE (p->perm, int);
}

int ILLprice_build_mpartial_info (lpinfo *lp, price_info *pinf, int pricetype)
{
   int  i = 0;
   int  rval = 0;
   int  extra = 0;
   int  nelems;
   mpart_info  *p;

   p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
   p->k = 50;
   p->cgroup  = 0;
   nelems = (pricetype == COL_PRICING) ? lp->nnbasic : lp->nrows;

   if (nelems % p->k)
      extra = nelems - p->k * (nelems / p->k);
   p->ngroups = nelems / p->k;
   if (extra != 0) p->ngroups ++;

   ILL_SAFE_MALLOC (p->gstart, p->ngroups, int);
   ILL_SAFE_MALLOC (p->gshift, p->ngroups, int);
   ILL_SAFE_MALLOC (p->gsize, p->ngroups, int);
   ILL_SAFE_MALLOC (p->bucket, 2 * p->k, int);
   ILL_SAFE_MALLOC (p->infeas, 2 * p->k, double);
   ILL_SAFE_MALLOC (p->perm, 2 * p->k, int);

   p->bsize = 0;

   if (extra != 0){
      p->gstart[0] = 0;
      p->gshift[0] = 1;
      p->gsize[0]  = extra;
      for (i=1; i<p->ngroups; i++){
         p->gstart[i] = extra + i-1;
         p->gshift[i] = p->ngroups - 1;
         p->gsize[i]  = p->k;
      }
   }
   else{
      for (i=0; i<p->ngroups; i++){
         p->gstart[i] = i;
         p->gshift[i] = p->ngroups;
         p->gsize[i]  = p->k;
      }
   }

 CLEANUP:
   if (rval) ILLprice_free_mpartial_info (p);
   ILL_RETURN (rval, "ILLprice_build_mpartial_info");
}

void ILLprice_init_mpartial_price (lpinfo *lp, price_info *pinf, int phase, int pricetype)
{
   int  i;
   mpart_info  *p;

   p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
   p->bsize = 0;
   i = p->cgroup;
   do{
      ILLprice_mpartial_group (lp, p, phase, i, pricetype);
      i = (i+1) % p->ngroups;
   } while (i != p->cgroup && p->bsize <= p->k);
   p->cgroup = i;
}

void ILLprice_update_mpartial_price (lpinfo *lp, price_info *pinf,
                                     int phase, int pricetype)
{
   int  i = 0;
   int  csize = 0;
   double  infeas;
   mpart_info  *p;
   price_res   pr;

   p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);

#ifdef MULTIP
   i = 0;
   while (i < p->bsize){
      if (pricetype == COL_PRICING){
         ILLprice_column (lp, p->bucket[i], phase, &pr);
         infeas = pr.dinfeas;
      }
      else{
         ILLprice_row (lp, p->bucket[i], phase, &pr);
         infeas = pr.pinfeas;
      }
      if (infeas == 0.0){
         p->bucket[i] = p->bucket[p->bsize-1];
         p->bsize --;
      }
      else{
         p->infeas[i] = infeas;
         i++;
      }
   }
   if (p->bsize > 0){
      for (i=0; i<p->bsize; i++) p->perm[i] = i;
      ILLutil_double_perm_quicksort(p->perm, p->infeas, p->bsize);

      csize = MIN (p->bsize, p->k);
      for (i=csize-1; i>=0; i--) lp->iwork[p->bucket[p->perm[i]]] = 1;

      for (i=0, csize=0; i<p->bsize; i++)
         if (lp->iwork[p->bucket[i]] == 1){
            p->infeas[csize] = p->infeas[i];
            p->bucket[csize] = p->bucket[i];
            csize ++;
         }
      p->bsize = csize;
   }
#else
   p->bsize = 0;
#endif

   i = p->cgroup;
   do{
      ILLprice_mpartial_group (lp, p, phase, i, pricetype);
      i = (i+1) % p->ngroups;
   } while (i != p->cgroup && p->bsize <= p->k);
   p->cgroup = i;

#ifdef MULTIP
   for (i=0; i<csize; i++) lp->iwork[p->bucket[i]] = 0;
#endif
}

void ILLprice_delete_onempart_price (lpinfo *lp, price_info *pinf,
                                     int indx, int pricetype)
{
   int  i = 0;
   mpart_info  *p;

   p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);

   for (i=0; i<p->bsize; i++)
      if (p->bucket[i] == indx){
         p->bucket[i] = p->bucket[p->bsize-1];
         p->infeas[i] = p->infeas[p->bsize-1];
         p->bsize --;
         break;
      }
}

void ILLprice_mpartial_group (lpinfo *lp, mpart_info *p, int phase, int g, int pricetype)
{
   int  i, ix;
   int  gstart = p->gstart[g];
   int  gsize  = p->gsize[g];
   int  gshift = p->gshift[g];
   double  infeas;
   price_res  pr;

   for (i=0, ix=gstart; i<gsize; i++, ix+=gshift){
#ifdef MULTIP
      if (lp->iwork[ix]) continue;
#endif
      if (pricetype == COL_PRICING){
         ILLprice_column (lp, ix, phase, &pr);
         infeas = pr.dinfeas;
      }
      else{
         ILLprice_row (lp, ix, phase, &pr);
         infeas = pr.pinfeas;
      }
      if (infeas != 0.0){
         p->infeas[p->bsize] = infeas;
         p->bucket[p->bsize] = ix;
         p->bsize ++;
      }
   }
}

void ILLprice_column (lpinfo *lp, int ix, int phase, price_res *pr)
{
   int     i;
   int     col;
   int     mcnt;
   int     mbeg;
   double  sum;

   pr->dinfeas = 0.0;
   col  = lp->nbaz[ix];
   if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
      return;
   sum  = 0.0;
   mcnt = lp->matcnt[col];
   mbeg = lp->matbeg[col];

   if (phase == PRIMAL_PHASEII){
      for (i=0; i<mcnt; i++)
         sum += lp->piz[lp->matind[mbeg+i]] * lp->matval[mbeg+i];
      lp->dz[ix]  = lp->cz[col] - sum;
      pr->dinfeas = compute_dualII_inf (lp, ix);
   }
   else{
      for (i=0; i<mcnt; i++)
         sum += lp->pIpiz[lp->matind[mbeg+i]] * lp->matval[mbeg+i];
      lp->pIdz[ix]  = -sum;
      pr->dinfeas = compute_dualI_inf (lp, ix);
   }
}

void ILLprice_row (lpinfo *lp, int ix, int phase, price_res *pr)
{
   if (phase == DUAL_PHASEII)
      pr->pinfeas = compute_primalII_inf (lp, ix);
   else
      pr->pinfeas = compute_primalI_inf (lp, ix);
}

int ILLprice_build_pdevex_norms (lpinfo *lp, p_devex_info *pdinfo, int reinit)
{
   int  j;
   int  rval = 0;

   if (reinit == 0){
      pdinfo->ninit = 0;
      ILL_SAFE_MALLOC (pdinfo->norms, lp->nnbasic, double);
      ILL_SAFE_MALLOC (pdinfo->refframe, lp->ncols, int);
   }

   if (reinit != 0) pdinfo->ninit ++;

   for (j=0; j<lp->ncols; j++){
      if (lp->vstat[j] == STAT_BASIC)
         pdinfo->refframe[j] = 0;
      else{
         pdinfo->norms[lp->vindex[j]] = 1.0;
         pdinfo->refframe[j] = 1;
      }
   }

 CLEANUP:
   if (rval){
      ILL_IFFREE (pdinfo->norms, double);
      ILL_IFFREE (pdinfo->refframe, int);
   }
   ILL_RETURN (rval, "ILLprice_build_pdevex_norms");
}

int ILLprice_update_pdevex_norms (lpinfo *lp, p_devex_info *pdinfo, int eindex,
                                  double yl)
{
   int  i, j;
   double  normj = 0.0;
   double  yl_sqr = yl*yl;
   double  zAj;

   for (i=0; i<lp->yjz.nzcnt; i++)
      if (pdinfo->refframe[lp->baz[lp->yjz.indx[i]]])
         normj += lp->yjz.coef[i] * lp->yjz.coef[i];

   if (pdinfo->refframe[lp->nbaz[eindex]])
      normj += 1.0;

   if ((normj > 1000.0*pdinfo->norms[eindex]) ||
       (normj < 0.001*pdinfo->norms[eindex]) )
      return ILLprice_build_pdevex_norms (lp, pdinfo, 1);

   for (i=0; i<lp->zA.nzcnt; i++){
      j   = lp->zA.indx[i];
      zAj = lp->zA.coef[i];
      pdinfo->norms[j] = MAX (pdinfo->norms[j], (normj*zAj*zAj / yl_sqr));
   }

   pdinfo->norms[eindex] = MAX (1.0, normj / yl_sqr);
   return 0;
}

int ILLprice_build_psteep_norms (lpinfo *lp, p_steep_info *psinfo)
{
   int  j;
   int  rval = 0;
   svector  yz;

   ILLsvector_init (&yz);
   rval = ILLsvector_alloc (&yz, lp->nrows);
   ILL_CLEANUP_IF (rval);
   ILL_SAFE_MALLOC (psinfo->norms, lp->nnbasic, double);

   for (j=0; j<lp->nnbasic; j++){
	  rval = ILLstring_report(NULL, &lp->O->reporter); 
	  ILL_CLEANUP_IF(rval); 
      ILLfct_compute_yz (lp, &yz, (svector *) NULL, lp->nbaz[j]);
      psinfo->norms[j] = 1.0 + ILLutil_norm_sqr (yz.coef, yz.nzcnt);
   }

 CLEANUP:
   ILLsvector_free (&yz);
   if (rval) ILL_IFFREE (psinfo->norms, double);
   ILL_RETURN (rval, "ILLprice_build_psteep_norms");
}

void ILLprice_update_psteep_norms (lpinfo *lp, p_steep_info *psinfo,
                                   svector *wz, int eindex, double yl)
{
   int  i, j, k;
   int  mcnt, mbeg;
   double  normj;
   double  zAj, wAj;
   double *v;

   normj = 1.0 + ILLutil_norm_sqr (lp->yjz.coef, lp->yjz.nzcnt);

#if 0
   Bico - remove warnings for dist 
   if (fabs ((normj - psinfo->norms[eindex]) / normj) > 1000.0 /* 0.01 */){
      printf ("warning: incorrect norm values\n");
      printf ("anorm = %.6f, pnorm = %.6f\n", normj, psinfo->norms[eindex]);
      fflush (stdout);
   }
#endif

   ILLfct_load_workvector (lp, wz);
   v = lp->work.coef;

   for (k=0; k<lp->zA.nzcnt; k++){
      j   = lp->zA.indx[k];
      zAj = lp->zA.coef[k];
      wAj  = 0.0;
      mcnt = lp->matcnt[lp->nbaz[j]];
      mbeg = lp->matbeg[lp->nbaz[j]];
      for (i=0; i<mcnt; i++)
         wAj += lp->matval[mbeg+i] * v[lp->matind[mbeg+i]];

      psinfo->norms[j] += (zAj * ((zAj * normj / yl)-(2.0 * wAj)))/yl;
      if (psinfo->norms[j] < 1.0)
         psinfo->norms[j] = 1.0;
   }

   psinfo->norms[eindex] = normj / (yl*yl);
   if (psinfo->norms[eindex] < 1.0)
      psinfo->norms[eindex] = 1.0;
   ILLfct_zero_workvector (lp);
}

int ILLprice_build_ddevex_norms (lpinfo *lp, d_devex_info *ddinfo, int reinit)
{
   int  i;
   int  rval = 0;

   if (reinit == 0){
      ddinfo->ninit = 0;
      ILL_SAFE_MALLOC (ddinfo->norms, lp->nrows, double);
      ILL_SAFE_MALLOC (ddinfo->refframe, lp->ncols, int);
   }
   if (reinit != 0) ddinfo->ninit ++;

   for (i=0; i<lp->ncols; i++)
      ddinfo->refframe[i] = (lp->vstat[i] == STAT_BASIC) ? 1 : 0;

   for (i=0; i<lp->nrows; i++) ddinfo->norms[i] = 1.0;

 CLEANUP:
   if (rval){
      ILL_IFFREE (ddinfo->norms, double);
      ILL_IFFREE (ddinfo->refframe, int);
   }
   ILL_RETURN (rval, "ILLprice_build_ddevex_norms");
}

int ILLprice_update_ddevex_norms (lpinfo *lp, d_devex_info *ddinfo, int lindex,
                                  double yl)
{
   int  i, r;
   double  normi = 0.0;
   double  yl_sqr = yl*yl;
   double  yr;

   for (i=0; i<lp->zA.nzcnt; i++)
      if (ddinfo->refframe[lp->nbaz[lp->zA.indx[i]]])
         normi += lp->zA.coef[i] * lp->zA.coef[i];

   if (ddinfo->refframe[lp->baz[lindex]])
      normi += 1.0;

   if ((normi > 1000.0*ddinfo->norms[lindex]) ||
       (normi < 0.001*ddinfo->norms[lindex]) )
      return ILLprice_build_ddevex_norms (lp, ddinfo, 1);

   for (i=0; i<lp->yjz.nzcnt; i++){
      r   = lp->yjz.indx[i];
      yr  = lp->yjz.coef[i];
      ddinfo->norms[r] = MAX (ddinfo->norms[r], (normi*yr*yr / yl_sqr));
   }

   ddinfo->norms[lindex] = MAX (1.0, normi / yl_sqr);
   return 0;
}

int ILLprice_build_dsteep_norms (lpinfo *lp, d_steep_info *dsinfo)
{
   int  i;
   int  rval = 0;
   svector z;

   ILLsvector_init (&z);
   rval = ILLsvector_alloc (&z, lp->nrows);
   ILL_CLEANUP_IF (rval);
   ILL_SAFE_MALLOC (dsinfo->norms, lp->nrows, double);

   for (i=0; i<lp->nrows; i++){
	  rval = ILLstring_report(NULL, &lp->O->reporter); 
	  ILL_CLEANUP_IF (rval);

      ILLfct_compute_zz (lp, &z, i);

      dsinfo->norms[i] = ILLutil_norm_sqr (z.coef, z.nzcnt);
      if (dsinfo->norms[i] < PARAM_MIN_DNORM)
         dsinfo->norms[i] = PARAM_MIN_DNORM;
   }

 CLEANUP:
   ILLsvector_free (&z);
   if (rval) ILL_IFFREE (dsinfo->norms, double);
   ILL_RETURN (rval, "ILLprice_build_dsteep_norms");
}

int ILLprice_get_dsteep_norms (lpinfo *lp, int count, int *rowind, double *norms)
{
   int  i;
   int  rval = 0;
   svector z;

   ILLsvector_init (&z);
   rval = ILLsvector_alloc (&z, lp->nrows);
   ILL_CLEANUP_IF (rval);

   for (i = 0; i < count; i++) {
      ILLfct_compute_zz (lp, &z, rowind[i]);
      norms[i] = ILLutil_norm_sqr (z.coef, z.nzcnt);
   }

 CLEANUP:
   ILLsvector_free (&z);
   ILL_RETURN (rval, "ILLprice_get_dsteep_norms");
}

void ILLprice_update_dsteep_norms (lpinfo *lp, d_steep_info *dsinfo,
                                   svector *wz, int lindex, double yl)
{
   int  i, k;
   double  yij;
   double  norml = 0.0;
   double  *v;

   norml = ILLutil_norm_sqr (lp->zz.coef, lp->zz.nzcnt);

#if 0
   Bico - remove warnings for dist 
   if (fabs ((norml - dsinfo->norms[lindex]) / norml) > 1000.0 /*0.01*/ ){
      printf ("warning: incorrect dnorm values\n");
      printf ("anorm = %.6f, pnorm = %.6f\n", norml, dsinfo->norms[lindex]);
      fflush (stdout);
   }
#endif

   ILLfct_load_workvector (lp, wz);
   v = lp->work.coef;

   for (k=0; k<lp->yjz.nzcnt; k++){
      i   = lp->yjz.indx[k];
      yij = lp->yjz.coef[k];
      dsinfo->norms[i] += (yij * ((yij * norml / yl) - (2.0 * v[i]))) / yl;
      if (dsinfo->norms[i] < PARAM_MIN_DNORM)
         dsinfo->norms[i] = PARAM_MIN_DNORM;
   }
   dsinfo->norms[lindex] = norml / (yl*yl);
   if (dsinfo->norms[lindex] < PARAM_MIN_DNORM)
      dsinfo->norms[lindex] = PARAM_MIN_DNORM;

   ILLfct_zero_workvector (lp);
}

static void update_d_scaleinf (price_info *p, heap *h, int j, double inf, int prule)
{
   if (inf == 0.0){
      p->d_scaleinf[j] = 0.0;
      if (h->hexist != 0 && h->loc[j] != -1) ILLheap_delete (h, j);
   }
   else {
      if (prule == QS_PRICE_PDANTZIG)
         p->d_scaleinf[j] = inf;
      else if (prule == QS_PRICE_PDEVEX)
         p->d_scaleinf[j] = (inf * inf) / p->pdinfo.norms[j];
      else if (prule == QS_PRICE_PSTEEP)
         p->d_scaleinf[j] = (inf * inf) / p->psinfo.norms[j];

      if (h->hexist != 0){
         if (h->loc[j] == -1) ILLheap_insert (h, j);
         else                 ILLheap_modify (h, j);
      }
   }
}

static double compute_dualI_inf (lpinfo *lp, int j)
{
   int     col  = lp->nbaz[j];
   int     vt   = lp->vtype[col];
   int     vs   = lp->vstat[col];
   double  dj   = lp->pIdz[j];
   double  inf  = 0.0;
   double  ftol = lp->tol->id_tol;

   if (vt != VARTIFICIAL && vt != VFIXED){
      if (dj < -ftol && (vs == STAT_LOWER || vs == STAT_ZERO))
         inf = -dj;
      else if (dj > ftol && (vs == STAT_UPPER || vs == STAT_ZERO))
         inf = dj;
   }
   return inf;
}

static double compute_dualII_inf (lpinfo *lp, int j)
{
   int     col  = lp->nbaz[j];
   int     vt   = lp->vtype[col];
   int     vs   = lp->vstat[col];
   double  dj   = lp->dz[j];
   double  inf  = 0.0;
   double  ftol = lp->tol->dfeas_tol;

   if (vt != VARTIFICIAL && vt != VFIXED){
      if (dj < -ftol && (vs == STAT_LOWER || vs == STAT_ZERO))
         inf = -dj;
      else if (dj > ftol && (vs == STAT_UPPER || vs == STAT_ZERO))
         inf = dj;
   }
   return inf;
}

void ILLprice_compute_dual_inf (lpinfo *lp, price_info *p, int *ix, int icnt, int phase)
{
   int     i;
   int     price;
   double  inf = 0.0;
   heap    *h  = &(p->h);

   price = (phase == PRIMAL_PHASEI) ? p->pI_price : p->pII_price;

   if (phase == PRIMAL_PHASEI){
      if (ix == NULL)
         for (i=0; i<lp->nnbasic; i++){
            inf = compute_dualI_inf (lp, i);
            update_d_scaleinf (p, h, i, inf, price);
         }
      else
         for (i=0; i<icnt; i++){
            inf = compute_dualI_inf (lp, ix[i]);
            update_d_scaleinf (p, h, ix[i], inf, price);
         }
   }
   else if (phase == PRIMAL_PHASEII){
      if (ix == NULL)
         for (i=0; i<lp->nnbasic; i++){
            inf = compute_dualII_inf (lp, i);
            update_d_scaleinf (p, h, i, inf, price);
         }
      else
         for (i=0; i<icnt; i++){
            inf = compute_dualII_inf (lp, ix[i]);
            update_d_scaleinf (p, h, ix[i], inf, price);
         }
   }
}

void ILLprice_primal (lpinfo *lp, price_info *pinf, price_res *pr, int phase)
{
   int  j, vs;
   double  d_e, d_max;
   double  ftol = lp->tol->dfeas_tol;
   heap *h = &(pinf->h);

   pr->eindex = -1;
   d_max      = 0.0;

#if USEHEAP > 0
   ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf, PRIMAL_SIMPLEX, 1);
#endif

   if (pinf->p_strategy == COMPLETE_PRICING){
      if (h->hexist){
         pr->eindex = ILLheap_findmin (h);
         if (pr->eindex != -1) ILLheap_delete (h, pr->eindex); 
      }
      else{
         for (j=0; j<lp->nnbasic; j++){
            if (d_max < pinf->d_scaleinf[j]){
               d_max      = pinf->d_scaleinf[j];
               pr->eindex = j;
            }
         }
      }
   }
   else if (pinf->p_strategy == MULTI_PART_PRICING){
      for (j=0; j<pinf->pmpinfo.bsize; j++){
         if (d_max < pinf->pmpinfo.infeas[j]){
            d_max      = pinf->pmpinfo.infeas[j];
            pr->eindex = pinf->pmpinfo.bucket[j];
         }
      }
   }

   if (pr->eindex < 0)
      pr->price_stat = PRICE_OPTIMAL;
   else{
      if (phase == PRIMAL_PHASEI) d_e = lp->pIdz[pr->eindex];
      else                        d_e = lp->dz[pr->eindex];
      vs = lp->vstat[lp->nbaz[pr->eindex]];

      pr->price_stat = PRICE_NONOPTIMAL;
      if (vs == STAT_UPPER || (vs == STAT_ZERO && d_e > ftol))
         pr->dir = VDECREASE;
      else
         pr->dir = VINCREASE;
   }
}

static void update_p_scaleinf (price_info *p, heap *h, int i, double inf, int prule)
{
   if (inf == 0.0){
      p->p_scaleinf[i] = 0.0;
      if (h->hexist != 0 && h->loc[i] != -1) ILLheap_delete (h, i);
   }
   else{
      if (prule == QS_PRICE_DDANTZIG)
         p->p_scaleinf[i] = inf;
      else if (prule == QS_PRICE_DSTEEP)
         p->p_scaleinf[i] = (inf * inf) / p->dsinfo.norms[i];
      else if (prule == QS_PRICE_DDEVEX)
         p->p_scaleinf[i] = (inf * inf) / p->ddinfo.norms[i];

      if (h->hexist != 0){
         if (h->loc[i] == -1) ILLheap_insert (h, i);
         else                 ILLheap_modify (h, i);
      }
   }
}

static double compute_primalI_inf (lpinfo *lp, int i)
{
   int     col  = lp->baz[i];
   double  x    = lp->xbz[i];
   double  l    = lp->lz[col];
   double  u    = lp->uz[col];
   double  inf  = 0.0;
   double  ftol = lp->tol->ip_tol;

   if (u != INFTY && x > ftol)
      inf = x;
   else if (l != NINFTY && x < - ftol)
      inf = -x;
   return inf;
}

static double compute_primalII_inf (lpinfo *lp, int i)
{
   int     col  = lp->baz[i];
   double  x    = lp->xbz[i];
   double  l    = lp->lz[col];
   double  u    = lp->uz[col];
   double  inf  = 0.0;
   double  ftol = lp->tol->pfeas_tol;

   if (u != INFTY && x > u + ftol)
      inf = x - u;
   else if (l != NINFTY && x < l - ftol)
      inf = l - x;
   return inf;
}

void ILLprice_compute_primal_inf (lpinfo *lp, price_info *p, int *ix, int icnt, int phase)
{
   int     i;
   int     price;
   double  inf = 0.0;
   heap    *h  = &(p->h);

   price = (phase == DUAL_PHASEI) ? p->dI_price : p->dII_price;

   if (phase == DUAL_PHASEI){
      if (ix == NULL)
         for (i=0; i<lp->nrows; i++){
            inf = compute_primalI_inf (lp, i);
            update_p_scaleinf (p, h, i, inf, price);
         }
      else
         for (i=0; i<icnt; i++){
            inf = compute_primalI_inf (lp, ix[i]);
            update_p_scaleinf (p, h, ix[i], inf, price);
         }
   }
   else if (phase == DUAL_PHASEII){
      if (ix == NULL)
         for (i=0; i<lp->nrows; i++){
            inf = compute_primalII_inf (lp, i);
            update_p_scaleinf (p, h, i, inf, price);
         }
      else
         for (i=0; i<icnt; i++){
            inf = compute_primalII_inf (lp, ix[i]);
            update_p_scaleinf (p, h, ix[i], inf, price);
         }
   }
}

void ILLprice_dual (lpinfo *lp, price_info *pinf, int phase, price_res *pr)
{
   int  i;
   double  p_max;
   double  ubound;
   double  ftol = lp->tol->pfeas_tol;
   heap  *h = &(pinf->h);

   pr->lindex = -1;
   p_max = 0.0;

#if USEHEAP > 0
   ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf, DUAL_SIMPLEX, 1);
#endif

   if (pinf->d_strategy == COMPLETE_PRICING){
      if (h->hexist){
         pr->lindex = ILLheap_findmin (h);
         if (pr->lindex != -1) ILLheap_delete (h, pr->lindex); 
      }
      else{
         for (i=0; i<lp->nrows; i++){
            if (p_max < pinf->p_scaleinf[i]){
               p_max = pinf->p_scaleinf[i];
               pr->lindex = i;
            }
         }
      }
   }
   else if (pinf->d_strategy == MULTI_PART_PRICING){
      for (i=0; i<pinf->dmpinfo.bsize; i++){
         if (p_max < pinf->dmpinfo.infeas[i]){
            p_max = pinf->dmpinfo.infeas[i];
            pr->lindex = pinf->dmpinfo.bucket[i];
         }
      }
   }

   if (pr->lindex < 0)
      pr->price_stat = PRICE_OPTIMAL;
   else{
      pr->price_stat = NONOPTIMAL;

      if (lp->uz[lp->baz[pr->lindex]] != INFTY){
         ubound = (phase == DUAL_PHASEI) ? 0.0 : lp->uz[lp->baz[pr->lindex]];
         if (lp->xbz[pr->lindex] > ubound + ftol)
            pr->lvstat = STAT_UPPER;
         else
            pr->lvstat = STAT_LOWER;
      }
      else
         pr->lvstat = STAT_LOWER;
   }
}

int ILLprice_get_rownorms (lpinfo *lp, price_info *pinf, double *rnorms)
{
   int  rval = 0;
   int  i;

   if (pinf->dsinfo.norms == NULL){
      rval = ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
      ILL_CLEANUP_IF (rval);
   }
   for (i=0; i<lp->nrows; i++)
      rnorms[i] = pinf->dsinfo.norms[i];

 CLEANUP:
   if (rval) ILL_IFFREE (pinf->dsinfo.norms, double);
   return rval;
}

int ILLprice_get_colnorms (lpinfo *lp, price_info *pinf, double *cnorms)
{
   int  rval = 0;
   int  i, j;

   if (pinf->psinfo.norms == NULL){
      rval = ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
      ILL_CLEANUP_IF (rval);
   }
   for (i=0; i<lp->nrows; i++)
      cnorms[lp->baz[i]] = 0.0;
   for (j=0; j<lp->nnbasic; j++)
      cnorms[lp->nbaz[j]] = pinf->psinfo.norms[j];

 CLEANUP:
   if (rval) ILL_IFFREE (pinf->psinfo.norms, double);
   return rval;
}

int ILLprice_get_newnorms (lpinfo *lp, int nelems, double *norms, int *matcnt,
                          int *matbeg, int *matind, double *matval, int option)
{
   int  i, j;
   int  rval = 0;
   double  n;
   svector a;
   svector y;

   ILLsvector_init (&y);
   rval = ILLsvector_alloc (&y, lp->nrows);
   ILL_CLEANUP_IF (rval);

   for (j=0; j<nelems; j++){
      a.nzcnt = matcnt[j];
      a.indx  = &(matind[matbeg[j]]);
      a.coef  = &(matval[matbeg[j]]);

      if (option == COLUMN_SOLVE)
         ILLbasis_column_solve (lp, &a, &y);
      else
         ILLbasis_row_solve (lp, &a, &y);

      n = 1.0;
      for (i=0; i<y.nzcnt; i++) n += y.coef[i] * y.coef[i];
      norms[j] = n;
   }

 CLEANUP:
   ILLsvector_free (&y);
   ILL_RETURN (rval, "ILLprice_get_newnorms");
}

int ILLprice_get_new_rownorms (lpinfo *lp, int newrows, double *rnorms,
                     int *rmatcnt, int *rmatbeg, int *rmatind, double *rmatval)
{
   return ILLprice_get_newnorms (lp, newrows, rnorms, rmatcnt, rmatbeg, rmatind, rmatval, ROW_SOLVE);
}

int ILLprice_get_new_colnorms (lpinfo *lp, int newrows, double *rnorms,
                         int *matcnt, int *matbeg, int *matind, double *matval)
{
   return ILLprice_get_newnorms (lp, newrows, rnorms, matcnt, matbeg, matind, matval, COLUMN_SOLVE);
}

int ILLprice_load_rownorms (lpinfo *lp, double *rnorms, price_info *pinf)
{
   int  i;
   int  rval = 0;

   ILL_IFFREE (pinf->dsinfo.norms, double);
   ILL_SAFE_MALLOC (pinf->dsinfo.norms, lp->nrows, double);
   for (i=0; i<lp->nrows; i++){
      pinf->dsinfo.norms[i] = rnorms[i];
      if (pinf->dsinfo.norms[i] < PARAM_MIN_DNORM)
         pinf->dsinfo.norms[i] = PARAM_MIN_DNORM;
   }

 CLEANUP:
   ILL_RETURN (rval, "ILLprice_load_rownorms");
}

int ILLprice_load_colnorms (lpinfo *lp, double *cnorms, price_info *pinf)
{
   int  j;
   int  rval = 0;

   ILL_IFFREE (pinf->psinfo.norms, double);
   ILL_SAFE_MALLOC (pinf->psinfo.norms, lp->nnbasic, double);
   for (j=0; j<lp->nnbasic; j++){
      pinf->psinfo.norms[j] = cnorms[lp->nbaz[j]];
      if (pinf->psinfo.norms[j] < 1.0) pinf->psinfo.norms[j] = 1.0;
   }

 CLEANUP:
   ILL_RETURN (rval, "ILLprice_load_colnorms");
}

#if PRICE_DEBUG > 0
void test_dsteep_norms (lpinfo *lp, price_info *p)
{
   int i, errn=0;
   double *pn = calloc (lp->nrows, sizeof (double));
   double err = 0.0;

   ILLprice_get_dsteep_norms (lp, lp->yjz.nzcnt,lp->yjz.indx, pn);
   for (i=0; i<lp->yjz.nzcnt; i++) if (fabs(pn[i] - p->dsinfo.norms[lp->yjz.indx[i]]) > 1e-5){
      errn ++; err += fabs(pn[i] - p->dsinfo.norms[lp->yjz.indx[i]]);
      p->dsinfo.norms[lp->yjz.indx[i]] = pn[i];  
   }
   if (errn > 0) printf ("%d: dnorm errn = %d, err = %.6f\n", lp->cnts->tot_iter,errn, err);
   free (pn);
}
#endif

