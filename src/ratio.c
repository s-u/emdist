/* "$RCSfile: ratio.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#include "stddefs.h"
#include "iqsutil.h"
#include "lpdefs.h"
#include "ratio.h"
#include "fct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void ILLratio_pI_test (lpinfo *lp, int eindex, int dir, ratio_res *rs) 
{
   int  i = 0, k = 0;
   int  col, ecol;
   int  cbnd, indx = 0;
   int  tctr = 0;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   double  t_i = 0.0, y_ij = 0.0, rcost;
   double  x, l, u;
   double  pivtol = lp->tol->pivot_tol;
   double  dftol  = lp->tol->id_tol; /*HHH*/
   double  *t = lp->upd.t;
   double  delta = 0.0; /*HHH*/

   rs->lindex     = -1;
   rs->tz         = 0.0;
   rs->ratio_stat = RATIO_FAILED;
   rs->pivotval   = 0.0;
   rs->lvstat     = -1;
   ecol           = lp->nbaz[eindex];

   if (lp->vtype[ecol] == VBOUNDED){
      t[0]  = lp->uz[ecol] - lp->lz[ecol];
      ix[0] = BBOUND;
      tctr ++;
   }

   for (k=0; k<lp->yjz.nzcnt; k++){
      y_ij = lp->yjz.coef[k];
      if (y_ij >= -pivtol && y_ij <= pivtol)
         continue;

      i   = lp->yjz.indx[k];
      x   = lp->xbz[i];
      col = lp->baz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      if ((dir == VINCREASE && y_ij > 0.0) ||
          (dir == VDECREASE && y_ij < 0.0)){
         if (y_ij < 0.0) y_ij = -y_ij;

         if (lp->bfeas[i] > 0){
            t[tctr]  = (x - u) / y_ij;
            ix[tctr] = 10*k + BATOUPPER;
            tctr ++;
            if (l != NINFTY){
               t[tctr]  = (x - l) / y_ij;
               ix[tctr] = 10*k + BATOLOWER;
               tctr ++;
            }
         }
         else if (lp->bfeas[i] == 0){
            if (l != NINFTY){
               t[tctr]  = (x - l) / y_ij;
               ix[tctr] = 10*k + BATOLOWER;
               tctr ++;
            }
         }
      }
      else if ((dir == VINCREASE && y_ij < 0.0) ||
               (dir == VDECREASE && y_ij > 0.0)){
         if (y_ij < 0.0) y_ij = -y_ij;

         if (lp->bfeas[i] < 0){
            t[tctr]  = (l - x) / y_ij;
            ix[tctr] = 10*k + BBTOLOWER;
            tctr ++;
            if (u != INFTY){
               t[tctr]  = (u - x) / y_ij;
               ix[tctr] = 10*k + BBTOUPPER;
               tctr ++;
            }
         }
         else if (lp->bfeas[i] == 0){
            if (u != INFTY){
               t[tctr]  = (u - x) / y_ij;
               ix[tctr] = 10*k + BBTOUPPER;
               tctr ++;
            }
         }
      }
   }
   if (tctr == 0){
      rs->ratio_stat = RATIO_FAILED;
      ILL_CLEANUP;
   }

   for (i=0; i<tctr; i++)  perm[i] = i;
   ILLutil_double_perm_quicksort (perm, t, tctr);

   lp->upd.c_obj = 0.0; /*HHH*/

   rcost = lp->pIdz[eindex];
   for (i=0; i<tctr; i++){
      t_i  = t[perm[i]];
      lp->upd.c_obj += (t_i-delta)*rcost; delta = t_i;/*HHH*/

      cbnd = ix[perm[i]] % 10;
      if (cbnd != BBOUND){
         k    = ix[perm[i]] / 10;
         y_ij = lp->yjz.coef[k];
         indx = lp->yjz.indx[k];
      }

      switch (cbnd){
      case BBOUND:
         rs->ratio_stat = RATIO_NOBCHANGE;
         rs->tz = (dir == VINCREASE) ? t_i : -t_i;
         ILL_CLEANUP;

      case BATOLOWER:
      case BATOUPPER: rcost += y_ij; break;
      case BBTOLOWER:
      case BBTOUPPER: rcost -= y_ij; break;
      }
      if ((dir == VINCREASE && rcost >= -dftol) ||
          (dir == VDECREASE && rcost <= dftol)){
         /* change 5 to -1 if t_i > 0 is required below */
         if (t_i < 0.0 && i > 5){
            /* printf ("pIhell %.5f %d\n", t_i, i); */
            t_i /= 2.0;
            rs->ratio_stat = RATIO_NEGATIVE;
            rs->tz = 0.0;
            ILL_CLEANUP;
         }
         rs->lindex = indx;
         rs->ratio_stat = RATIO_BCHANGE;
         if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
            rs->lvstat = STAT_LOWER;
         else
            rs->lvstat = STAT_UPPER;

         rs->pivotval = y_ij;
         rs->tz = (dir == VINCREASE) ? t_i : -t_i;
         ILL_CLEANUP;
      }
   }

 CLEANUP:
   ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
   lp->upd.tctr = tctr;
   lp->upd.i    = i;
   lp->upd.tz   = t_i;
   lp->upd.piv  = rs->pivotval;
   if (dir == VDECREASE)
      lp->upd.c_obj  = -lp->upd.c_obj;
   if (rs->lindex != -1)
      lp->upd.fs = lp->bfeas[rs->lindex];
}

void ILLratio_pII_test (lpinfo *lp, int eindex, int dir, ratio_res *rs)
{
   int  i, k, indx;
   int  col, ecol;
   double  x, l, u;
   double  y_ij, t_i, t_z;
   double  t_max, yi_max;
   double  pivtol = lp->tol->pivot_tol;
   double  pftol  = lp->tol->pfeas_tol; /*HHH*/

   rs->boundch    = 0;
   rs->lindex     = -1;
   rs->tz         = 0.0;
   rs->ratio_stat = RATIO_FAILED;
   rs->lvstat     = -1;
   rs->pivotval   = 0.0;
   rs->lbound     = 0.0;
   ecol           = lp->nbaz[eindex];

   for (k=0, t_max = INFTY; k<lp->yjz.nzcnt; k++){
      y_ij = lp->yjz.coef[k];
      if (y_ij >= -pivtol && y_ij <= pivtol)
         continue;

      t_i = INFTY;
      i   = lp->yjz.indx[k];
      x   = lp->xbz[i];
      col = lp->baz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      if ((dir == VINCREASE && y_ij > 0.0) || 
          (dir == VDECREASE && y_ij < 0.0)){
         if (l != NINFTY)
            t_i = (x - l + pftol) / fabs (y_ij);
      }
      else if ((dir == VINCREASE && y_ij < 0.0) ||
               (dir == VDECREASE && y_ij > 0.0)){
         if (u != INFTY)
            t_i = (u + pftol - x) / fabs (y_ij);
      }
      if (t_i == INFTY) continue;

      if (t_max > t_i){
         /*HHH tind = i; yval = fabs (y_ij); tval = t_i - pftol/fabs(y_ij);*/
         t_max  = t_i;
      }
   }

   if (lp->vtype[ecol] == VBOUNDED &&
       lp->uz[ecol] - lp->lz[ecol] <= t_max){

      t_max = lp->uz[ecol] - lp->lz[ecol];
      rs->ratio_stat = RATIO_NOBCHANGE;
      rs->tz = (dir == VINCREASE) ? t_max : -t_max;
      ILL_CLEANUP; 
   }

   if (t_max >= INFTY){
      rs->ratio_stat = RATIO_UNBOUNDED;
      ILL_CLEANUP; 
   }
   /* if (t_max < 0.0) printf ("pIIhell\n"); */

   indx   = -1;
   t_z    = 0.0;
   yi_max = 0.0;

   for (k=0; k<lp->yjz.nzcnt; k++){
      y_ij = lp->yjz.coef[k];
      if (y_ij >= -pivtol && y_ij <= pivtol)
         continue;

      t_i = INFTY;
      i   = lp->yjz.indx[k];
      x   = lp->xbz[i];
      col = lp->baz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      if ((dir == VINCREASE && y_ij > 0.0) ||
          (dir == VDECREASE && y_ij < 0.0)){
         if (l != NINFTY)
            t_i = (x - l) / fabs(y_ij);
      }
      else if ((dir == VINCREASE && y_ij < 0.0) ||
               (dir == VDECREASE && y_ij > 0.0)){
         if (u != INFTY)
            t_i = (u - x) / fabs(y_ij);
      }

      if (t_i <= t_max){
         if (fabs (y_ij) > fabs (yi_max)){
            yi_max  = y_ij;
            indx    = i;
            t_z     = t_i;
         }
      }
   }

   if (indx < 0){
      rs->ratio_stat = RATIO_FAILED;
   }
   else{
      /*
      if (tind != rs->lindex){
         HHHprintf ("tmax %e tval = %e yval = %e tind = %d\n", t_max, tval, yval, tind);
         HHHprintf ("h tval = %e yval = %e tind = %d\n",rs->tz, yi_max, rs->lindex);
      }
      */
      rs->lindex     = indx;
      rs->tz         = t_z;
      rs->pivotval   = yi_max;
      rs->ratio_stat = RATIO_BCHANGE;

      if (dir == VINCREASE)
         rs->lvstat = (yi_max > 0.0) ? STAT_LOWER : STAT_UPPER;
      else
         rs->lvstat = (yi_max > 0.0) ? STAT_UPPER : STAT_LOWER;

      if (rs->tz < 0.0){
         rs->tz      = fabs (t_max / 10.0);
         rs->boundch = 1;

         if (rs->lvstat == STAT_LOWER)
            rs->lbound = lp->xbz[rs->lindex] - rs->tz * fabs (yi_max);
         else
            rs->lbound = lp->xbz[rs->lindex] + rs->tz * fabs (yi_max);
      }
      if (dir == VDECREASE) rs->tz = -(rs->tz);
   }

   /*
   if (rs->lindex < 0){
      rs->ratio_stat = RATIO_FAILED;
   }
   else{
      rs->ratio_stat = RATIO_BCHANGE;
      if (dir == VINCREASE)
         rs->lvstat = (yi_max > 0.0) ? STAT_LOWER : STAT_UPPER;
      else
         rs->lvstat = (yi_max > 0.0) ? STAT_UPPER : STAT_LOWER;

      if (rs->tz < 0.0){
         rs->tz = 0.0;
         if (rs->lvstat == STAT_LOWER)
            rs->lbound = MIN (lp->xbz[rs->lindex],
                           lp->lz[lp->baz[rs->lindex]] - PARAM_BSHIFT * pftol);
         else
            rs->lbound = MAX (lp->xbz[rs->lindex],
                           lp->uz[lp->baz[rs->lindex]] + PARAM_BSHIFT * pftol);
         rs->boundch = 1;
      }
      rs->pivotval = yi_max;
      if (dir == VDECREASE) rs->tz = -(rs->tz);
   }
   */

CLEANUP:
   ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
}

#define GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
         x = -lp->dz[j]; \
         y = zAj; \
      } \
      else{ \
         x = lp->dz[j]; \
         y = -zAj; \
      } \
      if (lvstat == STAT_UPPER) \
         y = -y;


void ILLratio_dI_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs)
{
   int  j = 0, k;
   int  col;
   int  cbnd, indx;
   int  tctr  = 0;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   double  zAj;
   double  x, y, t_j = 0.0;
   double  theta, rcost;
   double  pftol = lp->tol->ip_tol;
   double  pivtol = lp->tol->pivot_tol;
   double  *t = lp->upd.t;
   double  delta = 0.0; /*HHH*/

   rs->eindex = -1;
   rs->tz = 0.0;
   rs->ratio_stat = RATIO_FAILED;

   for (k=0; k<lp->zA.nzcnt; k++){
      zAj = lp->zA.coef[k];
      if (zAj >= -pivtol && zAj <= pivtol)
         continue;

      t_j = INFTY;
      j   = lp->zA.indx[k];
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;

      GET_XY_DRATIOTEST;

      if (y < 0.0){
         if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO){
            t[tctr]  = x / y;
            ix[tctr] = 10*k + BBTOLOWER;
            tctr++;
         }
         else if (lp->vstat[col] == STAT_ZERO){
            if (lp->dfeas[j] < 0){
               t[tctr]  = x / y;
               ix[tctr] = 10*k + BBTOLOWER;
               tctr++;
            }
            if (lp->dfeas[j] <= 0){
               t[tctr]  = x / y;
               ix[tctr] = 10*k + BBTOUPPER;
               tctr++;
            }
         }
      }
      else{
         if (lp->dfeas[j] > 0){
            if (lp->vstat[col] == STAT_ZERO){
               t[tctr]  = x / y;
               ix[tctr] = 10*k + BATOUPPER;
               tctr++;
               t[tctr]  = x / y;
               ix[tctr] = 10*k + BATOLOWER;
               tctr++;
            }
         }
         else if (lp->dfeas[j] == 0){
            t[tctr]  = x / y;
            if (lp->vtype[col] == VBOUNDED)
               ix[tctr] = 10*k + BSKIP;
            else
               ix[tctr] = 10*k + BATOLOWER;
            tctr++;
         }
      }
   }

   if (tctr == 0){
      rs->ratio_stat = RATIO_FAILED;
      ILL_CLEANUP;
   }

   for (j=0; j<tctr; j++) perm[j] = j;
   ILLutil_double_perm_quicksort (perm, t, tctr);

   lp->upd.c_obj = 0.0; /*HHH*/

   rcost = (lvstat == STAT_LOWER) ? -lp->xbz[lindex] : lp->xbz[lindex];
   for (j=0; j<tctr; j++){
      cbnd = ix[perm[j]] % 10;
      if (cbnd == BSKIP) continue;

      t_j  = t[perm[j]];
      lp->upd.c_obj += (t_j-delta)*rcost; delta = t_j;/*HHH*/
      k    = ix[perm[j]] / 10;
      zAj  = lp->zA.coef[k];
      indx = lp->zA.indx[k];

      if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER || lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
         theta = -zAj;
      else
         theta = zAj;

      if (lvstat == STAT_UPPER)
         theta = -theta;

      switch (cbnd){
      case BATOLOWER:
      case BATOUPPER: rcost -= theta; break;
      case BBTOLOWER:
      case BBTOUPPER: rcost += theta; break;
      }
      if (rcost <= pftol){
         /* if (t_j < 0.0) printf ("dIhell\n"); */
         rs->eindex = indx;
         rs->tz = t_j;
         rs->pivotval = zAj;
         rs->ratio_stat = RATIO_BCHANGE;
         ILL_CLEANUP;
      }
   }

 CLEANUP:
   ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
   lp->upd.tctr = tctr;
   lp->upd.i    = j;
   lp->upd.tz   = fabs (t_j);
   lp->upd.piv  = rs->pivotval;
   if (rs->eindex != -1)
      lp->upd.fs = lp->dfeas[rs->eindex];
}

void ILLratio_dII_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs)
{
   int  j, k, indx;
   int  col, ecol;
   double  zAj;
   double  x, y;
   double  z_max;
   double  t_j, t_max, t_z;
   double  dftol  = lp->tol->dfeas_tol;
   double  pivtol = lp->tol->pivot_tol;

   rs->coeffch    = 0;
   rs->ecoeff    = 0.0;
   rs->eindex     = -1;
   rs->ratio_stat = RATIO_FAILED;
   lp->upd.tctr   = 0;
   lp->upd.dty    = 0.0;

   for (k=0, t_max = INFTY; k<lp->zA.nzcnt; k++){
      zAj = lp->zA.coef[k];
      if (zAj >= -pivtol && zAj <= pivtol)
         continue;

      t_j = INFTY;
      j   = lp->zA.indx[k];
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;

      GET_XY_DRATIOTEST;

      if (y > 0.0)
         t_j = (x + dftol) / y;
      else{
         if (lp->vstat[col] == STAT_ZERO)
            t_j = (x - dftol) / y;
      }
      if (t_j == INFTY) continue;

      if (t_max > t_j)
         t_max  = t_j;
   }

   if (t_max >= INFTY){
      rs->ratio_stat = RATIO_UNBOUNDED;
      ILL_CLEANUP; 
   }
   /* if (t_max < 0.0) printf ("dIIhell\n"); */

   indx  = -1;
   t_z   = 0.0;
   z_max = 0.0;

   for (k=0; k<lp->zA.nzcnt; k++){
      zAj = lp->zA.coef[k];
      if (zAj >= -pivtol && zAj <= pivtol)
         continue;

      t_j = INFTY;
      j   = lp->zA.indx[k];
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;

      GET_XY_DRATIOTEST;

      if (y > 0.0)
         t_j = x / y;
      else{
         if (lp->vstat[col] == STAT_ZERO)
            t_j = x / y;
      }

      if (t_j <= t_max){
         if (fabs (zAj) > fabs (z_max)){
            z_max = zAj;
            indx  = j;
            t_z   = t_j;
         }
      }
   }


   if (indx < 0){
      rs->ratio_stat = RATIO_FAILED;
   }
   else{
      rs->eindex     = indx;
      rs->tz         = t_z;
      rs->pivotval   = z_max;
      rs->ratio_stat = RATIO_BCHANGE;

      if (rs->tz < 0.0){
         rs->tz      = fabs (t_max / 20.0);
         rs->coeffch = 1;

         ecol = lp->nbaz[indx];
         if (lp->vstat[ecol] == STAT_LOWER)
            rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs (z_max);
         else if (lp->vstat[ecol] == STAT_UPPER)
            rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs (z_max);
         else{
            rs->tz = 0.0;
            rs->ecoeff = lp->cz[ecol] - lp->dz[indx];
         }
      }
   }

CLEANUP:
   ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
   lp->upd.piv = rs->pivotval;
}

void ILLratio_longdII_test (lpinfo *lp, int lindex, int lvstat, ratio_res *rs)
{
   int  j, k, indx = 0, tctr = 0;
   int  col, ecol;
   int  vs, bnd_exist = 0;
   int  *perm = lp->upd.perm;
   int  *ix   = lp->upd.ix;
   int  b_indx = -1;
   double  x, y, zAj = 0.0, z_max;
   double  xb, l, u;
   double  theta, rcost;
   double  t_j = 0.0, t_max, t_z;
   double  pftol  = lp->tol->pfeas_tol;
   double  dftol  = lp->tol->dfeas_tol;
   double  pivtol = lp->tol->pivot_tol;
   double  *t = lp->upd.t;
   double  delta = 0.0; /*HHH*/
   double  zb_val = 0.0;
   double  tb_val = NINFTY;

   rs->coeffch    = 0;
   rs->eindex     = -1;
   rs->ratio_stat = RATIO_FAILED;

   lp->upd.tctr   = 0;
   lp->upd.i      = 0;
   lp->upd.tz     = 0.0;
   lp->upd.piv    = 0.0;
   lp->upd.c_obj  = 0.0;
   lp->upd.dty    = 0.0;

   xb    = lp->xbz[lindex];
   col   = lp->baz[lindex];
   l     = lp->lz[col];
   u     = lp->uz[col];
   rcost = (lvstat == STAT_LOWER) ? l-xb : xb-u;

   for (k=0, t_max = INFTY; k<lp->zA.nzcnt; k++){
      zAj = lp->zA.coef[k];
      if (zAj >= -pivtol && zAj <= pivtol)
         continue;

      t_j = INFTY;
      j   = lp->zA.indx[k];
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
         continue;
      if (lp->vtype[col] == VBOUNDED){
         bnd_exist ++;
         continue;
      }

      GET_XY_DRATIOTEST;

      if (y > 0.0)
         t_j = (x + dftol) / y;
      else{
         if (lp->vstat[col] == STAT_ZERO)
            t_j = (x - dftol) / y;
      }
      if (t_j == INFTY) continue;

      if (t_max > t_j)
         t_max  = t_j;
   }
   if (t_max < 0.0){
      /*printf ("dIIhell, %.4f\n", t_max);*/
      rs->ratio_stat = RATIO_NEGATIVE;
      ILL_CLEANUP;
   }

   if (bnd_exist == 0 && t_max >= INFTY){
      rs->ratio_stat = RATIO_UNBOUNDED;
      /*
      printf ("x = %.8f, b = %.2f \n", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
      */
      ILL_CLEANUP; 
   }

   if (bnd_exist != 0){

      for (k=0; k<lp->zA.nzcnt; k++){
         zAj = lp->zA.coef[k];
         if (zAj >= -pivtol && zAj <= pivtol)
            continue;

         t_j = INFTY;
         j   = lp->zA.indx[k];
         col = lp->nbaz[j];

         if (lp->vtype[col] != VBOUNDED) continue;

         GET_XY_DRATIOTEST;

         if (y > 0.0){
            t_j = x / y;
            if (t_j <= t_max){
               t[tctr] = t_j;
               ix[tctr] = k;
               tctr ++;
            }
         }
      }
   }

   if (tctr != 0){
      for (j=0; j<tctr; j++) perm[j] = j;
      ILLutil_double_perm_quicksort (perm, t, tctr);

      for (j=0; j<tctr; j++){

         t_j  = t[perm[j]];
         lp->upd.c_obj += (t_j-delta)*rcost; delta = t_j;/*HHH*/
         k    = ix[perm[j]];
         zAj  = lp->zA.coef[k];
         indx = lp->zA.indx[k];
         col  = lp->nbaz[indx];
         l    = lp->lz[col];
         u    = lp->uz[col];
         vs   = lp->vstat[col];

         theta = (vs == STAT_UPPER) ? (l-u)*zAj : (u-l)*zAj;
         if (lvstat == STAT_LOWER)
            rcost += theta;
         else
            rcost -= theta;

         if (rcost <= pftol){
            rs->eindex = indx;
            rs->tz = t_j;
            rs->pivotval = zAj;
            rs->ratio_stat = RATIO_BCHANGE;

            if (rs->tz < 0.0){
               rs->tz = 0.0;
               rs->coeffch = 1;
               rs->ecoeff = lp->cz[col] - lp->dz[indx];
               lp->upd.c_obj += (rs->tz - delta)*rcost;
            }
            lp->upd.tctr = tctr;
            lp->upd.i    = j;
            lp->upd.tz   = rs->tz;
            ILL_CLEANUP;
         }
      }
      lp->upd.tctr = tctr;
      lp->upd.i    = tctr;
      lp->upd.tz   = t_j;
      zb_val = zAj;
      tb_val = t_j;
      b_indx = indx;
   }

   if (bnd_exist != 0 && t_max >= INFTY){
      rs->ratio_stat = RATIO_UNBOUNDED;
      /* printf ("rcost: %.8f\n", rcost); */
      ILL_CLEANUP; 
   }

   z_max = 0.0;
   indx = -1;
   t_z   = 0.0;

   for (k=0; k<lp->zA.nzcnt; k++){
      zAj = lp->zA.coef[k];
      if (zAj >= -pivtol && zAj <= pivtol)
         continue;

      t_j = INFTY;
      j   = lp->zA.indx[k];
      col = lp->nbaz[j];

      if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
          lp->vtype[col] == VBOUNDED)
         continue;

      GET_XY_DRATIOTEST;

      if (y > 0.0)
         t_j = x / y;
      else{
         if (lp->vstat[col] == STAT_ZERO)
            t_j = x / y;
      }

      if (t_j <= t_max){
         if (fabs (zAj) > fabs (z_max)){
            z_max = zAj;
            indx  = j;
            t_z   = t_j;
         }
      }
   }

   if (indx < 0){
      rs->ratio_stat = RATIO_FAILED;
      ILL_CLEANUP;
   }
   if ((tctr == 0) || (tb_val < 0.0) ||
       (tctr != 0 && t_z >= tb_val && fabs (z_max) >= fabs (zb_val))){
      lp->upd.c_obj += (t_z - delta) * rcost;
      delta          = t_z;
      rs->eindex     = indx;
      rs->tz         = t_z;
      rs->pivotval   = z_max;
      rs->ratio_stat = RATIO_BCHANGE;
   }
   else if ((tctr != 0) &&
       (fabs (zb_val) >= DBNDPIV_TOLER ||
        fabs (zb_val) >= DBNDPIV_RATIO * fabs (z_max))){
      rs->eindex     = b_indx;
      rs->tz         = tb_val;
      rs->pivotval   = zb_val;
      rs->ratio_stat = RATIO_BCHANGE;

      lp->upd.i     -= 1;
   }
   /* For now */
   else if (tctr != 0){
      rs->eindex     = b_indx;
      rs->tz         = tb_val;
      rs->pivotval   = zb_val;
      rs->ratio_stat = RATIO_BCHANGE;

      lp->upd.i     -= 1;
   }

   if (rs->tz < 0.0){
      /* if (tctr != 0) printf ("despite long step\n"); */

      rs->tz      = fabs (t_max / 20.0);
      rs->coeffch = 1;

      ecol = lp->nbaz[indx];
      if (lp->vstat[ecol] == STAT_LOWER)
         rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs (z_max);
      else if (lp->vstat[ecol] == STAT_UPPER)
         rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs (z_max);
      else{
         rs->tz = 0.0;
         rs->ecoeff = lp->cz[ecol] - lp->dz[indx];
      }
      lp->upd.c_obj += (rs->tz - delta) * rcost;
   }

CLEANUP:
   ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
   lp->upd.piv  = rs->pivotval;
}

void ILLratio_pivotin_test (lpinfo *lp, int *rlist, int rcnt, ratio_res *rs)
{
   int  i, k, col;
   double  x, l, u;
   double  y_ij, t_i;
   double  t_l, t_u;
   double  t_max, yi_max;
   double  pivtol = lp->tol->pivot_tol;

   if (rcnt <= 0 || rs == NULL) return;
   rs->boundch    = 0;
   rs->lindex     = -1;
   rs->tz         = 0.0;
   rs->ratio_stat = RATIO_FAILED;
   rs->lvstat     = -1;
   rs->pivotval   = 0.0;
   rs->lbound     = 0.0;

   for (i=0; i<rcnt; i++) lp->iwork[rlist[i]] = 1;

   for (k=0, t_max = INFTY; k<lp->yjz.nzcnt; k++){
      y_ij = lp->yjz.coef[k];
      if (y_ij >= -pivtol && y_ij <= pivtol)
         continue;

      i   = lp->yjz.indx[k];
      if (lp->iwork[lp->baz[i]] == 1) continue;
      x   = lp->xbz[i];
      col = lp->baz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];
      t_u = INFTY;
      t_l = NINFTY;

      if (l != NINFTY){
         t_l = (x - l) / y_ij;
         if (t_max > fabs (t_l)) t_max = fabs (t_l);
      }
      if (u != INFTY){
         t_u = (x - u) / y_ij;
         if (t_max > fabs (t_u)) t_max = fabs (t_u);
      }
   }

   if (t_max >= INFTY){
      rs->ratio_stat = RATIO_UNBOUNDED;
      ILL_CLEANUP; 
   }

   yi_max = 0.0;

   for (k=0; k<lp->yjz.nzcnt; k++){
      y_ij = lp->yjz.coef[k];
      if (y_ij >= -pivtol && y_ij <= pivtol)
         continue;

      i   = lp->yjz.indx[k];
      if (lp->iwork[lp->baz[i]] == 1) continue;
      x   = lp->xbz[i];
      col = lp->baz[i];
      l   = lp->lz[col];
      u   = lp->uz[col];

      t_u = INFTY;
      t_l = NINFTY;
      if (l != NINFTY) t_l = (x - l) / y_ij;
      if (u != INFTY)  t_u = (x - u) / y_ij;
      t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u;

      if (fabs (t_i) <= t_max + t_max*(1.0e-2)){
         if (fabs (y_ij) > fabs (yi_max)){
            yi_max     = y_ij;
            rs->lindex = i;
            rs->tz     = t_i;
            rs->lvstat = (t_i == t_l) ? STAT_LOWER : STAT_UPPER;
         }
      }
   }

   if (rs->lindex < 0){
      rs->ratio_stat = RATIO_FAILED;
   }
   else{
      rs->ratio_stat = RATIO_BCHANGE;
      rs->pivotval   = yi_max;
   }
 CLEANUP:
   for (i=0; i<rcnt; i++) lp->iwork[rlist[i]] = 0;
   return;
}

