/* RCSINFO $Id: lpdefs.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __QS_LPDEFS_H
#define __QS_LPDEFS_H

#include "qsopt.h"
#include "lpdata.h"
#include "factor.h"

/* infinity and negative infinity */
#define INFTY  (1e30)
#define NINFTY (-1e30)

/* storage type */
#define DENSE          0
#define SPARSE         1

/* type of vector */
#define ROW_SOLVE      1
#define COLUMN_SOLVE   2

/* direction of change in non-basic var */
#define VINCREASE      1
#define VDECREASE      2

/* status of variables */
#define STAT_BASIC     1
#define STAT_UPPER     2
#define STAT_LOWER     3
#define STAT_ZERO      4

#define BOUND_LOWER    1
#define BOUND_UPPER    2

/* type of variables */
#define VARTIFICIAL    1
#define VFIXED         2
#define VFREE          3
#define VUPPER         4
#define VLOWER         5
#define VBOUNDED       6

/* algo */
#define PRIMAL_SIMPLEX 1
#define DUAL_SIMPLEX   2
#define PRIMAL_OR_DUAL 3

/* phase */
#define PRIMAL_PHASEI  1
#define PRIMAL_PHASEII 2
#define DUAL_PHASEI    3
#define DUAL_PHASEII   4

/* number of phases */
#define PHASEI         1
#define PHASEII        2

/* type of pricing (all vars or some) */
#define COMPLETE_PRICING   1
#define PARTIAL_PRICING    2
#define MULTI_PART_PRICING 3

/* default pricing */

#define QS_DEFAULT_PRICE_PI  QS_PRICE_PSTEEP
#define QS_DEFAULT_PRICE_PII QS_PRICE_PSTEEP
#define QS_DEFAULT_PRICE_DI  QS_PRICE_DSTEEP
#define QS_DEFAULT_PRICE_DII QS_PRICE_DSTEEP

/* lp sol status */
#define ILL_LP_SOLVED            1
#define ILL_LP_UNSOLVED          2
#define ILL_MAX_ITER             3
#define ILL_MAX_TIME             4
#define ILL_BND_REACHED          5
#define ILL_PPHASEI_ERROR        6
#define ILL_PPHASEII_ERROR       7
#define ILL_DPHASEI_ERROR        8
#define ILL_DPHASEII_ERROR       9
#define ILL_LP_ABORTED          10

/* basis status */
#define OPTIMAL                  1
#define NONOPTIMAL               2
#define PRIMAL_FEASIBLE          3
#define PRIMAL_INFEASIBLE        4
#define PRIMAL_UNBOUNDED         5
#define DUAL_FEASIBLE            7
#define DUAL_INFEASIBLE          8
#define DUAL_UNBOUNDED           9

/* type of ratio test */
#define RATIOTEST_NORMAL 1
#define RATIOTEST_HARRIS 2

/* control parameters */
#define PARAM_PRATIOTESTS        10
#define PARAM_DRATIOTESTS        20
#define PARAM_PRIMAL_REFACTORGAP 50
#define PARAM_PRIMAL_RESOLVEGAP  25
#define PARAM_DUAL_REFACTORGAP   100
#define PARAM_DUAL_RESOLVEGAP    25
#define PARAM_MAX_NOSOLVE        500
#define PARAM_MAX_NOPROG         300
#define PARAM_NOPROG_FACTOR      15

/* numerical parameters */       
#define PARAM_BSHIFT             10   
#define PARAM_CSHIFT             10   
#define PARAM_IBASIS_RPIVOT      0.98 
#define PARAM_IBASIS_RTRIANG     0.01 
#define PARAM_MIN_DNORM          1e-24

/* general constants */
#define PARAM_HEAP_UTRIGGER      10
#define PARAM_HEAP_RATIO         4.0

/* errors */
#define E_GENERAL_ERROR          1
#define E_INV_LINSOLVE_OPTION    2
#define E_NO_MEMORY              3
#define E_INVALID_OPTION         4
#define E_NULL_ARGUMENT          5
#define E_SIMPLEX_ERROR          6
#define E_BASIS_SINGULAR         7

/* tolerances */
#define PFEAS_TOLER              1e-6
#define DFEAS_TOLER              1e-6
#define PIVOT_TOLER              1e-10
#define SZERO_TOLER              1e-15
#define PIVZ_TOLER               1e-12
#define OBJBND_TOLER             1e-2
#define DBNDPIV_TOLER            1e-3
#define DBNDPIV_RATIO            1e-2
#define PNSTEP_TOLER             1e-9
#define DNSTEP_TOLER             1e-9
#define ALTPIV_TOLER             1e-8

/* structure for statistics */
typedef struct{
   int ynz_cnt;           /* nz in entering columns */
   int num_y;
   double y_ravg;         /* weighted avg. of current & prior y */
   int znz_cnt;           /* nz in ith row of B^{-1}, ie z_i */
   int num_z;
   double z_ravg;         /* weighted avg. of current & prior z */
   int zanz_cnt;          /* nz in z^TA */
   int num_za;
   double za_ravg;        /* weighted avg. of current & prior za */
   int pnorm_cnt;         /* nz in columns for primal norms */
   int dnorm_cnt;         /* nz in rows for dual norms */
   int pinz_cnt;          /* nz in phase II pi (solve) */
   int num_pi;            /* # of pi solves */
   int pi1nz_cnt;         /* nz in phase I pi (solve) */
   int num_pi1;           /* # of phase I pi solves */
   int upnz_cnt;          /* nz in ftran update vector */
   int num_up;            /* # of ftran_updates */
   int pupv_cnt;          /* nz in primal steep updates */
   int dupv_cnt;          /* nz in dual steep updates */

   int start_slacks;      /* # slacks in beginning */
   int final_slacks;      /* # slacks in the end */
   int start_art;         /* # arts in beginning */
   int final_art;         /* # arts in the end */

   int pI_iter;           /* primal phase I iterations */
   int pII_iter;
   int dI_iter;           /* dual phase I iterations */
   int dII_iter;
   int tot_iter;

   int pivpI[10];         /* sizes of pivots */
   int pivpII[10];
   int pivdI[10];
   int pivdII[10];
} count_struct;

/* structure for tolerances */
typedef struct{
   double pfeas_tol;
   double dfeas_tol;
   double pivot_tol;
   double szero_tol;
   double ip_tol;         /* inner primal & dual feas toler */
   double id_tol;
} tol_struct;

/* bound information */
typedef struct bndinfo{
   double pbound;
   double cbound;
   int btype;
   int varnum;
   struct bndinfo *next;
} bndinfo;

/* bound information */
typedef struct coefinfo{
   double pcoef;
   double ccoef;
   int varnum;
   struct coefinfo *next;
} coefinfo;

/* feasibility info */
typedef struct feas_info{
   int pstatus;
   int dstatus;
   double totinfeas;
} feas_info;

typedef struct lp_status_info{
   char optimal;
   char primal_feasible;
   char primal_infeasible;
   char primal_unbounded;
   char dual_feasible;
   char dual_infeasible;
   char dual_unbounded;
} lp_status_info;

typedef struct pI_uinfo{
   int tctr;
   int i;
   int *perm;
   int *ix;
   int fs;
   double piv;
   double *t;
   double dty;
   double c_obj;
   double tz;
} pI_uinfo;

extern void ILLlp_status_info_init (lp_status_info *ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct lpinfo {

   double objval;         /* obj info */
   lp_status_info probstat;   /* final status */
   lp_status_info basisstat;  /* final status */
   
   double pobjval;        /* intermediate status info */
   double dobjval;
   double pinfeas;
   double dinfeas;
   double objbound;

   int    nrows;          /* input info follows; given in col format*/
   int    ncols;
   int    *matcnt;
   int    *matbeg;
   int    *matind;
   double *matval;
   int     matfree;
   int     matsize;
   double *bz;
   double *lz;
   double *uz;
   double *cz;
   int    localrows;      /* set to 1 if these are created locally */
   int    *rowcnt;        /* row info follows, copy of col info */
   int    *rowbeg;
   int    *rowind;
   double *rowval;  

   double *xbz;           /* output info x, pi, reduced cost*/
   double *piz;
   double *dz;
   double *pIxbz;         /* output info (phase I) x, pi, reduced cost*/
   double *pIpiz;
   double *pIdz;

   int final_phase;       /* final phase, inf & unboundedness info */ 
   int infub_ix;

   int basisid;           /* basis and variable info follows */
   int nnbasic;
   int *baz;
   int *nbaz;
   int *vstat;
   int *vindex;
   int fbasisid;
   factor_work *f;
   int *vtype;            /* internal var info */

   svector zz;            /* local ILLfactor_update vectors z, yj, za */
   svector yjz;
   svector zA;
   svector work;          /* local work vector */
   svector srhs;          /* local vectors for lin. eq. solves */
   svector ssoln;
   int *iwork;            /* local work vector */
   pI_uinfo upd;          /* phase I update info */
   int *bfeas;            /* primal and dual infeasibility info */
   int *dfeas;

   tol_struct *tol;       /* tolerances */
   count_struct *cnts;    /* counts */
   int nbchange;          /* # bound shifts */
   int ncchange;          /* # obj coef shifts */
   bndinfo *bchanges;     /* list of bound shifts */
   coefinfo *cchanges;    /* list of coef shifts */
   int pIratio;           /* ratio tests */
   int pIIratio;
   int dIratio;
   int dIIratio;

   int maxiter;
   int iterskip;
   double maxtime;
   double starttime;
   struct ILLlpdata *O;
   ILLrandstate rstate;

} lpinfo;

/* pricing structures */
typedef struct{
   int ninit;
   double *norms;
   int *refframe;
} p_devex_info;

typedef struct{
   double *norms;
} p_steep_info;

typedef struct{
   int  k;
   int  cgroup;
   int  ngroups;
   int  *gstart;
   int  *gshift;
   int  *gsize;
   int  bsize;
   int  *bucket;
   int  *perm;
   double  *infeas;
} mpart_info;

typedef struct{
   int ninit;
   double *norms;
   int *refframe;
} d_devex_info;

typedef struct{
   double *norms;
} d_steep_info;

/* pricing information */
typedef struct price_info {
   int p_strategy;
   int d_strategy;
   int pI_price;
   int pII_price;
   int dI_price;
   int dII_price;
   int cur_price;
   double *p_scaleinf;
   double *d_scaleinf;
   p_devex_info   pdinfo;
   p_steep_info   psinfo;
   mpart_info     pmpinfo;
   d_devex_info   ddinfo;
   d_steep_info   dsinfo;
   mpart_info     dmpinfo;
   heap   h;
   double htrigger;
   int    hineff;
} price_info;

#endif /* __QS_LPDEFS_H */
