/* RCSINFO $Id: factor.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __QS_FACTOR_H_
#define __QS_FACTOR_H_

#define FALSE 0
#define TRUE  1

typedef char QSbool;

typedef struct uc_info {
    int cbeg;
    int nzcnt;
    int next;
    int prev;
    int delay;
} uc_info;

typedef struct ur_info {
    double max;
    int rbeg;
    int nzcnt;
    int pivcnt;
    int next;
    int prev;
    int delay;
} ur_info;

typedef struct lc_info {
    int cbeg;
    int nzcnt;
    int c;
    int crank;
    int delay;
} lc_info;

typedef struct lr_info {
    int rbeg;
    int nzcnt;
    int r;
    int rrank;
    int delay;
} lr_info;

typedef struct er_info {
    int rbeg;
    int nzcnt;
    int r;
} er_info;

#define QS_FACTOR_MAX_K        1
#define QS_FACTOR_P            2
#define QS_FACTOR_ETAMAX       3
#define QS_FACTOR_FZERO_TOL    4
#define QS_FACTOR_SZERO_TOL    5
#define QS_FACTOR_PARTIAL_TOL  6
#define QS_FACTOR_UR_SPACE_MUL 7
#define QS_FACTOR_UC_SPACE_MUL 8
#define QS_FACTOR_LC_SPACE_MUL 9
#define QS_FACTOR_LR_SPACE_MUL 10
#define QS_FACTOR_ER_SPACE_MUL 11
#define QS_FACTOR_GROW_MUL     12
#define QS_FACTOR_MAXMULT      13
#define QS_FACTOR_MINMULT      14
#define QS_FACTOR_UPDMAXMULT   15
#define QS_FACTOR_DENSE_FRACT  16
#define QS_FACTOR_DENSE_MIN    17

typedef struct factor_work {
    int       max_k;
    double    fzero_tol;
    double    szero_tol;
    double    partial_tol;
    double    ur_space_mul;
    double    uc_space_mul;
    double    lc_space_mul;
    double    lr_space_mul;
    double    er_space_mul;
    double    grow_mul;
    int       p;
    int       etamax;
    double    minmult;
    double    maxmult;
    double    updmaxmult;
    double    dense_fract;
    int       dense_min;

    double    maxelem_orig;
    int       nzcnt_orig;
    double    maxelem_factor;
    int       nzcnt_factor;
    double    maxelem_cur;
    int       nzcnt_cur;

    double    partial_cur;

    int       dim;
    int       stage;
    int       nstages;
    int       etacnt;
    double   *work_coef;
    int      *work_indx;
    uc_info  *uc_inf;
    ur_info  *ur_inf;
    lc_info  *lc_inf;
    lr_info  *lr_inf;
    er_info  *er_inf;
    int      *ucindx;          /* row index for column data */
    int      *ucrind;          /* index of column in row data */
    double   *uccoef;          /* coefficient for column data */
    int      *urindx;          /* col index for row data */
    int      *urcind;          /* index of row in column data */
    double   *urcoef;          /* coefficient for row data */
    int      *lcindx;          /* row index for L data */
    double   *lccoef;          /* coefficient for L row data */
    int      *lrindx;          /* col index for L data */
    double   *lrcoef;          /* coefficient for L col data */
    int      *erindx;          /* col index for eta data */
    double   *ercoef;          /* coefficient for eta data */
    int      *rperm;
    int      *rrank;
    int      *cperm;
    int      *crank;
    svector   xtmp;
    int       ur_freebeg;
    int       ur_space;
    int       uc_freebeg;
    int       uc_space;
    int       lc_freebeg;
    int       lc_space;
    int       lr_freebeg;
    int       lr_space;
    int       er_freebeg;
    int       er_space;

    int      *p_nsing;
    int     **p_singr;
    int     **p_singc;
    
    double   *dmat;
    int       drows;
    int       dcols;
    int       dense_base;
} factor_work;

void
    ILLfactor_init_factor_work (factor_work *f),
    ILLfactor_free_factor_work (factor_work *f),
    ILLfactor_ftran (factor_work *f, svector *a, svector *x),
    ILLfactor_ftran_update (factor_work *f, svector *a, svector *upd, svector *x),
    ILLfactor_btran (factor_work *f, svector *a, svector *x);

int
    ILLfactor_create_factor_work (factor_work *f, int dim),
    ILLfactor_set_factor_iparam (factor_work *f, int param, int val),
    ILLfactor_set_factor_dparam (factor_work *f, int param, double val),
    ILLfactor (factor_work *f, int *basis, int *cbeg, int *clen, int *cindx,
        double *ccoef, int *p_nsing, int **p_singr, int **p_singc),
    ILLfactor_update (factor_work *f, svector *a, int col, int *p_refact);

#define E_CHECK_FAILED 6
#define E_NO_PIVOT 7
#define E_FACTOR_BLOWUP 8
#define E_UPDATE_NOSPACE 9
#define E_UPDATE_SINGULAR_ROW 10
#define E_UPDATE_SINGULAR_COL 11
#define E_SING_NO_DATA 12
#define E_SINGULAR_INTERNAL 13

#endif /* __QS_FACTOR_H_ */

