/*  $RCSfile: qstruct.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __QS_QSTRUCT_H
#define __QS_QSTRUCT_H

typedef struct qsdata {
    struct ILLlpdata   *qslp;
    struct lpinfo      *lp;  
    struct price_info  *pricing;
    struct ILLlp_basis *basis; 
    struct ILLlp_cache *cache; 
    char               *name;
    int                qstatus;    /* QS_LP_UNSOLVED or an opt status  */
    int                factorok;   /* set to 0 if factorization is old */
    int                simplex_display;  /* 0 off, 1 on */
    int                simplex_scaling;  /* 0 off, 1 on */
} QSdata; 

typedef struct qsbasis {
    int    nstruct;
    int    nrows;
    char  *cstat;
    char  *rstat;
} QSbasis;

#endif /* __QS_QSTRUCT_H */

