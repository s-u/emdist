/* RCS_INFO = "$RCSfile: lpdata.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
/****************************************************************************/
/*                                                                          */
/*               Routines for Manipulating and Writing LPdata               */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLlpdata_buildrows (ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
/*            int **rowind, double **rowval, int include_logicals)          */
/*      - include_logicals:  if nonzero, then logical variables will be     */
/*          included in the row data                                        */
/*                                                                          */
/*                                                                          */
/*    All _init routines initialize fields of allocated structure to        */
/*    appropiate default values                                             */
/*    The _free routines free structures contained in pareameter structure  */ 
/*    but not the parameter itself.                                         */
/*    The _alloc routines check whether given parameter is NULL; they either*/
/*    print an error message or fill structure with default values or the   */
/*    given paremeter values.                                               */ 
/*                                                                          */
/*    void ILLlpdata_init (ILLlpdata *lp)                                   */
/*    void ILLlpdata_free (ILLlpdata *lp)                                   */
/*                                                                          */
/*    void ILLlp_basis_init (ILLlp_basis *B)                                */
/*    void ILLlp_basis_free (ILLlp_basis *B)                                */
/*    int ILLlp_basis_alloc (ILLlp_basis *B, int nstruct, int nrows)        */
/*                                                                          */
/*    void ILLlp_cache_init (ILLlp_cache *C)                                */
/*    void ILLlp_cache_free (ILLlp_cache *C)                                */
/*    int ILLlp_cache_alloc (ILLlp_cache *C, int nstruct, int nrows)        */
/*                                                                          */
/*    void ILLlp_sinfo_init (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_sinfo_free (ILLlp_sinfo *sinfo)                            */
/*                                                                          */
/*    int ILLlp_rows_init(ILLlp_rows *lprows, ILLlpdata *lp,                */
/*                                           int include_logicals)          */
/*                                                                          */
/****************************************************************************/

#include "iqsutil.h"
#include "lpdata.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lp.h"
#include "mps.h"
#include "rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0; 

QSdata* ILLread(qsline_reader *file, const char* fname, int isMps) 
{
    int rval = 0;
    QSdata *p = (QSdata *) NULL;
    ILLlpdata *lp;
    rawlpdata rawlp;

    ILL_FAILfalse(file != NULL, NULL);
    ILL_FAILfalse(fname != NULL, NULL);

    p = QScreate_prob (fname, QS_MIN); 
    ILL_CHECKnull(p, NULL);
    ILL_IFFREE(p->qslp->probname, char);
    lp = p->qslp;

    ILLinit_rawlpdata(&rawlp, file->error_collector);
    ILLlpdata_init(lp); 

    if (isMps != 0) {
        rval = ILLread_mps(file, fname, &rawlp);
    } else {
        rval = ILLread_lp(file, fname, &rawlp); 
    }
    ILL_CLEANUP_IF(rval); 

    rval = ILLrawlpdata_to_lpdata (&rawlp, lp);
    ILL_CLEANUP_IF(rval); 

CLEANUP:
    ILLfree_rawlpdata(&rawlp); 
    if (rval != 0) {
        QSfree_prob (p);
        p = (QSdata *) NULL;
    }
    return p;
}

void ILLlpdata_init (ILLlpdata *lp) 
{
    if (lp) {
        lp->nrows      = 0;
        lp->ncols      = 0;
        lp->nstruct    = 0;
        lp->nzcount    = 0;
        lp->rowsize    = 0;
        lp->colsize    = 0;
        lp->structsize = 0;
        lp->objsense   = ILL_MIN;
        lp->sense      = (char *) NULL;
        lp->obj        = (double *) NULL;
        lp->rhs        = (double *) NULL;
        lp->rangeval   = (double *) NULL;
        lp->lower      = (double *) NULL;
        lp->upper      = (double *) NULL;

        ILLmatrix_init (&lp->A);
        ILLmatrix_init (&lp->sos);
        lp->rA         = (ILLlp_rows *) NULL;
        lp->is_sos_mem = NULL; 
        lp->refrowname = NULL; 
        lp->refind     = -1; 

        lp->colnames   = (char **) NULL;
        ILLsymboltab_init(&lp->coltab); 
        lp->rownames   = (char **) NULL;
        ILLsymboltab_init(&lp->rowtab); 
        lp->objname    = (char *) NULL;

        lp->probname   = (char *) NULL;
        lp->intmarker  = (char *) NULL;
        lp->structmap  = (int *) NULL;
        lp->rowmap     = (int *) NULL;
        lp->basis      = (ILLlp_basis *) NULL;
        /*lp->presolve   = (ILLlp_predata *) NULL; */
        lp->sinfo      = (ILLlp_sinfo *) NULL;

		ILLstring_reporter_init(&lp->reporter, ILL_fprintf, stdout); 
    }
}

void ILLlpdata_free (ILLlpdata *lp)
{
    int i;

    if (lp) {
        ILL_IFFREE (lp->sense, char);
        ILL_IFFREE (lp->obj, double);
        ILL_IFFREE (lp->rhs, double);
        ILL_IFFREE (lp->rangeval, double);
        ILL_IFFREE (lp->lower, double);
        ILL_IFFREE (lp->upper, double);
        ILLmatrix_free (&lp->A);
        if (lp->rA) {
            ILLlp_rows_clear(lp->rA); 
            ILL_IFFREE (lp->rA, ILLlp_rows);
        }
        ILL_IFFREE (lp->is_sos_mem, int);
        ILL_IFFREE (lp->refrowname, char);
        ILLmatrix_free (&lp->sos);
        if (lp->colnames) {
            for (i = 0; i < lp->nstruct; i++) {
                ILL_IFFREE (lp->colnames[i], char);
            }
            ILL_IFFREE (lp->colnames, char *);
        }
        ILLsymboltab_free(&lp->coltab); 
        if (lp->rownames) {
            for (i = 0; i < lp->nrows; i++) {
                ILL_IFFREE (lp->rownames[i], char);
            }
            ILL_IFFREE (lp->rownames, char *);
        }
        ILLsymboltab_free(&lp->rowtab); 
        ILL_IFFREE (lp->objname, char);
        ILL_IFFREE (lp->probname, char);
        ILL_IFFREE (lp->intmarker, char);
        ILL_IFFREE (lp->structmap, int);
        ILL_IFFREE (lp->rowmap, int);
        if (lp->sinfo) {
            ILLlp_sinfo_free (lp->sinfo);
            ILL_IFFREE (lp->sinfo, ILLlp_sinfo);
        }
        ILLlpdata_init (lp);
    }
}

void ILLlp_basis_init (ILLlp_basis *B)
{
    if (B) {
        B->cstat = (char *) NULL;
        B->rstat = (char *) NULL;
        B->rownorms = (double *) NULL;
        B->colnorms = (double *) NULL;
        B->nstruct = 0;
        B->nrows   = 0;
    }
}

void ILLlp_basis_free (ILLlp_basis *B)
{
    if (B) {
        ILL_IFFREE (B->cstat, char);
        ILL_IFFREE (B->rstat, char);
        ILL_IFFREE (B->rownorms, double);
        ILL_IFFREE (B->colnorms, double);
        B->nstruct = 0;
        B->nrows   = 0;
    }
}

int ILLlp_basis_alloc (ILLlp_basis *B, int nstruct, int nrows)
{
    int rval = 0;

    ILL_FAILtrue(B == NULL, "ILLlp_basis_alloc called without a basis");

    B->nstruct = nstruct;
    B->nrows   = nrows;

    if (nstruct > 0) {
        ILL_SAFE_MALLOC(B->cstat, nstruct, char);
    }

    if (nrows > 0) {
        ILL_SAFE_MALLOC(B->rstat, nrows, char);
    }

CLEANUP:

    if (rval) {
        ILLlp_basis_free (B);
    }

    ILL_RETURN(rval, "ILLlp_basis_alloc");
}

void ILLlp_cache_init (ILLlp_cache *C)
{
    if (C) {
        C->x     = (double *) NULL;
        C->rc    = (double *) NULL;
        C->pi    = (double *) NULL;
        C->slack = (double *) NULL;
        C->val     = 0.0;
        C->nstruct = 0;
        C->nrows   = 0;
        C->status  = 0;
    }
}

void ILLlp_cache_free (ILLlp_cache *C)
{
    if (C) {
        ILL_IFFREE (C->x, double);
        ILL_IFFREE (C->rc, double);
        ILL_IFFREE (C->pi, double);
        ILL_IFFREE (C->slack, double);
        C->nstruct = 0;
        C->nrows   = 0;
        C->status  = 0;
    }
}

int ILLlp_cache_alloc (ILLlp_cache *C, int nstruct, int nrows)
{
    int rval = 0;

    ILL_FAILtrue(C == NULL, "ILLlp_cache_alloc called without a cache");

    C->nstruct = nstruct;
    C->nrows   = nrows;

    if (nstruct > 0) {
        ILL_SAFE_MALLOC(C->x, nstruct, double); 
        ILL_SAFE_MALLOC(C->rc, nstruct, double);
    }

    if (nrows > 0) {
        ILL_SAFE_MALLOC(C->pi, nrows, double);
        ILL_SAFE_MALLOC(C->slack, nrows, double);
    }

CLEANUP:

    if (rval) {
        ILLlp_cache_free (C);
    }

    ILL_RETURN(rval, "ILLlp_cache_alloc");
}


int ILLlp_rows_init(ILLlp_rows *lprows, ILLlpdata *lp, int include_logicals) 
{
    int rval = 0;
    int i, k, st;
    int *beg, *cnt, *ind;
    double *val;
    ILLmatrix *A;
    char *hit = (char *) NULL;
    int *inv_structmap = (int *) NULL;

    /* If logicals are not included, then the columns are ordered as in */
    /* lp->structmap.  Otherwise, the columns are ordered as in the     */
    /* matrix structure.                                                */

    if (lprows != NULL) { 
        lprows->rowbeg = (int *) NULL;
        lprows->rowcnt = (int *) NULL;
        lprows->rowind = (int *) NULL;
        lprows->rowval = (double *) NULL;
    }

    ILL_FAILfalse((lp != NULL) && (lprows != NULL), 
                  "called with a NULL pointer");

    A = &lp->A;

    if (lp->nrows > 0)  {
	if (include_logicals == 0) {
	    ILL_FAILtrue(lp->rowmap == NULL, "Programming error.");
	    ILL_SAFE_MALLOC (hit, lp->ncols, char);

	    for (i = 0; i < lp->ncols; i++) {
		hit[i] = 0;
	    }
	    for (i = 0; i < lp->nrows; i++) {
		hit[lp->rowmap[i]] = 1;
	    }

	    ILL_SAFE_MALLOC (inv_structmap, lp->ncols, int);

	    for (i = 0; i < lp->nstruct; i++) {
		inv_structmap[lp->structmap[i]] = i;
	    }
	}

	ILL_SAFE_MALLOC(lprows->rowbeg, lp->nrows, int);
	ILL_SAFE_MALLOC(lprows->rowcnt, lp->nrows, int);

	if ( ((include_logicals != 0) && lp->nzcount > 0) ||
	     ((include_logicals == 0) && lp->nzcount > lp->nrows)) {
	    if (include_logicals != 0) {
		ILL_SAFE_MALLOC(lprows->rowind, lp->nzcount, int);
		ILL_SAFE_MALLOC(lprows->rowval, lp->nzcount, double);
	    } else {
		ILL_SAFE_MALLOC(lprows->rowind, lp->nzcount - lp->nrows, int);
		ILL_SAFE_MALLOC(lprows->rowval, lp->nzcount - lp->nrows, double);
	    }
	}

	beg = lprows->rowbeg;
	cnt = lprows->rowcnt;
	ind = lprows->rowind;
	val = lprows->rowval;

	for (i = 0; i < lp->nrows; i++) {
	    cnt[i] = 0;
	}

	for (i = 0; i < lp->ncols; i++) {
	    if ((include_logicals != 0) || hit[i] == 0) {
		k  =     A->matbeg[i];
		st = k + A->matcnt[i];
		for (; k < st; k++) {
		    cnt[A->matind[k]]++;
		}
	    }
	}

	for (i = 0, k = 0; i < lp->nrows; i++) {
	    beg[i] = k;
	    k += cnt[i];
	}

	for (i = 0; i < lp->ncols; i++) {
	    if ((include_logicals != 0) || hit[i] == 0) {
		k  =     A->matbeg[i];
		st = k + A->matcnt[i];
		for (; k < st; k++) {
		    if (include_logicals != 0) {
			ind[beg[A->matind[k]]] = i;
		    } else {
			ind[beg[A->matind[k]]] = inv_structmap[i];
		    }
		    val[beg[A->matind[k]]] = A->matval[k];
		    beg[A->matind[k]]++;
		}
	    }
	}

	for (i = 0, k = 0; i < lp->nrows; i++) {
	    beg[i] = k;
	    k += cnt[i];
	}
    }
CLEANUP:

    if (rval) {
        ILLlp_rows_clear(lprows); 
    }
    ILL_IFFREE (hit, char);
    ILL_IFFREE (inv_structmap, int);

    ILL_RETURN(rval, "ILLlp_rows_init");
}

void ILLlp_rows_clear(ILLlp_rows *lprows) 
{
    if (lprows != NULL) { 
	ILL_IFFREE(lprows->rowbeg, int); 
	ILL_IFFREE(lprows->rowcnt, int); 
	ILL_IFFREE(lprows->rowind, int); 
	ILL_IFFREE(lprows->rowval, double); 
    }
}

static int wr_line(ILLlpdata *lp, const char *format, va_list argptr)
{
	char buffer[256];
	int rval = 0; 
	rval = vsprintf(buffer, format, argptr);
	if (rval > 0) { 
		rval = ILLstring_report(buffer, &lp->reporter); 
	}
	return rval;
}

int ILLprint_report(ILLlpdata *lp, const char *format, ...) 
{
   va_list marker;
   int rval = 0;

   va_start( marker, format);     /* ANSI style */
   rval = wr_line(lp, format, marker); 
   va_end( marker );              /* Reset variable arguments.      */
   return rval; 
}


