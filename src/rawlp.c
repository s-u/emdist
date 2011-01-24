/* RCS_INFO = "$RCSfile: rawlp.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
/****************************************************************************/
/* DataStructure and routines to deal with raw lp information as read       */
/* from mps or lp files.                                                    */ 
/****************************************************************************/

#include "iqsutil.h"
#include "rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
static int TRACE = 0; 

ILL_PTRWORLD_ROUTINES(colptr, colptralloc, colptr_bulkalloc, colptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE(colptr, colptr_listfree, colptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (colptr, colptr_check_leaks, this, int)

const int ILL_SOS_TYPE1 = 1;
const int ILL_SOS_TYPE2 = 2;

static void ILLprt_double(FILE *f, double d) 
{
    if (ILL_MAXDOUBLE == d) {
	fprintf(f, "MAX_DOUBLE");  
    } else { 
        if (-ILL_MAXDOUBLE == d) {
	    fprintf(f, "-MAX_DOUBLE");  
        } else {
	    fprintf(f, "%f", d);  
        }
    }
}

static int ILLraw_check_bounds(rawlpdata *lp);
 
void ILLinit_rawlpdata(rawlpdata *lp, qserror_collector *collector)
{
    if (lp) {
        lp->name  = (char *) NULL;
        lp->ncols = 0;
        lp->nrows = 0;
        lp->cols        = (colptr **) NULL;
        ILLsymboltab_init(&lp->coltab); 
        ILLsymboltab_init(&lp->rowtab); 
        lp->rowsense    = (char *) NULL;
        lp->rhsname     = (char *) NULL;
        lp->rhs         = (double *) NULL;
        lp->rhsind      = (char *) NULL;
        lp->rangesname  = (char *) NULL;
        lp->rangesind   = (char *) NULL;
        lp->ranges      = (colptr *) NULL;
        lp->boundsname  = (char *) NULL;
	lp->lbind       = (char *) NULL;
	lp->ubind       = (char *) NULL;
        lp->lower       = (double *) NULL;
        lp->upper       = (double *) NULL;
        lp->intmarker   = (char *) NULL;
        lp->colsize     = 0;
        lp->sensesize   = 0;
        lp->intsize     = 0;
        lp->objindex    = -1;
        lp->rhssize     = 0;
        lp->objsense    = ILL_MIN;

        lp->refrow = NULL; 
        lp->refrowind = -1; /* undefined */
	lp->is_sos_size = 0; 
	lp->is_sos_member = NULL; 
	lp->nsos_member = 0; 
	lp->sos_weight_size= 0; 
	lp->sos_weight = NULL; 
	lp->sos_col_size= 0; 
	lp->sos_col = NULL; 
	lp->nsos = 0;
	lp->sos_setsize = 0; 
	lp->sos_set = NULL; 
        ILLptrworld_init (&lp->ptrworld);
	lp->error_collector = collector; 
    }
}

void ILLraw_clear_matrix(rawlpdata *lp)
{
    int i;
    if ((lp != NULL) && (lp->cols != NULL)) {
	for (i = 0; i < lp->ncols; i++) {
	    colptr_listfree (&lp->ptrworld, lp->cols[i]);
            lp->cols[i] = NULL; 
	}
    }
}

void ILLfree_rawlpdata (rawlpdata *lp)
{
    int total, onlist;

    if (lp) {
        ILL_IFFREE (lp->name, char);
        ILLsymboltab_free(&lp->rowtab); 
        ILLsymboltab_free(&lp->coltab); 
        ILL_IFFREE (lp->rowsense, char);
        ILLraw_clear_matrix(lp); 
	ILL_IFFREE (lp->cols, colptr *);
        colptr_listfree (&lp->ptrworld, lp->ranges);
        if (colptr_check_leaks (&lp->ptrworld, &total, &onlist)) {
            fprintf (stderr, "WARNING: %d outstanding colptrs\n",
                     total - onlist);
        }
        ILLptrworld_delete (&lp->ptrworld);
        ILL_IFFREE (lp->rhsname, char);
        ILL_IFFREE (lp->rhs, double);
        ILL_IFFREE (lp->rhsind, char);
        ILL_IFFREE (lp->rangesname, char);
        ILL_IFFREE (lp->rangesind, char);
        ILL_IFFREE (lp->boundsname, char);
        ILL_IFFREE (lp->lbind, char);
        ILL_IFFREE (lp->ubind, char);
        ILL_IFFREE (lp->lower, double);
        ILL_IFFREE (lp->upper, double);
        ILL_IFFREE (lp->intmarker, char);
        ILL_IFFREE (lp->refrow, char ); 
	ILL_IFFREE (lp->is_sos_member, int ); 
	ILL_IFFREE (lp->sos_weight, double ); 
	ILL_IFFREE (lp->sos_col, int ); 
	ILL_IFFREE (lp->sos_set, sosptr ); 
        ILLinit_rawlpdata (lp, NULL);
    }
}

const char* ILLraw_rowname(rawlpdata *lp, int i)
{
    const char* name = NULL; 
    ILL_FAILfalse_no_rval((i >= 0) && (i < lp->nrows), "index out of range"); 
    ILL_FAILfalse_no_rval(lp->nrows == lp->rowtab.tablesize, 
		  "tab and lp must be in synch");
    name = ILLsymboltab_get(&lp->rowtab, i); 
CLEANUP:
    return name;
} 
const char* ILLraw_colname (rawlpdata *lp, int i) 
{
    const char* name = NULL; 
    ILL_FAILfalse_no_rval((i >= 0) && (i < lp->ncols), "index out of range"); 
    ILL_FAILfalse_no_rval(lp->ncols == lp->coltab.tablesize, 
		  "tab and lp must be in synch");
    name = ILLsymboltab_get(&lp->coltab, i); 
CLEANUP:
    return name;
} 

int ILLraw_add_col(rawlpdata *lp, const char* name, int intmarker) 
{
    int rval =0; 
    int pindex, hit;

    rval = ILLsymboltab_register(&lp->coltab, name, -1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF(rval); 
    if (lp->ncols >= lp->colsize) {
	rval = rval || ILLutil_reallocrus_scale((void **) &lp->cols,
						&lp->colsize, lp->ncols + 1, 1.3,
						sizeof (colptr *)); 
    }
    if (lp->ncols >= lp->intsize) {
	rval = rval || ILLutil_reallocrus_scale((void **) &lp->intmarker,
						&lp->intsize, lp->ncols + 1, 
						1.3, sizeof (char)); 
    }
    if (lp->ncols >= lp->is_sos_size) {
	rval = rval || ILLutil_reallocrus_scale((void **) &lp->is_sos_member,
						&lp->is_sos_size, lp->ncols + 1, 
						1.3, sizeof (int)); 
    }
    ILL_CLEANUP_IF(rval); 
    lp->cols[lp->ncols] = (colptr *) NULL;
    lp->is_sos_member[lp->ncols] = -1; 
    lp->intmarker[lp->ncols] = intmarker; 
    lp->ncols++; 
CLEANUP:
    ILL_RETURN(rval, "ILLraw_add_col"); 
}
  
int ILLraw_init_rhs(rawlpdata *lp) 
{ 
    int i, rval =0;

    ILL_FAILfalse(lp->rhsind == NULL, "Should be called exactly once"); 
    if (lp->nrows > 0) {
	ILL_SAFE_MALLOC(lp->rhsind, lp->nrows, char);
	for (i = 0; i < lp->nrows; i++) {
	    lp->rhsind[i] = (char) 0;
	}
    }
CLEANUP:
    ILL_RETURN(rval, "ILLraw_init_rhs"); 
}

int ILLraw_init_ranges(rawlpdata *lp) 
{
    int i, rval = 0; 

    ILL_FAILfalse(lp->rangesind == NULL, "Should be called exactly once"); 
    if (lp->nrows > 0) {
	ILL_SAFE_MALLOC(lp->rangesind, lp->nrows, char);
	for (i = 0; i < lp->nrows; i++) {
	    lp->rangesind[i] = (char) 0;
	}
    }
CLEANUP: 
    ILL_RETURN(rval, "ILLraw_init_ranges"); 
}

int ILLraw_add_col_coef(rawlpdata *lp, int colind, int rowind, double coef)
{
    colptr *cp = ILLcolptralloc (&lp->ptrworld);
    if (!cp) {  return 1; } 
    cp->this = rowind;
    cp->coef = coef;
    cp->next = lp->cols[colind];
    lp->cols[colind] = cp;
    return 0; 
} 


int ILLraw_add_ranges_coef(rawlpdata *lp, int rowind, double coef)
{
    colptr *cp = ILLcolptralloc (&lp->ptrworld);
    if (!cp) { return 1; } 
    cp->this = rowind;
    cp->coef = coef;
    cp->next = lp->ranges;
    lp->ranges = cp;
    lp->rangesind[rowind] = (char) 1;
    return 0; 
}

int ILLraw_add_sos(rawlpdata *lp, int tp) 
{
    int rval = 0; 
    sosptr *sos, *bef; 

    if (lp->nsos >= lp->sos_setsize) {
	if (ILLutil_reallocrus_scale ((void **) &lp->sos_set,
				      &lp->sos_setsize, lp->nsos + 1, 1.3,
				      sizeof (sosptr *))) { 
            ILL_CLEANUP_IF(rval);
	}
    } 
    sos = lp->sos_set + lp->nsos; 
    sos->nelem = 0;  
    sos->type = tp; 
    if (lp->nsos == 0) {
	sos->first = 0; 
    }  else {
	bef = &(lp->sos_set[lp->nsos-1]); 
	sos->first = bef->first + bef->nelem; 
    } 
    lp->nsos++; 
CLEANUP:
    ILL_RETURN(rval, "ILLraw_add_sos"); 
}

int ILLraw_is_mem_other_sos(rawlpdata *lp, int colind) 
{
    return  (lp->is_sos_member[colind] >= 0) &&
        (lp->is_sos_member[colind] != (lp->nsos-1)); 
}

int ILLraw_add_sos_member(rawlpdata *lp, int colind) 
{
    int rval = 0; 
    ILL_FAILfalse(lp->nsos > 0, 
		  "we should have called ILLraw_add_sos earlier"); 
    ILL_FAILtrue(ILLraw_is_mem_other_sos(lp, colind), 
		 "colind is member of another sos set");
    
    if (lp->is_sos_member[colind] == -1) { 
	if (lp->nsos_member >= lp->sos_weight_size) { 
	    if (ILLutil_reallocrus_scale((void **) &lp->sos_weight,
					 &lp->sos_weight_size, 
					 lp->nsos_member + 1, 
					 1.3, sizeof(double))) { 
                ILL_CLEANUP_IF(rval);
	    }
	}
	if (lp->nsos_member >= lp->sos_col_size) { 
	    if (ILLutil_reallocrus_scale((void **) &lp->sos_col,
					 &lp->sos_col_size, 
					 lp->nsos_member + 1, 
					 1.3, sizeof(int))) { 
                ILL_CLEANUP_IF(rval);
	    }
	}
	lp->sos_col[lp->nsos_member]= colind; 
	lp->sos_set[lp->nsos-1].nelem++; 
	lp->is_sos_member[colind] = lp->nsos-1; 
	lp->nsos_member++; 
    }
CLEANUP:
    ILL_RETURN(rval, "ILLraw_add_sos_member"); 
}


int ILLraw_add_row(rawlpdata *lp, const char* name, char sense, double rhs) 
{
    int pindex, hit, rval = 0; 

    rval = ILLsymboltab_register(&lp->rowtab, name, -1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF(rval); 
    if (lp->nrows >= lp->sensesize) {
	if (ILLutil_reallocrus_scale ((void **) &lp->rowsense,
				      &lp->sensesize, lp->nrows + 1, 
				      1.3, sizeof (char))) {
	    ILL_CLEANUP_IF(rval);
	}
    }
    if (lp->nrows >= lp->rhssize) {
	if (ILLutil_reallocrus_scale ((void **) &lp->rhs,
				      &lp->rhssize, lp->nrows + 1, 1.3, 
				      sizeof (double))) { 
	    ILL_CLEANUP_IF(rval);
	}
    }
    lp->rowsense[lp->nrows] = sense;
    lp->rhs[lp->nrows] = rhs;
    lp->nrows++; 

CLEANUP:
    ILL_RETURN(rval, "ILLraw_add_row"); 
} 

static int ILLcheck_rawlpdata(rawlpdata *lp)
{
    int i, col, rval = 0; 
    int si, *perm = NULL; 
    const char *c1, *c2;
    sosptr *set; 

    ILL_FAILfalse(lp, "lp must not be NULL");  

    /* check    
     *          *) that there is at least one variable 
     *          *) that the weights in all SOS sets are distinct 
     *          *) all sos members are non integer variables     
     *          *) sos set members have distint weights 
     *          *) objindex is not -1
     *          *) INVARIANT: rowname[objindex] != NULL
     *          *) INVARIANT: upper/lower arrays are filled in 
     *          *) INVARIANT: if col or rownames != NULL then 
     *                        all their elements are not NULL
     */
    if (lp->ncols < 1) { 
        return ILLdata_error(lp->error_collector, "There are no variables.");
    }
    if (lp->objindex == -1) {
        return ILLdata_error(lp->error_collector, "There is no objective fct.");
    } 
    ILL_FAILfalse(ILLraw_rowname(lp, lp->objindex) != NULL, 
		  "must have objective name");
    if (lp->nsos_member > 1) { 
	ILL_SAFE_MALLOC(perm, lp->nsos_member, int); 
	for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si; 
            for (i=0; i < set->nelem; i++) {
                col = lp->sos_col[set->first + i]; 
		if (lp->intmarker[col]) { 
		    rval = ILLdata_error(lp->error_collector, 
					"SOS set member \"%s\" is an %s.\n", 
					 ILLraw_colname(lp, col), 
					 "integer/binary variable");
		} 
            } 
	    if (set->nelem > 1) {
                for (i=0; i < set->nelem; i++) {
		    perm[i] = set->first + i; 
		}
		ILLutil_double_perm_quicksort(perm, lp->sos_weight, 
                                              set->nelem); 
                for (i=1; i < set->nelem; i++) {
		    if (lp->sos_weight[perm[i-1]] == 
                        lp->sos_weight[perm[i]]) {
                        c1 = ILLraw_colname (lp, lp->sos_col[perm[i]]);
                        c2 = ILLraw_colname (lp, lp->sos_col[perm[i-1]]);
			ILLdata_error(lp->error_collector, 
					"\"%s\" and \"%s\" both have %s %f.\n", c1, c2,
				    "SOS weight", 
				    lp->sos_weight[perm[i]]); 
			rval = 1; 
		    } 
                } 
	    } 
	}
    }
    for (i = 0; i < lp->ncols; i++) {
	ILL_CHECKnull(ILLraw_colname(lp, i), "There is a NULL col name"); 
    } 
    for (i = 0; i < lp->nrows; i++) {
	ILL_CHECKnull(ILLraw_rowname(lp, i), "There is a NULL row name"); 
    } 
    ILL_FAILtrue((lp->upper == NULL) | (lp->lower == NULL), 
		 "Upper/Lower arrays must be filled in.");
   
    rval += ILLraw_check_bounds(lp); 
CLEANUP:
    ILL_IFFREE(perm, int); 
    ILL_RESULT(rval, "ILLcheck_rawlpdata"); 
}

int ILLraw_init_bounds(rawlpdata *lp) 
{
    int i, rval = 0; 

    ILL_FAILfalse(lp->upper == NULL, "Should be called exactly once"); 
    ILL_FAILfalse(lp->lower == NULL, "Should be called exactly once"); 
    ILL_FAILfalse(lp->lbind == NULL, "Should be called exactly once"); 
    ILL_FAILfalse(lp->ubind == NULL, "Should be called exactly once"); 

    ILL_SAFE_MALLOC(lp->upper, lp->ncols, double);
    ILL_SAFE_MALLOC(lp->lower, lp->ncols, double);
    ILL_SAFE_MALLOC(lp->lbind, lp->ncols, char);
    ILL_SAFE_MALLOC(lp->ubind, lp->ncols, char);

    for (i = 0; i < lp->ncols; i++) {
	lp->lbind[i] = (char) 0;
	lp->ubind[i] = (char) 0;
	lp->lower[i] = 0.0;
    }
CLEANUP:
    ILL_RETURN(rval, "ILLraw_init_bounds");
}

const char* ILLraw_set_lowerBound(rawlpdata *lp, int i, double bnd) 
{
    ILL_FAILtrue_no_rval(i >= lp->ncols, "proper colind"); 
    if (lp->lbind[i]) {
	return "Using previous bound definition."; 
    } 
    lp->lower[i] = bnd;
    lp->lbind[i] = (char) 1; 
CLEANUP:
    return NULL; 
}

const char* ILLraw_set_upperBound(rawlpdata *lp, int i, double bnd) 
{
   ILL_FAILtrue_no_rval(i >= lp->ncols, "proper colind"); 
    if (lp->ubind[i]) {
	return "Using previous bound definition."; 
    } 
    lp->upper[i] = bnd;
    lp->ubind[i] = (char) 1;
    if  (lp->lower[i] == 0.0 && bnd == 0.0) {
	return "0.0 upper bound fixes variable.";
    }
CLEANUP:
    return NULL;
} 

const char* ILLraw_set_fixedBound(rawlpdata *lp, int i, double bnd) 
{
    ILL_FAILtrue_no_rval(i >= lp->ncols, "proper colind"); 
    if (lp->ubind[i] || lp->lbind[i]) {
	return "Using previous bound definition."; 
    } 
    lp->lower[i] = bnd;
    lp->lbind[i] = (char) 1;
    lp->upper[i] = bnd;
    lp->ubind[i] = (char) 1;
CLEANUP:
    return NULL;
}

const char* ILLraw_set_unbound(rawlpdata *lp, int i) 
{
   ILL_FAILtrue_no_rval(i >= lp->ncols, "proper colind"); 
   if (lp->lbind[i] || lp->ubind[i]) {
	return "Using previous bound definition."; 
   }
   lp->lower[i] = -ILL_MAXDOUBLE; 
   lp->upper[i] = ILL_MAXDOUBLE;
   lp->lbind[i] = 1; 
   lp->ubind[i] = 1;
CLEANUP: 
   return NULL;
}

const char* ILLraw_set_binaryBound(rawlpdata *lp, int i) 
{
   ILL_FAILtrue_no_rval(i >= lp->ncols, "proper colind"); 
   if (lp->lbind[i] || lp->ubind[i]) {
	return "Using previous bound definition."; 
   }
   lp->lower[i] = 0.0;
   lp->upper[i] = 1.0; 
   lp->lbind[i] = 1; 
   lp->ubind[i] = 1;
CLEANUP: 
   return NULL;
}

int ILLraw_fill_in_bounds(rawlpdata *lp) 
{
    int rval = 0, i; 
    if (lp->lbind == NULL) {
        ILLraw_init_bounds(lp);
    } 
    ILL_FAILtrue(lp->upper == NULL, "must all be there now");
    ILL_FAILtrue(lp->lower == NULL, "must all be there now");
    ILL_FAILtrue(lp->lbind == NULL, "must all be there now");
    ILL_FAILtrue(lp->ubind == NULL, "must all be there now");
    for (i = 0; i < lp->ncols; i++) {
        if (!lp->lbind[i]) { 
            if (lp->ubind[i] && (lp->upper[i] < 0.0)) {
	        lp->lower[i] = -ILL_MAXDOUBLE; 
            } 
	}
        if (!lp->ubind[i]) { 
            /* int vars without bounds are binary                        */
            /* all, also int vars                                        */
            /*          with explicit lower bound 0.0 are in [0.0,+inf]  */
	    if (((lp->intmarker != NULL) && lp->intmarker[i]) &&         
                 !lp->lbind[i]) { 
		lp->upper[i] = 1.0;
	    } else {
		lp->upper[i] = ILL_MAXDOUBLE; 
	    }
	}
    }

CLEANUP: 
    if (rval) { 
	ILL_IFFREE(lp->lower, double); 
	ILL_IFFREE(lp->upper, double); 
    }
    ILL_RETURN(rval, "ILLraw_fill_in_bounds"); 
} 

static int ILLraw_check_bounds(rawlpdata *lp) 
{
    int rval = 0, i; 
    ILL_FAILtrue(lp->upper == NULL, "must all be there now");
    ILL_FAILtrue(lp->lower == NULL, "must all be there now");
    ILL_FAILtrue(lp->lbind == NULL, "must all be there now");
    ILL_FAILtrue(lp->ubind == NULL, "must all be there now");
    for (i = 0; i < lp->ncols; i++) {
       if (lp->upper[i] < lp->lower[i])  {
           rval += ILLdata_error(lp->error_collector, 
						"Lower bound is bigger than %s \"%s\".\n", 
                        "upper bound for", ILLraw_colname(lp, i)); 
       } 
    }
    ILL_RESULT(rval, "ILLraw_check_bounds");
CLEANUP:
    ILL_RETURN(rval, "ILLraw_check_bounds");
} 

int ILLraw_first_nondefault_bound(ILLlpdata *lp) 
{
    int ri = lp->nstruct, i; 
    ILL_FAILtrue_no_rval(lp->lower == NULL || lp->upper == NULL, 
		 "Should not call write_bounds when lower or upper are NULL"); 
    for (ri = 0; ri < lp->nstruct; ri++) {
        i = lp->structmap[ri];
	if (!ILLraw_default_lower(lp, i) || !ILLraw_default_upper(lp, i)) 
	    break;
    }
CLEANUP: 
    return ri;
} 

int ILLraw_default_lower(ILLlpdata *lp, int i) 
{
    ILL_FAILtrue_no_rval(lp->lower == NULL || lp->upper == NULL, 
		 "Should not call write_bounds when lower or upper are NULL"); 
    ILL_FAILfalse_no_rval(lp->ncols > i, "i is not col index"); 
    if ((lp->lower[i] == 0.0) && (lp->upper[i] >= 0.0)) {
      return 1; 
    } 
    if ((lp->lower[i] == -ILL_MAXDOUBLE) && (lp->upper[i] < 0.0)) {
      return 1;
    } 
CLEANUP: 
    return 0; 
} 

int ILLraw_default_upper(ILLlpdata *lp, int i) 
{
    int isInt; 

    ILL_FAILtrue_no_rval(lp->lower == NULL || lp->upper == NULL, 
		 "Should not call write_bounds when lower or upper are NULL"); 
    ILL_FAILfalse_no_rval(lp->ncols >= i, "i is not col index"); 
    isInt = (lp->intmarker != NULL) && lp->intmarker[i];
    if (isInt) { 
       if (lp->lower[i] == 0.0)  {
           return (lp->upper[i] == 1.0); 
       }
    }
    
    if (lp->upper[i] == ILL_MAXDOUBLE) { 
       return 1; 
    }
CLEANUP:
    return 0;
}

int ILLraw_fill_in_rownames(rawlpdata *lp)
{
    int i, rval = 0; 
    char uname[ILL_namebufsize];
    ILLsymboltab *rowtab; 
    char first = 1; 

    rowtab = &lp->rowtab; 
    ILL_FAILtrue(lp->nrows != rowtab->tablesize,
                 "must have same #entries"); 
    for (i = 0; (rval == 0) && i < lp->nrows; i++) {
        if (ILLsymboltab_get(rowtab, i) == NULL) { 
	    if (first) {
		ILLdata_warn(lp->error_collector, 
			    "Generating names for unnamed rows.");
                first = 0; 
	    }

            ILLsymboltab_unique_name(rowtab, i, "c",  uname); 
	    rval = ILLsymboltab_rename(rowtab, i, uname); 
	    ILL_CLEANUP_IF(rval); 
	} 
    }
CLEANUP:
    ILL_RESULT(rval, "ILLraw_fill_in_rownames"); 
} 

static int whichColsAreUsed(rawlpdata *raw, ILLlpdata *lp, int *colindex) 
{
    int rval = 0; 
    int i, objind = raw->objindex; 
    colptr *cp; 
    char *colUsed = NULL; 

    /* colUsed[i]  variable raw->colnames[i] is used in obj fct 
     * and/or equation(s) */
    ILL_SAFE_MALLOC(colUsed, raw->ncols, char);
    for (i=0; i < raw->ncols; i++) {
	colUsed[i] = 0; 
    }
    for (i = 0; i < raw->ncols; i++) {
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    if ((cp->this == objind) || 
		(raw->rowsense[cp->this] != 'N')) { 
		colUsed[i] = 1; 
		break;
	    }
	}
    }

    /* colindex[i] = -1 for undefined, 0, 1, ... lp->ncol-1 
     * lp->ncols <= raw->ncols */
    for (i = 0; i < raw->ncols; i++) {
        if (colUsed[i]) { 
	    colindex[i] = lp->ncols++; 
        } else { 
	    colindex[i] = -1; 
            ILLdata_warn(raw->error_collector, 
				"\"%s\" is used in non objective 'N' rows only.",
			 ILLraw_colname(raw, i)); 
	}
    }
    if (lp->ncols < 1) { 
        rval = ILLdata_error(raw->error_collector, 
			"There are no variables.");
        ILL_CLEANUP_IF(rval); 
    }
CLEANUP: 
    ILL_IFFREE(colUsed, char); 
    ILL_RESULT(rval, "whichColsAreUsed"); 
}

static int whichRowsAreUsed(rawlpdata *raw, ILLlpdata *lp, 
                             int *rowindex) 
{
    int i, rval = 0;

    /* only use non 'N' rows */
    for (i = 0; i < raw->nrows; i++) {
        if (raw->rowsense[i] != 'N') {
            rowindex[i] = lp->nrows++;
        } else {
            rowindex[i] = -1; 
        }
    }
    if (lp->nrows == 0) { 
        rval = ILLdata_error(raw->error_collector, "There are no constraints.");
    }
    ILL_RESULT(rval, "whichRowsAreUsed"); 
}


static int transferObjective(rawlpdata *raw, ILLlpdata *lp, int *colindex) 
{
    int rval = 0, i, ci, objind = raw->objindex;
    colptr *cp;
    int* coefWarn = NULL; 

    /* transfer objective fct */
    ILL_SAFE_MALLOC(lp->obj, lp->ncols, double);
    ILL_SAFE_MALLOC(coefWarn, lp->ncols, int);
    for (i = 0; i < lp->ncols; i++) {
	lp->obj[i] = 0.0;
	coefWarn[i] = 0; 
    }
    for (i = 0; i < raw->ncols; i++) {
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    if (cp->this == objind) { 
		ci = colindex[i]; 
		ILL_FAILfalse(ci != -1, 
			      "all vars in obj fct should be marked as useful"); 
		coefWarn[ci]++; 
		lp->obj[ci] += cp->coef;
		if (coefWarn[ci] == 2) {
		    ILLdata_warn(raw->error_collector, 
				"Multiple coefficients for \"%s\" in %s.",
				 ILLraw_colname(raw, i), 
				 "objective function"); 
		} 
	    }
	}
    }
CLEANUP: 
    ILL_IFFREE(coefWarn, int); 
    ILL_RETURN(rval, "transferObjective"); 
}

static int transferColNamesLowerUpperIntMarker(rawlpdata *raw, ILLlpdata *lp,
					       int *colindex) 
{
    int i, ci, ind, pre, rval = 0; 
    int hasIntVar; 
    ILL_SAFE_MALLOC(lp->colnames, lp->ncols, char*);
    if (raw->upper)  { 
	ILL_SAFE_MALLOC(lp->upper, lp->ncols, double);
    }
    if (raw->lower) {
	ILL_SAFE_MALLOC(lp->lower, lp->ncols, double);
    }
    ILL_SAFE_MALLOC(lp->intmarker, lp->ncols, char);
    hasIntVar = 0; 
    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i]; 
	if (ci != -1) { 
	    ILL_FAILfalse((ci >= 0) && (ci < lp->ncols), 
			  "colindex problem"); 
	    ILL_UTIL_STR(lp->colnames[ci], ILLraw_colname(raw, i));  
            rval = ILLsymboltab_register(&lp->coltab, 
					 lp->colnames[ci], -1, &ind, &pre); 
            ILL_FAILfalse((rval == 0) && (ind == ci) && (pre == 0), 
			  "should have new entry"); 
	    if (raw->upper) {
		lp->upper[ci] = raw->upper[i];
	    } 
	    if (raw->lower) {
		lp->lower[ci] = raw->lower[i];
	    } 
	    lp->intmarker[ci] = raw->intmarker[i]; 
	    hasIntVar = hasIntVar || lp->intmarker[ci];
            ILL_IFDOTRACE {
                if (lp->lower) { 
		    ILLprt_double(stdout, lp->lower[ci]); 
		    ILL_IFTRACE(" <= "); 
                }
		ILL_IFTRACE("%s", lp->colnames[ci]);
                if (lp->upper) { 
		    ILL_IFTRACE(" <= "); 
		    ILLprt_double(stdout, lp->upper[ci]); 
                }
                if (lp->intmarker[ci]) {
		    ILL_IFTRACE(" INTEGER "); 
                } 
		ILL_IFTRACE("\n"); 
            }
	}
    } 
    if (!hasIntVar) {
	ILL_IFFREE(lp->intmarker, char); 
    } 
CLEANUP: 
    ILL_RETURN(rval, "transferColNamesLowerUpperIntMarker"); 
}

static void safeRegister(ILLsymboltab *tab, const char* name, int i) 
{
    int ind, pre, rval;
    rval  = ILLsymboltab_register(tab, name, -1, &ind, &pre); 
    ILL_FAILfalse((rval == 0) && (ind == i) && (pre == 0), 
			      "Pgming Error: should have new entry"); 
CLEANUP: 
    return;
} 

static int transferSenseRhsRowNames(rawlpdata *raw, ILLlpdata *lp,
				    int *rowindex) 
{ 
    int i, ri, rval = 0; 
    int objind = raw->objindex; 

    /* transfer sense/rhs/rownames */
    if (lp->nrows > 0) {
        ILL_SAFE_MALLOC(lp->sense, lp->nrows, char);
        ILL_SAFE_MALLOC(lp->rhs, lp->nrows, double);
        ILL_SAFE_MALLOC(lp->rownames, lp->nrows, char*);

        ILL_FAILfalse(ILLraw_rowname(raw, raw->objindex), "NULL objname"); 
        safeRegister(&lp->rowtab, ILLraw_rowname(raw, raw->objindex), 0);

        ri = 0;
        for (i = 0; i < raw->nrows; i++) {
            ri = rowindex[i]; 
            if (i == raw->refrowind) {
                ILL_UTIL_STR(lp->refrowname, ILLraw_rowname(raw, i)); 
		lp->refind = ri; 
            } 
            if (raw->rowsense[i] != 'N') {
                ILL_FAILfalse(ILLraw_rowname(raw, i) != NULL,
			      "all rownames should be non NULL");
                ILL_UTIL_STR(lp->rownames[ri], ILLraw_rowname(raw, i)); 
                safeRegister(&lp->rowtab, lp->rownames[ri], ri + 1); 
                lp->sense[ri] = raw->rowsense[i];
                lp->rhs[ri] = raw->rhs[i];
            } else if (i == objind) {
                ILL_FAILfalse(lp->objname == NULL, "objname == NULL"); 
                ILL_UTIL_STR(lp->objname, ILLraw_rowname(raw, i));
            } else {
                /* unused 'N' row */
            }
        }
        ILL_FAILfalse((lp->nrows + 1) == lp->rowtab.tablesize, 
                                       "problem with rowtab structure");
    }
CLEANUP: 
    ILL_RETURN(rval, "transferSenseRhsRowNames"); 
}

static int buildMatrix(rawlpdata *raw, ILLlpdata *lp, 
		       int *rowindex, int *colindex) 
{
    int i, ri, ci, k, nempty = 0,  rval = 0; 
    int *nRowsUsed = (int *) NULL;
    int *coefSet = (int *) NULL;
    int *coefWarn = (int*) NULL; 
    ILLmatrix *A = &lp->A; 
    colptr *cp = NULL; 

    /* put subjective fcts into matrix */ 
    ILL_SAFE_MALLOC(A->matcnt, lp->ncols, int);
    ILL_SAFE_MALLOC(A->matbeg, lp->ncols, int);
    ILL_SAFE_MALLOC(nRowsUsed, lp->nrows, int);

    ILL_SAFE_MALLOC(coefWarn, lp->ncols, int);
    for (i=0; i < lp->ncols; i++) {
	coefWarn[i] = 0; 
    }
    for (i=0; i < lp->nrows; i++) {
	nRowsUsed[i] = -1; 
    }
    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i]; 
	if (ci == -1) continue;
	k = 0;
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    ri = rowindex[cp->this]; 
	    if (ri >= 0) {
		if (nRowsUsed[ri] != i) { 
		    nRowsUsed[ri] = i; 
		    k++;
		} else {
		    if (! coefWarn[ci]) { 
			ILLdata_warn(raw->error_collector, 
				"Multiple coefficients for \"%s\" %s.",
				     lp->colnames[i], "in a row"); 
			coefWarn[ci] = 1;
		    }
		} 
	    }
	}
	A->matcnt[ci] = k;
	A->matbeg[ci] = lp->nzcount + nempty;   /* mark empty cols */
	lp->nzcount += k;
	if (k == 0) nempty++;
    }

    A->matrows    = lp->nrows;
    A->matcols    = lp->ncols;
    A->matcolsize = lp->ncols;
    A->matsize = lp->nzcount + nempty + 1;
    A->matfree = 1;
    ILL_SAFE_MALLOC(A->matind, A->matsize, int);
    ILL_SAFE_MALLOC(A->matval, A->matsize, double);
    ILL_SAFE_MALLOC(coefSet, lp->nrows, int);

    for (k=0; k < lp->nrows; k++) {
        coefSet[k] = -1;  
    } 

    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i];
	if (ci == -1) continue;  /* unused variable */
	k = A->matbeg[ci];
	if (A->matcnt[ci] == 0) {
	    A->matind[k] = 1;  /* Used in addcols and addrows */
	} else {
	    for (cp = raw->cols[i]; cp; cp = cp->next) {
		ri = rowindex[cp->this] ;
		if (ri >= 0) {
		    if (coefSet[ri] == -1) { 
			A->matind[k] = ri;
			A->matval[k] = cp->coef;
			coefSet[ri] = k; 
			k++;
		    } else {
			A->matval[coefSet[ri]] += cp->coef;
		    } 
		}
	    }
	    if (k != A->matbeg[ci] + A->matcnt[ci]) { 
		ILL_ERROR(rval, "problem with matrix"); 
	    }
            for (k--; k >= A->matbeg[ci]; k--) {
                coefSet[A->matind[k]] = -1;
            }
	}
    }
    A->matind[lp->nzcount + nempty] = -1;
CLEANUP: 
    ILL_IFFREE (nRowsUsed, int);
    ILL_IFFREE (coefWarn, int); 
    ILL_IFFREE (coefSet, int);
    ILL_RETURN(rval, "buildMatrix"); 
}

static int transferRanges(rawlpdata *raw, ILLlpdata *lp, int *rowindex)
{
    int i, ri, rval = 0; 
    colptr *cp;

    /*****************************************************/
    /*                                                   */
    /*  Interpretation of RANGE values in MPS files      */
    /*                                                   */
    /*    G    rhs           <= row <= rhs + |range|     */
    /*    L    rhs - |range| <= row <= rhs               */
    /*    E +  rhs           <= row <= rhs + range       */
    /*    E -  rhs + range   <= row <= rhs               */
    /*                                                   */
    /*     - where + and - refer to the sign of range    */
    /*       and the letters refer to sense of the row.  */
    /*                                                   */
    /*    We will store ranged rows as                   */
    /*                                                   */
    /*       rhs  <= row  <= rhs + range                 */
    /*                                                   */
    /*****************************************************/


    ILL_SAFE_MALLOC(lp->rangeval, lp->nrows, double);
    for (i = 0; i < lp->nrows; i++) {
	lp->rangeval[i] = 0.0;
    }
    for (cp = raw->ranges; cp; cp = cp->next) {
        i = cp->this; 
	ri = rowindex[cp->this];
	switch (raw->rowsense[i]) {
	case 'N':
	    ILLdata_error(raw->error_collector, "No range for N-row.\n"); 
	    rval = 1; goto CLEANUP;
	case 'G':
	    lp->sense[ri] = 'R';
	    lp->rangeval[ri] = ILL_OURABS (cp->coef);
	    break;
	case 'L':
	    lp->sense[ri] = 'R';
	    lp->rhs[ri] -= ILL_OURABS (cp->coef);
	    lp->rangeval[ri] = ILL_OURABS (cp->coef); 
	    break;
	case 'E':
	    lp->sense[ri] = 'R';
	    if (cp->coef >= 0.0) {
		lp->rangeval[ri] = cp->coef;
	    } else {
		lp->rhs[ri] += cp->coef;
		lp->rangeval[ri] = -cp->coef;
	    }
	    break;
	}
    }
CLEANUP: 
    ILL_RETURN(rval, "transferRanges"); 
}

static int  initStructmap(ILLlpdata *lp)
{
    int i, rval = 0; 

    /* all vars are structural */
    ILL_SAFE_MALLOC(lp->structmap, lp->nstruct, int);
    for (i = 0; i < lp->nstruct; i++) {
        lp->structmap[i] = i;
    }

CLEANUP: 
    ILL_RETURN(rval, "initStructmap"); 
}
  
static int buildSosInfo(rawlpdata *raw, ILLlpdata *lp, int *colindex) 
{
    int i, ci, set, rval = 0; 
    int nSosMem, nSetMem; 

    /* build sos info */ 
    /* see comment in lpdata.h about ILLlpdata's sos and is_sos_mem 
     * fields and section of ILLprint_rawlpdata that prints SOS sets */

    ILL_SAFE_MALLOC(lp->is_sos_mem, lp->ncols, int); 
    nSosMem = 0; 
    for (i =0; i < raw->ncols; i++) {
	ci = colindex[i]; 
	if (ci !=  -1) { 
	    lp->is_sos_mem[ci] = raw->is_sos_member[i]; 
	    if (raw->is_sos_member[i] != -1) nSosMem++; 
	}
    }
    if (nSosMem > 0) { 
	lp->sos.matsize = nSosMem; 
	lp->sos.matcols = raw->nsos; 
	lp->sos.matcolsize = raw->nsos; 
	lp->sos.matrows = lp->ncols; 
	lp->sos.matfree = 0; 
	ILL_SAFE_MALLOC(lp->sos.matval, nSosMem, double); 
	ILL_SAFE_MALLOC(lp->sos.matind, nSosMem, int); 
	ILL_SAFE_MALLOC(lp->sos.matbeg, raw->nsos, int); 
	ILL_SAFE_MALLOC(lp->sos.matcnt, raw->nsos, int); 
	ILL_SAFE_MALLOC(lp->sos_type, raw->nsos, char); 
        nSosMem = 0; 
	for (set = 0; set < raw->nsos; set++) {
            lp->sos_type[set] = raw->sos_set[set].type;
            lp->sos.matbeg[set] = nSosMem;
	    nSetMem = 0; 
            for (i = raw->sos_set[set].first; 
                 i < raw->sos_set[set].first + raw->sos_set[set].nelem; 
                 i++) {
		ci = colindex[raw->sos_col[i]]; 
		if (ci != -1) { 
		    lp->sos.matind[nSosMem + nSetMem] = ci;
		    lp->sos.matval[nSosMem + nSetMem] = raw->sos_weight[i]; 
                    nSetMem++; 
		}
	     } 
             lp->sos.matcnt[set] = nSetMem; 
             nSosMem += nSetMem; 
	} 
    }
CLEANUP: 
    ILL_RETURN(rval, "buildSosInfo"); 
}

static int convert_rawlpdata_to_lpdata(rawlpdata *raw, ILLlpdata *lp)
/* 
 * only raw's non 'N' rows are converted to matrix entries in lp
 * columns that are used in non objective 'N' rows only are not 
 * converted. That is they don't end up in lp's matrix, row/colnames, 
 * upper/lower bounds or SOS information.
 */
{
    int rval = 0;
    int *rowindex = (int *) NULL;
    int *colindex = (int *) NULL;

    ILL_FAILfalse((raw && lp), "rawlpdata_to_lpdata called without input");
    if (raw->name == NULL) {
        ILLdata_warn(raw->error_collector, "Setting problem name to \"unnamed\"."); 
        ILL_UTIL_STR(raw->name, "unnamed"); 
    }
    rval = ILLcheck_rawlpdata(raw); 
    ILL_CLEANUP_IF(rval); 

    ILL_FAILfalse(raw->objindex != -1, "rawlpdata must have objective fct."); 
    ILLlpdata_init(lp);

    ILL_IFFREE (lp->probname, char);
    lp->probname = raw->name;
    raw->name = (char *) NULL;

    /* MINIMIZE or MAXIMIZE ? */
    lp->objsense = raw->objsense;
    if (lp->objsense != ILL_MIN && lp->objsense != ILL_MAX) {
        ILLdata_error(raw->error_collector, "Bad objsense.\n"); 
        rval = 1; goto CLEANUP;
    }

    ILL_SAFE_MALLOC(colindex, raw->ncols, int);
    ILL_SAFE_MALLOC(rowindex, raw->nrows, int);
    rval = whichColsAreUsed(raw, lp, colindex) || 
	whichRowsAreUsed(raw, lp, rowindex); 
    ILL_CLEANUP_IF(rval); 
    ILL_FAILtrue(lp->ncols == 0 || lp->nrows == 0, "we need rows and cols"); 

    /* array sizes */
    lp->rowsize = lp->nrows;
    lp->colsize = lp->ncols;
    lp->nstruct = lp->ncols;
    lp->structsize = lp->ncols;
    ILLsymboltab_create(&lp->rowtab, lp->nrows); 
    ILLsymboltab_create(&lp->coltab, lp->ncols); 

    rval = transferObjective(raw, lp, colindex); 
    rval = rval || transferColNamesLowerUpperIntMarker(raw, lp, colindex); 
    rval = rval || buildMatrix(raw, lp, rowindex, colindex); 
    rval = rval || buildSosInfo(raw, lp, colindex); 
    ILL_CLEANUP_IF(rval); 
    ILL_IFDOTRACE { ILLmatrix_prt(stdout, &lp->A); }  

    rval = transferSenseRhsRowNames(raw, lp, rowindex); 
    if ((lp->nrows > 0) && raw->ranges) {
	rval = rval || transferRanges(raw, lp, rowindex); 
    }
    ILL_CLEANUP_IF(rval); 

    rval = initStructmap(lp); 
    ILL_CLEANUP_IF(rval); 

CLEANUP:

    ILL_IFFREE (rowindex, int);
    ILL_IFFREE (colindex, int);
    ILLfree_rawlpdata (raw);

    ILL_RESULT(rval, "convert_rawlpdata_to_lpdata");
}

int ILLrawlpdata_to_lpdata (rawlpdata *raw, ILLlpdata *lp)
{
    int rval = 0; 

    ILL_IFDOTRACE { 
	ILLprint_rawlpdata (raw);
    }
    rval = convert_rawlpdata_to_lpdata(raw, lp); 
    if (rval == 0) { 
        rval = ILLlp_add_logicals(lp); 
    }
    ILL_RESULT(rval, "ILLrawlpdata_to_lpdata"); 
} 

static int set_field_name(char **field, const char* name, int *skip) 
{
    int rval = 0; 
    /* name is bounds/rhs/rangesname field from rawlpdata */
    *skip = 0; 
    if (!*field) { 
	ILL_UTIL_STR(*field, name);
    }

    if (strcmp(*field, name)) {
        /* not first specified RHS/BOUNDS - skip it */
	*skip = 1; 
    }
CLEANUP:
    ILL_RETURN(rval, "set_field_name"); 
}

int ILLraw_set_rhs_name(rawlpdata *lp, const char* name, int *skip) 
{
    return set_field_name(&lp->rhsname, name, skip) ;
}
int ILLraw_set_bounds_name(rawlpdata *lp, const char* name, int *skip) 
{
    return set_field_name(&lp->boundsname, name, skip) ;
}
int ILLraw_set_ranges_name(rawlpdata *lp, const char* name, int *skip) 
{
    return set_field_name(&lp->rangesname, name, skip) ;
}

void ILLprint_rawlpdata (rawlpdata *lp)
{
    int i, cnt, si, m;
    char c; 
    double d; 
    colptr *cp;
    sosptr *set;

    if (lp) {
        if (lp->name) {
            printf ("PROBLEM  %s\n", lp->name);
            fflush (stdout);
        }
        if (lp->rowsense && lp->rhs) {
            printf ("Subject To\n");
            for (i = 0; i < lp->nrows; i++) {
		switch (lp->rowsense[i] ) {
		case 'E': c = '='; break;
		case 'L': c = '<'; break;
		case 'G': c = '>'; break;
		default: c = '?'; break;
		} 
                printf ("%s: %c %f\n", ILLraw_rowname(lp, i),  c, lp->rhs[i]);
            }
            printf ("\n"); fflush (stdout);
        }
        if (lp->ncols > 0) {
            printf ("Columns\n");
            for (i = 0; i < lp->ncols; i++) {
                for (cp = lp->cols[i]; cp; cp = cp->next) {
                    printf ("%s: ", ILLraw_rowname(lp, cp->this)); 
                    printf("%c ", (cp->coef < 0.0) ? '-' : '+'); 
                    d = ILL_OURABS(cp->coef); 
                    if (d != 1.0) { 
			printf(" %f ", d); 
                    }
                    printf ("%s\n", ILLraw_colname(lp, i));
                }
                printf ("\n"); fflush (stdout);
            }
        }
        if (lp->rangesname) {
            printf ("RANGES %s\n", lp->rangesname);
            for (cp = lp->ranges; cp; cp = cp->next) {
                printf ("(%s, %f) ", ILLraw_rowname(lp, cp->this), cp->coef);
            }
            printf ("\n"); fflush (stdout);
        }
        if (lp->boundsname) {
            printf ("BOUNDS %s\n", lp->boundsname); fflush (stdout);
        } else {
            printf ("BOUNDS \n"); fflush (stdout);
        } 
        if (lp->lower && lp->upper) {
            for (i = 0; i < lp->ncols; i++) {
                ILLprt_double(stdout, lp->lower[i]); 
                printf(" <= %s <= ", ILLraw_colname(lp, i)); 
                ILLprt_double(stdout, lp->upper[i]); 
                printf("\n"); 
            }
        } 
        if (lp->intmarker) {
	    printf("Integer\n"); 
	    cnt = 0; 
	    for (i = 0; i < lp->ncols; i++) {
		if (lp->intmarker[i]) {
		    printf ("%s", ILLraw_colname(lp, i));
		    cnt++; 
		    if (cnt == 8) {
			printf("\n    ");  
			cnt = 0;
		    } 
		}
	    }
	    printf ("\n"); fflush (stdout);
        }
        printf ("SOS-SETS\n"); 
        for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si; 
	    printf("SOS-SET %d: %s; nelem=%d; first=%d;\n", 
		   si, ((set->type == ILL_SOS_TYPE1) ? "TYPE1" : "TYPE2"), 
		   set->nelem, set->first); 
	    printf("\t"); 
	    for (m = set->first; m < set->first + set->nelem; m++) {
                printf(" %s %f; ", ILLraw_colname(lp, lp->sos_col[m]), 
		       lp->sos_weight[m]); 
	    } 
	    printf("\n"); 
	}
        printf ("\n"); fflush (stdout);
    } 
}

static int ILLmsg(qserror_collector *collector, 
				  int isError, const char* format, va_list args)
{
	const char* pre;
	int slen, errtype;  
	qsformat_error error; 
	char error_desc[256]; 

	vsprintf(error_desc, format, args); 
	slen = strlen(error_desc); 
	if ((slen > 0) && error_desc[slen-1] != '\n') {
		error_desc[slen] = '\n'; 
		error_desc[slen+1] = '\0'; 
	}

	if (collector != NULL) { 
		errtype = (isError) ? QS_DATA_ERROR : QS_DATA_WARN; 
		ILLformat_error_create(&error, errtype, error_desc, 
			-1, NULL, -1); 
		ILLformat_error(collector, &error); 
		ILLformat_error_delete(&error); 
	} else { 
		pre = (isError) ? "Data Error" : "Data Warning";
		fprintf (stderr, "%s: %s", pre, error_desc); 
	}
	return 1; 
}

int ILLdata_error(qserror_collector *collector, const char* format, ...) 
{
	va_list args;
	va_start(args, format);
	return ILLmsg(collector, TRUE, format, args);
}
void ILLdata_warn(qserror_collector *collector, const char* format, ...) 
{
	va_list args;
	va_start(args, format);
	(void) ILLmsg(collector, FALSE, format, args);
}

colptr* ILLcolptralloc(ILLptrworld *p) 
{
    return colptralloc(p); 
}  
