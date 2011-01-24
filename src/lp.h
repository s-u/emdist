/* RCSINFO $Id: lp.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef LP_H
#define LP_H

#include "readline.h" 

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading and Writing LP Files           */
/*                                                                          */
/****************************************************************************/
/* 
 * -) anything after '\' is comment 
 * -) Problem is optional and comes first
 * -) Minimize, Maximize, ... comes after Problem
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "lpdata.h"
#include "rawlp.h"
#include "read_lp.h"
#include "write_lp.h"

extern int ILLread_lp(qsline_reader *file, 
		      const char *fname, rawlpdata *lp) ;
extern int ILLwrite_lp(ILLlpdata *l, qserror_collector *collector) ; 
			/* write using current lp->reporter */
extern int ILLis_lp_name_char(char c, int pos);  

extern int ILLread_constraint_name(ILLread_lp_state *state, char** rowname); 
extern int ILLread_one_constraint(ILLread_lp_state *state, const char* rowname,
				  rawlpdata *lp, int allowNewColsAddRow);
extern int ILLread_constraint_expr(ILLread_lp_state *state, rawlpdata *lp, 
			        int rowind, int allowNew);

#endif 
