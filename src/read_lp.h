/*  RCS_INFO = "$RCSfile: read_lp_state.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef READ_LP_STATE_H
#define READ_LP_STATE_H

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading LP Files                       */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "iqsutil.h"
#include "readline.h"

typedef struct ILLread_lp_state {
    qsline_reader *file; 
    const char* file_name; 
    int   interactive; 
    int  line_num;
    char realline[ILL_namebufsize];
    char line[ILL_namebufsize];
    char field[ILL_namebufsize+1];
    char fieldOnFirstCol; 
    char *p; 
    char eof; 

    /* these fields aare set as sideeffects to reader routines */
    char sense_val; 
    double bound_val; 
    int  column_index;
} ILLread_lp_state; 

extern int  ILLread_lp_state_init(ILLread_lp_state *state, 	
				  qsline_reader *file, 
				  const char* fname, int interactve);
extern int ILLread_lp_state_next_line(ILLread_lp_state *state); 
extern int  ILLread_lp_state_next_var(ILLread_lp_state *state); 
extern int  ILLread_lp_state_keyword(ILLread_lp_state *state, const char **kwd);
extern int  ILLread_lp_state_bad_keyword(ILLread_lp_state *state); 
extern int  ILLtest_lp_state_keyword(ILLread_lp_state *state, const char* kwd[]);
extern int  ILLread_lp_state_next_field(ILLread_lp_state *state);  
extern int  ILLread_lp_state_next_field_on_line(ILLread_lp_state *state);  
extern void ILLread_lp_state_prev_field(ILLread_lp_state *state);  
extern int  ILLread_lp_state_sign(ILLread_lp_state *state, double *sign); 
extern int  ILLread_lp_state_possible_coef(ILLread_lp_state *state, 
					double *coef, double defValue); 
                                        /* returns 1 iff found a number 
                                         * otherwise 0 */
extern int  ILLread_lp_state_possible_bound_value(ILLread_lp_state *state);
                                               /* returns 1 iff found a number 
                                                * otherwise 0 */
extern int  ILLread_lp_state_colon(ILLread_lp_state *state); 
extern int  ILLread_lp_state_has_colon(ILLread_lp_state *state); 
extern int  ILLread_lp_statxe_has_colon(ILLread_lp_state *state); 
extern int  ILLread_lp_state_next_constraint(ILLread_lp_state *state); 
extern int  ILLread_lp_state_sense(ILLread_lp_state *state); 
extern int  ILLtest_lp_state_sense(ILLread_lp_state *state, char all);
extern void  ILLtest_lp_state_bound_sense(ILLread_lp_state *state); 
extern int  ILLread_lp_state_value(ILLread_lp_state *state, double *d); 
extern int  ILLtest_lp_state_next_is(ILLread_lp_state *state, const char *str);
extern int  ILLread_lp_state_skip_blanks(ILLread_lp_state *state, char wrapLines); 

extern int  ILLcheck_subject_to(ILLread_lp_state *state); 

/*---------------------------------------------------------------------------*/
/* errors and warnings 
 */
extern int  ILLlp_error(ILLread_lp_state *state, const char* format, ...) ; 
extern void ILLlp_warn(ILLread_lp_state *state, const char* format, ...) ; 

/*---------------------------------------------------------------------------*/
/* shared with read_mps_state.c 
 */
extern int ILLget_value (char *line, double *coef) ;

#endif 
