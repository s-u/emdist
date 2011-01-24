/*  RCS_INFO = "$RCSfile: wr_lp.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef WRITE_LP_STATE_H
#define WRITE_LP_STATE_H

/****************************************************************************/
/*                                                                          */
/*               Routines to support writing of LP files                    */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */
#include "iqsutil.h" 

typedef struct ILLwrite_lp_state {
    char buf[ILL_namebufsize];
    char *p; 
    int startlen; 
    int total; 
} ILLwrite_lp_state; 

extern void ILLwrite_lp_state_init(ILLwrite_lp_state *line, const char* str) ;
extern void ILLwrite_lp_state_append(ILLwrite_lp_state *line, const char *str);
extern void ILLwrite_lp_state_append_coef(ILLwrite_lp_state *line, 
                                         double v, int cnt) ;
      /* append number sign ('+', '-') iff cnt > 0 or v < 0.0  
       * append number iff v != 1.0, v != -1.0  
       */
extern void ILLwrite_lp_state_append_number(ILLwrite_lp_state *line, double v); 
extern void ILLwrite_lp_state_save_start(ILLwrite_lp_state *line) ;
extern void ILLwrite_lp_state_start(ILLwrite_lp_state *line) ;

#endif
