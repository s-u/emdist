/* RCSINFO $Id: format.h,v 1.3 2003/11/05 16:59:48 meven Exp $ */
#ifndef QS_FORMAT_ERROR_H
#define QS_FORMAT_ERROR_H

#include <stdio.h> 
#include "qsopt.h"

/****************************************************************************/
/*
   The LP/MPS readers, writers, 
       ILLrawlpdata_to_lpdata, and 
   use ILLformat_error to report problems with their input iff
       the line reader used in reading the problem  or 
       the  qserror_collector pointer passed to ILLwrite_lp_file
   is not NULL.

   The QSgui code uses this feature to collect qsformat_error instances 
   which it uses after reading is done to insert error messages into the 
   input window. 
*/
/****************************************************************************/

/* 
for error type USE: 
          QS_DATA_ERROR			
          QS_DATA_WARN			
          QS_MPS_FORMAT_ERROR		
          QS_MPS_FORMAT_WARN		
          QS_LP_FORMAT_ERROR		
          QS_LP_FORMAT_WARN		
          QS_LP_OBJ_WARN			
          QS_GENERIC_ERROR		
*/ 

typedef struct qsformat_error { 
    int type;              
    char* desc; 
    int lineNumber;   /* 1 based line counting */
    char* theLine; 
    int at; 
    struct qsformat_error *next;
} qsformat_error; 

extern int 
    ILLformat_error_create(qsformat_error *error, int mode, const char* desc,
			   int lineNum, const char* theLine, int atPos); 
extern void 
    ILLformat_error_delete(qsformat_error *error); 

extern void
    ILLformat_error_print(FILE *out, qsformat_error *e); 



/*****************************************************************************
 * collecting error messages 
 * either with defining own qsad_error_fct and corresponding data structure 
 * or by using predefined ILLadd_error_to_memory fct with qserror_memory
 */ 

typedef int (*qsadd_error_fct)(void* dest, const qsformat_error *error); 

typedef struct qserror_collector 
{
    qsadd_error_fct add_error; 
    void *dest; 
}  qserror_collector;

typedef struct qserror_memory
{
    char has_error[QS_INPUT_NERROR]; 
    qsformat_error *error_list; 
    unsigned int nerror; 
    char hasErrorLines; 
} qserror_memory;


extern qserror_collector*
    ILLerror_collector_new(qsadd_error_fct  fct, void* dest);

qserror_collector*
    ILLerror_memory_collector_new(qserror_memory *dest);

extern void 
    ILLerror_collector_free(qserror_collector *c); 

#define ILLformat_error(collector, error)  \
	((collector)->add_error((collector)->dest, error))


extern int ILLadd_error_to_memory(void *dest, const qsformat_error *error); 

extern qserror_memory* ILLerror_memory_create(char takeErrorLines);
extern void ILLerror_memory_free(qserror_memory * mem); 

#endif 
