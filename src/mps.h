/* RCSINFO $Id: mps.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef MPS_H
#define MPS_H

#include "readline.h" 
#include "format.h" 

/****************************************************************************/
/*                                                                          */
/*              Routines to support Reading and Writing MPS Files           */
/*                                                                          */
/****************************************************************************/

typedef enum { ILL_MPS_NAME, ILL_MPS_OBJSENSE, ILL_MPS_OBJNAME, 
               ILL_MPS_ROWS, ILL_MPS_COLS, ILL_MPS_RHS, ILL_MPS_RANGES,
               ILL_MPS_BOUNDS, ILL_MPS_REFROW, ILL_MPS_ENDATA, 
               ILL_MPS_NONE} ILLmps_section ; 
#define ILL_MPS_N_SECTIONS ILL_MPS_NONE

extern const char* ILLmps_section_name[ILL_MPS_N_SECTIONS + 2];

#include "lpdata.h"
#include "rawlp.h"
#include "read_mps.h"

extern int ILLread_mps (qsline_reader *file,
			const char *filename, rawlpdata *lp);

extern int ILLwrite_mps(ILLlpdata *lp, qserror_collector *collector);	
				/* use lp->reporter for output */

#endif 
