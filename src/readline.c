/* RCS_INFO = "$Id: line_reader.c,v 1.2 2003/11/05 16:49:52 meven Exp $"; */

#include <stdio.h> 
#include "iqsutil.h"
#include "readline.h"


/* #ifdef _WINDOWS*/
qsline_reader* 
ILLline_reader_new(qsread_line_fct fct, void* data_src)
{
    qsline_reader* reader; 
    int rval = 0; 
    ILL_NEW(reader, qsline_reader); 
    if (reader != NULL) { 
	reader->read_line_fct = fct; 
	reader->data_src = data_src;
	reader->error_collector = NULL; 
    } 
CLEANUP:
    return reader;
}

void 
ILLline_reader_free(qsline_reader *reader)
{
   ILL_IFFREE(reader, qsline_reader);
} 
/* #endif */

