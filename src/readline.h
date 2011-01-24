/* RCSINFO $Id: rdline.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef LINE_READER_FILE_H
#define LINE_READER_FILE_H

#include "qsopt.h" 

/* #ifdef _WINDOWS */
typedef char* (*qsread_line_fct)(char* s, int size, void *src); 

typedef struct qsline_reader
{
    qsread_line_fct read_line_fct; 
    void *data_src; 
    struct qserror_collector *error_collector; 
} qsline_reader;

qsline_reader* ILLline_reader_new(qsread_line_fct fct, void* data_src); 
void ILLline_reader_free(qsline_reader *reader); 

#define ILLline_reader_get(s, size, reader)  \
	(reader)->read_line_fct(s, size, (reader)->data_src) 
	                 /* used by parsers to retrieve next input line */
/* #else  
 *
 * typedef FILE qsline_reader; 
 *
 * #define ILLline_reader_new(fct, data) ((FILE*) (data)) 
 * #define ILLline_reader_free(reader)  
 * #define ILLline_reader_get(s, size, reader) fgets(s,size,reader) 
 * #endif 
 */

#endif  
