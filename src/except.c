/* RCSINFO $Id: exception.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#include "except.h" 
#include <stdio.h> 
#include <string.h> 
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


void ILL_report(const char* msg,  
	       const char* fct, const char* file, unsigned int line, 
	       int with_src_info)
{
    if (msg != NULL) { 
	fprintf(stderr, "FAILURE: %s", msg);
	if (msg[strlen(msg)-1] != '\n') {
	    fprintf(stderr, "\n"); 
	} 
	if (with_src_info == 1) { 
	    fprintf(stderr, "\t"); 
	    if (fct != NULL) {
		fprintf(stderr, "in function %s ", fct); 
	    }
	    fprintf(stderr, "in file %s line %d", file, line); 
	}
	fprintf(stderr, ".\n"); 
    }
}

