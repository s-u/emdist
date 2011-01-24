/* RCS_INFO = "$Id: reporter.c,v 1.1 2003/11/05 16:49:52 meven Exp $"; */

#include <stdio.h> 
#include "iqsutil.h"
#include "reporter.h"

int ILL_fprintf(void *dest, const char *s)
{
	if (s != NULL) return fprintf((FILE*) dest, s); 
	return 0;
}

void
ILLstring_reporter_init(qsstring_reporter *reporter, 
                        qsreport_string_fct fct, void* dest)
{
    int rval = 0; 
    ILL_FAILfalse(reporter != NULL, "Must get non NULL reporter");
    if (reporter != NULL) { 
	reporter->report_fct = fct; 
	reporter->dest = dest;
    } 
CLEANUP: 
    return;
}

void 
ILLstring_reporter_copy(qsstring_reporter *dest,
                             qsstring_reporter *src)
{
   *dest = *src; 
} 


static int 
string_reporter_main(int ac, char**av) 
{
    int i = 0; 
    qsstring_reporter reporter; 
	
    ILLstring_reporter_init(&reporter, ILL_fprintf, stdout); 
    for (i = 0; i < ac; i++) {
         (void)	ILLstring_report(av[i], &reporter); 
	 (void) ILLstring_report("\n", &reporter); 
    } 
    (void) ILLstring_report(NULL, &reporter); 
	
    return 0;
}

#ifdef REPORTER_MAIN 
int 
main(int ac, char** av) 
{
    return string_reporter_main(ac, av); 
} 
#endif 
