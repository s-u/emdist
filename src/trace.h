/* RCSINFO $Id: trace.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_trace_h 
#define ILL_trace_h 

#include <stdio.h> 

/* users of these macros must declare a static int TRACE variable */
#ifndef NDEBUG
#define ILL_IFTRACE        if (TRACE) printf 
#define ILL_IFTRACE2       if (TRACE > 1) printf 
#define ILL_IFDOTRACE      if (TRACE) 
#else 
/* the optimizer will take care of this */
#define ILL_IFTRACE        if (0) printf 
#define ILL_IFTRACE        if (0) printf 
#define ILL_IFDOTRACE      if (0)
#endif 

#endif 
