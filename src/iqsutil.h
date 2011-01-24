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
/* RCSINFO $Id: config.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef __CONFIG_H
#define __CONFIG_H

/* INCLUDE/config.h.  Generated automatically by configure.  */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/


/* Define if your compiler is missing the appropriate function prototype */

/* #undef ILL_PROTO_PRINTF */
/* #undef ILL_PROTO_GETHOSTNAME */
/* #undef ILL_PROTO_GETRUSAGE */

/* Define if you want to use posix threads */
/* #undef ILL_POSIXTHREADS */

/* Define to empty if the keyword `const' does not work.  */
/* #undef const */

/* Define to `int' if <sys/types.h> doesn't define.  */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define to `unsigned char' if <sys/types.h> doesn't define.  */
/* #undef u_char */

/* Define to `int' if the builtin type `void' does not work.  */
/* #undef void */

/* Define if you have the gethostname function.  */

/* Define if you have the socket function.  */

/* Define if you have the strdup function.  */

/* Define if you have the getrusage function.  */
#ifndef WIN32
#ifndef HAVE_GETRUSAGE
#define HAVE_GETRUSAGE
#endif
#endif 

/* Define if you have the times function.  */

/* Define if you have the clock function.  */
#ifndef HAVE_CLOCK
#define HAVE_CLOCK
#endif

/* Define if you have the sleep function.  */

/* Define if you have the <stdio.h> header file.  */
#define HAVE_STDIO_H 1

/* Define if you have the <stdarg.h> header file.  */
#define HAVE_STDARG_H 1

/* Define if you have the <stdlib.h> header file.  */
#define HAVE_STDLIB_H 1

/* Define if you have the <math.h> header file.  */
#define HAVE_MATH_H 1

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1


/* Define if you have the <strings.h> header file.  */
#ifndef WIN32
#define HAVE_STRINGS_H 1
#endif

/* Define if you have the <errno.h> header file.  */
#define HAVE_ERRNO_H 1

/* Define if you have the <assert.h> header file.  */
#define HAVE_ASSERT_H 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */

/* Define if you have the <sys/time.h> header file.  */

/* Define if you have the <time.h> header file.  */
#ifndef HAVE_TIME_H
#define  HAVE_TIME_H
#endif

/* Define if you have the <stddef.h> header file.  */
#define HAVE_STDDEF_H 1

/* Define if you have the <unistd.h> header file.  */
#ifndef WIN32
#ifndef HAVE_UNISTD_H
#define HAVE_UNISTD_H
#endif
#endif

/* Define if you have the <malloc.h> header file.  */
/* #define HAVE_MALLOC_H 1 */

/* Define if you have the <sys/types.h> header file.  */

/* Define if you have the <sys/stat.h> header file.  */

/* Define if you have the <sys/resource.h> header file.  */
#ifndef WIN32
#ifndef HAVE_SYS_RESOURCE_H
#define HAVE_SYS_RESOURCE_H
#endif
#endif

/* Define if you have the <fcntl.h> header file.  */

/* Define if you have the <signal.h> header file.  */

/* Define if you have the <sys/socket.h> header file.  */

/* Define if you have the <netdb.h> header file.  */

/* Define if you have the <netinet/in.h> header file.  */

/* Define if you have the <sys/param.h> header file.  */

/* Define if you have the <sys/times.h> header file.  */

/* Define if your compiler supports attribute modifiers  */
/* such as __attribute__ ((unused)) (gcc 2.8.1 does)     */
#define ILL_ATTRIBUTE 1

/* Some machine (o/s) specific problems */

/* Define if unistd.h uses __vfork but does not prototype it */
/* This happens under Irix 6 */
/* #undef ILL_PROTO___VFORK */

#ifdef WIN32
#define HAVE_DOS_TIME
#endif 

/* __CONFIG_H */
#endif 
/* RCSINFO $Id: machdefs.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __MACHDEFS_H
#define __MACHDEFS_H

#ifdef HAVE_STDIO_H
# include <stdio.h>
#endif

#ifdef HAVE_STDARG_H
# include <stdarg.h>
#endif
#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_STRING_H
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
#  include <strings.h>
#endif
#ifdef HAVE_ERRNO_H
# include <errno.h>
#endif
#ifdef HAVE_ASSERT_H
# include <assert.h>
#endif
#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  ifdef HAVE_TIME_H
#   include <time.h>
#  endif
# endif
#endif
#ifdef HAVE_STDDEF_H
# include <stddef.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_MALLOC_H
# include <malloc.h>
#endif
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h> 
#endif
#ifdef HAVE_FCNTL_H
# include <fcntl.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
# include <sys/socket.h>
#endif
#ifdef HAVE_NETDB_H
# include <netdb.h>
#endif
#ifdef HAVE_NETINET_IN_H
# include <netinet/in.h>
#endif

#ifdef HAVE_SOCKET
#define ILL_NETREADY
#endif

#ifdef USEDMALLOC
#endif

#ifdef ILL_ATTRIBUTE
#define ILL_UNUSED __attribute__ ((unused))
#else
#define ILL_UNUSED
#endif

#ifdef ILL_PROTO_PRINTF
/* assume that if you're missing printf, you're missing a bunch */
extern int
    printf (const char *, ...),
    fprintf (FILE *, const char *, ...),
    fflush (FILE *),
    scanf (const char *, ...),
    sscanf (const char *, const char *, ...),
    fscanf (FILE *, const char *, ...),
    fclose (FILE *),
    ungetc (int, FILE *),
    _filbuf (FILE *),
    time (int *);
#ifdef ILL_NETREADY
extern int
    socket (int, int, int),
    connect (int, const struct sockaddr *, int),
    accept (int, struct sockaddr *, int *),
    bind (int, const struct sockaddr *, int),
    listen (int, int);
#endif
extern void
   *memset (void *, int, size_t),
    perror (const char *);
#endif

#ifdef ILL_PROTO_RENAME
extern int
    rename (const char *, const char *);
#endif

#ifdef ILL_PROTO_GETHOSTNAME
extern int
    gethostname (char *, int);
#endif

#ifdef ILL_PROTO_GETRUSAGE
extern int
    getrusage (int, struct rusage *);
#endif

#ifdef ILL_PROTO___VFORK
extern pid_t
    __vfork (void);
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifndef INT_MAX
#define INT_MAX ((int) (~(((unsigned) 1) << ((8*sizeof(int))-1))))
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif  /* __MACHDEFS_H */
/* RCSINFO $Id: util.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_UTIL_H
#define ILL_UTIL_H 


#ifdef _USRDLL

#ifdef QSLIB_EXPORTS
#define QSLIB_INTERFACE __declspec(dllexport)
#else
#define QSLIB_INTERFACE __declspec(dllimport)
#endif

#else 

#define QSLIB_INTERFACE extern 

#endif 

#ifdef WIN32 
#define strcasecmp(s1, s2) 	stricmp(s1, s2)
#define strncasecmp(s1, s2, n) 	strnicmp(s1, s2, n)
#endif 

/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  ILL_SWAP(a,b,t)                                                          */
/*    swaps a and b, using t as temporary space.  a, b, and t should all    */
/*    be the same type.                                                     */
/*                                                                          */
/*  ILL_OURABS(a)                                                            */
/*    returns the absolute value of a.                                      */
/*                                                                          */
/****************************************************************************/
typedef char ILLbool;
#define FALSE 0
#define TRUE  1

#define ILL_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define ILL_OURABS(a) (((a) >= 0) ? (a) : -(a))


/****************************************************************************/
typedef struct ILLrandstate {
    int a;
    int b;
    int arr[55];
} ILLrandstate;

/****************************************************************************/
/*                                                                          */
/*                             sortrus.c                                    */
/*                                                                          */
/****************************************************************************/
void
    ILLutil_int_array_quicksort (int *len, int n),
    ILLutil_int_perm_quicksort (int *perm, int *len, int n),
    ILLutil_double_perm_quicksort (int *perm, double *len, int n),
    ILLutil_str_perm_quicksort (int *perm, char**len, int n),
    ILLutil_rselect (int *arr, int l, int r, int m, double *coord,
        ILLrandstate *rstate);

char
    *ILLutil_linked_radixsort (char *data, char *datanext, char *dataval,
        int valsize);

/****************************************************************************/
/*                                                                          */
/*                             allocrus.c                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*  see allocrus.c                                                          */
/*                                                                          */
/****************************************************************************/

extern int ILLTRACE_MALLOC; 

#ifndef USEDMALLOC
#define ILL_UTIL_SAFE_MALLOC(nnum,type,varname)                             \
    (((ILLTRACE_MALLOC) ? printf("%s.%d: %s: ILL_UTIL_SAFE_MALLOC: %s = %d * %s\n", __FILE__, __LINE__, __DEV_FUNCTION__, #varname, nnum, #type) : 0), \
     (type *) ILLutil_allocrus (((size_t) (nnum)) * sizeof (type))) 
#else
#define ILL_UTIL_SAFE_MALLOC(nnum,type,varname)                             \
    (((ILLTRACE_MALLOC) ? printf("%s.%d: %s: ILL_UTIL_SAFE_MALLOC: %s = %d * %s\n", __FILE__, __LINE__, __DEV_FUNCTION__, #varname, nnum, #type) : 0), \
     (type *) malloc (((size_t) (nnum)) * sizeof (type)))
#endif

#define ILL_IFFREE(object,type) {                                           \
    if ((object)) {                                                        \
       ILLutil_freerus ((void *) (object));                                 \
       object = (type *) NULL;                                             \
    }}

#define ILL_PTRWORLD_ALLOC_ROUTINE(type, ptr_alloc_r, ptr_bulkalloc_r)        \
                                                                             \
static int ptr_bulkalloc_r (ILLptrworld *world, int nalloc)                   \
{                                                                            \
    ILLbigchunkptr *bp;                                                       \
    int i;                                                                   \
    int count = ILL_BIGCHUNK / sizeof ( type );                               \
    type *p;                                                                 \
                                                                             \
    while (nalloc > 0) {                                                     \
        bp = ILLutil_bigchunkalloc ();                                        \
        if (bp == (ILLbigchunkptr *) NULL) {                                  \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return 1;                                                        \
        }                                                                    \
        bp->next = world->chunklist ;                                        \
        world->chunklist = bp;                                               \
                                                                             \
        p = ( type * ) bp->this_one;                                         \
        for (i=count-2; i>=0; i--) {                                         \
            p[i].next = &p[i+1];                                             \
        }                                                                    \
        p[count - 1].next = (type *) world->freelist;                        \
        world->freelist = (void *) p;                                        \
        nalloc -= count;                                                     \
    }                                                                        \
    return 0;                                                                \
}                                                                            \
                                                                             \
static type *ptr_alloc_r (ILLptrworld *world)                                 \
{                                                                            \
    type *p;                                                                 \
                                                                             \
    if (world->freelist == (void *) NULL) {                                  \
        if (ptr_bulkalloc_r (world, 1)) {                                    \
            fprintf (stderr, "ptr alloc failed\n");                          \
            return ( type * ) NULL;                                          \
        }                                                                    \
    }                                                                        \
    p = (type *) world->freelist ;                                           \
    world->freelist = (void *) p->next;                                      \
                                                                             \
    return p;                                                                \
}

#define ILL_PTRWORLD_FREE_ROUTINE(type, ptr_free_r)                           \
                                                                             \
static void ptr_free_r (ILLptrworld *world, type *p)                          \
{                                                                            \
    p->next = (type *) world->freelist ;                                     \
    world->freelist = (void *) p;                                            \
}

#define ILL_PTRWORLD_LISTADD_ROUTINE(type, entrytype, ptr_listadd_r, ptr_alloc_r) \
                                                                             \
static int ptr_listadd_r (type **list, entrytype x, ILLptrworld *world)       \
{                                                                            \
    if (list != (type **) NULL) {                                            \
        type *p = ptr_alloc_r (world);                                       \
                                                                             \
        if (p == (type *) NULL) {                                            \
            fprintf (stderr, "ptr list add failed\n");                       \
            return 1;                                                        \
        }                                                                    \
        p->this = x;                                                         \
        p->next = *list;                                                     \
        *list = p;                                                           \
    }                                                                        \
    return 0;                                                                \
}

#define ILL_PTRWORLD_LISTFREE_ROUTINE(type, ptr_listfree_r, ptr_free_r)       \
                                                                             \
static void ptr_listfree_r (ILLptrworld *world, type *p)                      \
{                                                                            \
    type *next;                                                              \
                                                                             \
    while (p != (type *) NULL) {                                             \
        next = p->next;                                                      \
        ptr_free_r (world, p);                                               \
        p = next;                                                            \
    }                                                                        \
}

#define ILL_PTRWORLD_LEAKS_ROUTINE(type, ptr_leaks_r, field, fieldtype)       \
                                                                             \
static int ptr_leaks_r (ILLptrworld *world, int *total, int *onlist)          \
{                                                                            \
    int count = ILL_BIGCHUNK / sizeof ( type );                               \
    int duplicates = 0;                                                      \
    type * p;                                                                \
    ILLbigchunkptr *bp;                                                       \
                                                                             \
    *total = 0;                                                              \
    *onlist = 0;                                                             \
                                                                             \
    for (bp = world->chunklist ; bp; bp = bp->next)                          \
        (*total) += count;                                                   \
                                                                             \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        (*onlist)++;                                                         \
        p-> field = ( fieldtype ) 0;                                         \
    }                                                                        \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        if ((unsigned long) p-> field == (unsigned long) (size_t) 1)                           \
            duplicates++;                                                    \
        else                                                                 \
            p-> field = ( fieldtype ) (size_t) 1;                            \
    }                                                                        \
    if (duplicates) {                                                        \
        fprintf (stderr, "WARNING: %d duplicates on ptr free list \n",       \
                 duplicates);                                                \
    }                                                                        \
    return *total - *onlist;                                                 \
}

#define ILL_PTRWORLD_ROUTINES(type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r) \
ILL_PTRWORLD_ALLOC_ROUTINE (type, ptr_alloc_r, ptr_bulkalloc_r)               \
ILL_PTRWORLD_FREE_ROUTINE (type, ptr_free_r)

#define ILL_PTRWORLD_LIST_ROUTINES(type, entrytype, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r, ptr_listadd_r, ptr_listfree_r) \
ILL_PTRWORLD_ROUTINES (type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r)        \
ILL_PTRWORLD_LISTADD_ROUTINE (type, entrytype, ptr_listadd_r, ptr_alloc_r)    \
ILL_PTRWORLD_LISTFREE_ROUTINE (type, ptr_listfree_r, ptr_free_r)

#define ILL_BIGCHUNK ((int) ((1<<16) - sizeof (ILLbigchunkptr) - 16))

struct ILLbigchunk;

typedef struct ILLbigchunkptr {
    void                 *this_one;
    struct ILLbigchunk    *this_chunk;
    struct ILLbigchunkptr *next;
} ILLbigchunkptr;


typedef struct ILLptrworld {
    int refcount;
    void *freelist;
    ILLbigchunkptr *chunklist;
} ILLptrworld;



void
   *ILLutil_allocrus (size_t size),
   *ILLutil_reallocrus (void *ptr, size_t size),
    ILLutil_freerus (void *p),
    ILLutil_bigchunkfree (ILLbigchunkptr *bp),
    ILLptrworld_init (ILLptrworld *world),
    ILLptrworld_add (ILLptrworld *world),
    ILLptrworld_delete (ILLptrworld *world);

int
    ILLutil_reallocrus_scale (void **pptr, int *pnnum, int count, double scale,
        size_t size),
    ILLutil_reallocrus_count (void **pptr, int count, size_t size);

ILLbigchunkptr
    *ILLutil_bigchunkalloc (void);


/****************************************************************************/
/*                                                                          */
/*                             urandom.c                                    */
/*                                                                          */
/****************************************************************************/

/* since urandom's generator does everything modulo ILL_PRANDMAX, if two
 * seeds are congruent mod x and x|ILL_PRANDMAX, then the resulting numbers
 * will be congruent mod x.  One example was if ILL_PRANDMAX = 1000000000 and
 * urandom is used to generate a point set from a 1000x1000 grid, seeds
 * congruent mod 1000 generate the same point set.
 *
 * For this reason, we use 1000000007 (a prime)
 */
#define ILL_PRANDMAX 1000000007

void
   ILLutil_sprand (int seed, ILLrandstate *r);

int
   ILLutil_lprand (ILLrandstate *r);


/****************************************************************************/
/*                                                                          */
/*                             zeit.c                                       */
/*                                                                          */
/****************************************************************************/

typedef struct ILLutil_timer {
    double  szeit;
    double  cum_zeit;
    char    name[40];
    int     count;
} ILLutil_timer;


double
    ILLutil_zeit (void),
    ILLutil_real_zeit (void),
    ILLutil_stop_timer (ILLutil_timer *t, int printit),
    ILLutil_total_timer (ILLutil_timer *t, int printit);


void
    ILLutil_init_timer (ILLutil_timer *t, const char *name),
    ILLutil_start_timer (ILLutil_timer *t),
    ILLutil_suspend_timer (ILLutil_timer *t),
    ILLutil_resume_timer (ILLutil_timer *t);


/****************************************************************************/
/*                                                                          */
/*                             util.c                                       */
/*                                                                          */
/****************************************************************************/
#define ILL_UTIL_STR(new, str) \
    { new = ILLutil_str(str); \
      if (str != NULL) { ILL_CHECKnull(new, "out of memeory"); } }

extern char* ILLutil_str(const char* str) ; 
   /* allocates and returns a copy of s */

extern int ILLutil_array_index(char* list[], int n, const char* name); 
   /* returns index of name in list or -1  */

extern int ILLutil_index(const char* list[], const char* name); 
   /* returns index of name in list or -1  */

extern unsigned int ILLutil_nextprime (unsigned int x);

extern const char *ILLutil_strchr (const char *s, int c);

extern int ILLutil_strcasecmp(const char *s1, const char *s2); 
extern int ILLutil_strncasecmp(const char *s1, const char *s2, size_t n); 


extern int
    ILLutil_our_gcd (int a, int b),
    ILLutil_our_lcm (int a, int b),
    ILLutil_our_log2 (int a);

double
    ILLutil_our_floor (double x),
    ILLutil_our_ceil (double x),
    ILLutil_our_frac (double x),
    ILLutil_norm_sqr (double *v, int len);


/****************************************************************************/
/*                                                                          */
/*                             bgetopt.c                                    */
/*                                                                          */
/****************************************************************************/
int
    ILLutil_bix_getopt (int argc, char **argv, const char *def, int *p_optind,
        char **p_optarg);


#define ILL_BIX_GETOPT_UNKNOWN -3038


/****************************************************************************/
/*                                                                          */
/*                             dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct ILLdheap {
    double  *key;
    int     *entry;
    int     *loc;
    int     total_space;
    int     size;
} ILLdheap;

void
    ILLutil_dheap_free (ILLdheap *h),
    ILLutil_dheap_delete (ILLdheap *h, int i),
    ILLutil_dheap_changekey (ILLdheap *h, int i, double newkey),
    ILLutil_dheap_findmin (ILLdheap *h, int *i),
    ILLutil_dheap_deletemin (ILLdheap *h, int *i);

int
    ILLutil_dheap_init (ILLdheap *h, int k),
    ILLutil_dheap_resize (ILLdheap *h, int newsize),
    ILLutil_dheap_insert (ILLdheap *h, int i);


/****************************************************************************/
/*                                                                          */
/*                             priority.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct ILLpriority {
    ILLdheap   heap;
    union ILLpri_data {
        void *data;
        int  next;
    } *pri_info;
    int     space;
    int     freelist;
} ILLpriority;

void
    ILLutil_priority_free (ILLpriority *pri),
    ILLutil_priority_delete (ILLpriority *pri, int handle),
    ILLutil_priority_changekey (ILLpriority *pri, int handle, double newkey),
    ILLutil_priority_findmin (ILLpriority *pri, double *keyval, void **en),
    ILLutil_priority_deletemin (ILLpriority *pri, double *keyval, void **en);

int
    ILLutil_priority_init (ILLpriority *pri, int k),
    ILLutil_priority_insert (ILLpriority *pri, void *data, double keyval,
        int *handle);


#endif  /* ILL_UTIL_H */
/* RCSINFO $Id: except.h,v 1.3 2003/11/05 17:02:10 meven Exp $ */
#ifndef ILL_except
#define ILL_except

/* Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__'
   which contains the name of the function currently being defined.
#  define __DEV_FUNCTION__     __PRETTY_FUNCTION__
   This is broken in G++ before version 2.6.
   C9x has a similar variable called __func__, but prefer the GCC one since
   it demangles C++ function names.  */
# ifdef __GNUC__
#  if __GNUC__ > 2 || (__GNUC__ == 2 \
                       && __GNUC_MINOR__ >= (defined __cplusplus ? 6 : 4))
#   define __DEV_FUNCTION__    __PRETTY_FUNCTION__
#  else
#   define __DEV_FUNCTION__    ((__const char *) 0)
#  endif
# else
#  if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#   define __DEV_FUNCTION__    __func__
#  else
#   define __DEV_FUNCTION__    ((const char *) 0)
#  endif
# endif


/* put debugger breakpoint here */
extern void ILL_report(const char* msg,  
		      const char* fct, const char* file, unsigned int line, 
		      int with_source_info);
/* printed message looks as follows
 *
 * with_source_info == 0:   <msg> + "\n"
 *
 * with_source_info == 1:   if (fct != NULL) 
 *                              <msg> + " in function <fct>\n";
 *                          else 
 *                              <msg> + " in file <file> line <line>\n"; 
 */

#define ILL_GENERAL_ERROR -1
#define ILL_NO_MEMORY 2
#define ILL_NULL_PTR  3

#define ILL_REPORT(msg,with) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, with)
#ifdef NDEBUG 
#define ILL_REPRT(msg) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, 0)
#else 
#define ILL_REPRT(msg) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, 1)
#endif 

#define ILL_RESULT(expr, msg) \
{								\
	if (TRACE > 0) { ILL_RETURN(expr, msg); } \
    return expr; \
}

#define ILL_RETURN_PTR(ptr, msg) \
    { void *ILL_RETURN_p = ptr; \
      if (ILL_RETURN_p == NULL) {  \
        if (TRACE > 0) ILL_REPRT(msg); \
      } \
      return ILL_RETURN_p;  \
    }

#ifdef NDEBUG 
#define ILL_RETURN(expr, msg)           \
{										\
	if (expr != 0) {					\
			if (TRACE > 0) ILL_REPRT(msg); \
    }									\
    return expr;						\
}

#else 
#define ILL_RETURN(expr, msg)           \
	{									\
		if (expr != 0) {                \
			ILL_REPRT(msg);             \
		}                               \
		ILL_IFTRACE("%s: returning %d\n", __DEV_FUNCTION__, expr); \
		return expr;					\
	}
#endif

#define ILL_CHECKnull(expr, msg) \
    { if ((expr) == NULL)  {  \
         ILL_REPRT(msg); \
         rval = ILL_NULL_PTR; \
         goto CLEANUP;  \
      } }

#define ILL_FAILtrue(expr, msg) \
    { if (expr)  {  \
         ILL_REPRT(msg); \
         rval = ILL_GENERAL_ERROR; \
         goto CLEANUP;  \
      } }

#define ILL_FAILtrue_no_rval(expr, msg) \
    { if (expr)  {  \
         ILL_REPRT(msg); \
         goto CLEANUP;  \
      } }


#define ILL_FAILfalse(expr, msg)  ILL_FAILtrue(!(expr), msg) 
#define ILL_FAILfalse_no_rval(expr, msg)  ILL_FAILtrue_no_rval(!(expr), msg) 

#define ILL_ERROR(rval, msg)      {									\
									fprintf(stderr, "%s\n", msg);	\
									rval = 1; goto CLEANUP;			\
								  }
#define ILL_CLEANUP_IF(rval)      { if ((rval) != 0) { goto CLEANUP; } }
#define ILL_CLEANUP               goto CLEANUP

#define ILL_SAFE_MALLOC(lhs, n, type) \
    { lhs = ILL_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         ILL_REPRT("Out of memory"); \
         rval = ILL_NO_MEMORY; \
         goto CLEANUP;  \
      }}

#define ILL_SAFE_MALLOC_no_rval(lhs, n, type) \
    { lhs = ILL_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         ILL_REPRT("Out of memory"); \
         goto CLEANUP;  \
      }}


#define ILL_NEW(ptr, type) ILL_SAFE_MALLOC(ptr, 1, type)
#define ILL_NEW_no_rval(ptr, type) ILL_SAFE_MALLOC_no_rval(ptr, 1, type)

#endif
/* RCSINFO $Id: symtab.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_SYMTAB_H
#define ILL_SYMTAB_H

#define ILL_namebufsize 2048

typedef struct ILLsymbolent {
    int symbol;
    int index;
    int next;
} ILLsymbolent;

typedef struct ILLsymboltab {
    int *hashtable;
    ILLsymbolent *nametable;
    char *namelist;
    int tablesize;
    int strsize;
    int hashspace;
    int namespace;
    int strspace;
    int freedchars; 
    int the_hash; 
    int the_index; 
    int the_prev_index; 
    int index_ok;
} ILLsymboltab;
/* 
 * hashtable[stringhash(entry) % hashspace] either NO_INDEX or some hash number
 * nametable[hash number] = { next:    another index for nametable 
 *                            symbol:  index into namelist where string resides
 *                          }
 * tablesize:  number of entries   (tablesize <= namespace)
 * namespace:  length of nametable and indexlist
 * hashspace:  length of hashtable   nextprime(namespace)
 * strsize:    number of chars used in namelist 
 * strspace:   length of namelist
 * indexlist:  LP col/row indices for the table entries
 * indexlist_ok:  1 if column indices in indexlist are up-to-date, 0 otherwise
 *
 * Deletion of entries affects their ordering in symboltab->nametable. 
 * Strings may move around within symboltab->namelist. 
 */

    
#define ILL_SYM_NOINDEX (-1)
extern void
    ILLsymboltab_init (ILLsymboltab *h),
    ILLsymboltab_free (ILLsymboltab *h),
    ILLsymboltab_size (const ILLsymboltab *h, int *p_size),
    ILLsymboltab_prt (FILE *fd, ILLsymboltab *h) ; 

extern int
    ILLsymboltab_create (ILLsymboltab *h, int init_size),
    ILLsymboltab_copy (ILLsymboltab *src, ILLsymboltab *dst),
    ILLsymboltab_register (ILLsymboltab *h, const char *s, int itemindex,
        int *p_index, int *p_existed),
    ILLsymboltab_lookup (ILLsymboltab *h, const char *s, int *p_index),
    ILLsymboltab_index_ok (ILLsymboltab *h),
    ILLsymboltab_index_reset (ILLsymboltab *h, int icount, char **names),
    ILLsymboltab_getindex (ILLsymboltab *h, const char *name, int *hindex),
    ILLsymboltab_contains(ILLsymboltab *h, const char *s), 
    ILLsymboltab_delete(ILLsymboltab *h, const char *s),
    ILLsymboltab_uname(ILLsymboltab *h, 
		       char name[ILL_namebufsize],
		       const char *try_prefix1, const char* try_prefix2);

extern void ILLsymboltab_unique_name(ILLsymboltab *tab, int i, const char* pref,  
		       char uname[ILL_namebufsize] ) ;

extern const char* ILLsymboltab_get(const ILLsymboltab *tab, int i); 
extern int ILLsymboltab_rename(ILLsymboltab *h, int i, const char *new_name); 

#endif /* __SYMTAB_H */
/* RCSINFO $Id: names.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_NAMES_H
#define ILL_NAMES_H

extern void ILLfree_names (char **names, int count); 
extern int ILLgenerate_names(char prefix, int nnames, char ***names); 

#endif 
