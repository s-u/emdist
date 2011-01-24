/* RCSINFO $Id: util.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_UTIL_H
#define ILL_UTIL_H 

#include "machdefs.h" 

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
