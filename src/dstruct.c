/* RCS_INFO = "$RCSfile: dstruct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0; 

#include "iqsutil.h"
#include "dstruct.h"
#include "qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/****************************************************************************/
/*                                                                          */
/*                            svector                                       */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

void ILLsvector_init (svector *s)
{
    s->nzcnt = 0;
    s->indx = (int *) NULL;
    s->coef = (double *) NULL;
}

void ILLsvector_free (svector *s)
{
    s->nzcnt = 0;
    ILL_IFFREE (s->indx, int);
    ILL_IFFREE (s->coef, double);
}

int ILLsvector_alloc (svector *s, int nzcnt)
{
    int rval = 0;
    
    s->nzcnt = nzcnt;
    if (nzcnt == 0) {
        s->indx = (int *) NULL;
        s->coef = (double *) NULL;
    } else {
        ILL_SAFE_MALLOC (s->indx, nzcnt, int);
        ILL_SAFE_MALLOC (s->coef, nzcnt, double);
    }
    return 0;
 CLEANUP:
    ILL_IFFREE (s->indx, int);
    ILL_IFFREE (s->coef, double);
    ILL_RETURN (rval, "ILLsvector_alloc");
}

int ILLsvector_copy (const svector *s_in, svector *s_out)
{
    int i;
    int nzcnt = s_in->nzcnt;
    int rval = 0;
    
    rval = ILLsvector_alloc (s_out, nzcnt);
    ILL_CLEANUP_IF (rval);
    for (i=0; i<nzcnt; i++) {
        s_out->indx[i] = s_in->indx[i];
        s_out->coef[i] = s_in->coef[i];
    }

 CLEANUP:
    ILL_RETURN (rval, "ILLsvector_copy");
}

/****************************************************************************/
/*                                                                          */
/*                            heap                                          */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

#define DEBUG_HEAP 0

#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)

static int
   siftup (heap *h, int hloc, int ix),
   siftdown (heap *h, int hloc, int ix),
   maxchild (heap *h, int hloc);

static int siftup (heap *h, int hloc, int ix)
{
   int  i = hloc;
   int  p = HEAP_UP (i);
   double  val = h->key[ix];

   while (i > 0 && val > h->key[h->entry[p]]){
      h->entry[i] = h->entry[p];
      h->loc[h->entry[i]] = i;
      i = p;
      p = HEAP_UP (p);
   }
   h->entry[i] = ix;
   h->loc[ix] = i;
   return i;
}

static int siftdown (heap *h, int hloc, int ix)
{
   int  i = hloc;
   int  c = maxchild (h, i);
   double  val = h->key[ix];

   while (c != -1 && val < h->key[h->entry[c]]){
      h->entry[i] = h->entry[c];
      h->loc[h->entry[i]] = i;
      i = c;
      c = maxchild (h, c);
   }
   h->entry[i] = ix;
   h->loc[ix] = i;
   return i;
}

static int maxchild (heap *h, int hloc)
{
   int  i;
   int  mc = -1;
   int  hmin = HEAP_D * hloc + 1;
   int  hmax = HEAP_D * hloc + HEAP_D;
   double  mv = -QS_MAXDOUBLE;

   for (i = hmin; i <= hmax && i < h->size; i++)
      if (mv < h->key[h->entry[i]]) { mv = h->key[h->entry[i]]; mc = i; }
   return mc;
}

#if DEBUG_HEAP > 0

static void printheap (heap *h)
{
   int i;

   printf ("entry: ");
   for (i=0; i<h->size; i++) printf ("%d ", h->entry[i]);
   printf ("\n loc: ");
   for (i=0; i<h->maxsize; i++) printf ("%d ", h->loc[i]);
   printf ("\n key: ");
   for (i=0; i<h->maxsize; i++) printf ("%.2f ", h->key[i]);
   printf ("\n key(sorted): ");
   for (i=0; i<h->size; i++) printf ("%.2f ", h->key[h->entry[i]]);
   printf ("\n");
}

static void heapcheck (heap *h)
{ 
   int i, tcnt = 0;

   for (i=0; i<h->maxsize; i++){
      if (h->loc[i] < -1) printf("error in heap\n");
      else if (h->loc[i] > -1) tcnt++;
   }
   if (tcnt != h->size) printf("error 3 in heap\n");

   for (i=0; i<h->size; i++){
      if (h->loc[h->entry[i]] != i) printf ("error 1 in heap\n");
      if (h->key[h->entry[i]] == 0.0) printf ("error 2 in heap\n");
      if (h->key[h->entry[i]] > h->key[h->entry[HEAP_UP(i)]])
                               printf ("error 4 in heap\n");
   }
}

#endif

void ILLheap_insert (heap *h, int ix)
{
   int  i = h->size;
    
   i = siftup (h, i, ix);
   h->size ++;

#if DEBUG_HEAP > 0
   heapcheck (h);
#endif
#if DEBUG_HEAP > 1
   printheap (h);
#endif
}

void ILLheap_modify (heap *h, int ix)
{
   int  i  = h->loc[ix];
   int  pi = i;

   if (h->loc[ix] == -1) return;
   i = siftup (h, i, ix);
   if (pi == i)
      i = siftdown (h, i, ix);

#if DEBUG_HEAP > 0
   heapcheck (h);
#endif
#if DEBUG_HEAP > 1
   printheap (h);
#endif
}

void ILLheap_delete (heap *h, int ix)
{
   int  i   = h->loc[ix];
   int  pi  = i;
   int  nix = h->entry[h->size-1];

   h->loc[ix] = -1;
   h->size--;
   if (nix == ix){
#if DEBUG_HEAP > 0
   heapcheck (h);
#endif
#if DEBUG_HEAP > 1
   printheap (h);
#endif
      return;
   }

   h->entry[i] = nix;
   h->loc[nix] = i;

   i = siftup (h, i, nix);
   if (pi == i)
      siftdown (h, i, nix);

#if DEBUG_HEAP > 0
   heapcheck (h);
#endif
#if DEBUG_HEAP > 1
   printheap (h);
#endif
}

int ILLheap_findmin (heap *h)
{
   if (h->hexist == 0 || h->size <= 0)
      return -1;
   return h->entry[0];
}

void ILLheap_init (heap *h)
{
   h->entry  = NULL;
   h->loc    = NULL;
   h->key    = NULL;
   h->hexist = 0;
}

int ILLheap_build (heap *h, int nelems, double *key)
{
   int  rval = 0;
   int  i, n = 0;

   h->hexist   = 1;
   h->size     = 0;
   h->maxsize  = nelems;
   h->key      = key;
   ILL_SAFE_MALLOC (h->entry, nelems, int);
   ILL_SAFE_MALLOC (h->loc, nelems, int);

   for (i=0; i<nelems; i++){
      if (key[i] > 0.0) {
         h->entry[n] = i; h->loc[i] = n; n++;
      }
      else
         h->loc[i] = -1;
   }
   h->size = n;
   for (i=n-1; i>= 0; i--)
      siftdown (h, i, h->entry[i]);

#if DEBUG_HEAP > 0
   heapcheck (h);
#endif
#if DEBUG_HEAP > 1
   printheap (h);
#endif

 CLEANUP:
   if (rval) ILLheap_free (h);
   ILL_RETURN (rval, "ILLheap_init");
}

void ILLheap_free (heap *h)
{
   if (h->hexist){
      ILL_IFFREE (h->entry, int);
      ILL_IFFREE (h->loc, int);
      h->hexist = 0;
      h->maxsize = 0;
      h->size = 0;
   }
}


/****************************************************************************/
/*                                                                          */
/*                          matrix                                          */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

void ILLmatrix_init (ILLmatrix *A) 
{
    if (A) {
        A->matval     = (double *) NULL;
        A->matcnt     = (int *) NULL;
        A->matbeg     = (int *) NULL;
        A->matind     = (int *) NULL;
        A->matcols    = 0;
        A->matcolsize = 0;
        A->matrows    = 0;
        A->matsize    = 0;
        A->matfree    = 0;
    }
}

void ILLmatrix_free (ILLmatrix *A)
{
    if (A) {
        ILL_IFFREE (A->matval, double);
        ILL_IFFREE (A->matcnt, int);
        ILL_IFFREE (A->matbeg, int);
        ILL_IFFREE (A->matind, int);
        ILLmatrix_init (A);
    }
}

void ILLmatrix_prt (FILE *fd, ILLmatrix *A)
{
    int j, k; 
    if (A == NULL) {
	fprintf(fd, "Matrix %p: empty\n", (void*) A);  
    } else {
	fprintf(fd, "Matrix %p: nrows = %d ncols = %d\n",
		(void *) A, A->matrows, A->matcols); 
	for (j = 0; j < A->matcols; j++) {
	    fprintf(fd, "col %d: ", j); 
	    for (k = A->matbeg[j]; k < A->matbeg[j] + A->matcnt[j]; k++) {
		fprintf(fd, "row %d=%.3f ", A->matind[k], A->matval[k]); 
	    } 
	    fprintf(fd, "\n"); 
        }
    }
}
