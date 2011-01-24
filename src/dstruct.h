/* RCSINFO $Id: dstruct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                           svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef __SVECTOR_H
#define __SVECTOR_H

typedef struct svector {
    int     nzcnt;
    int    *indx;
    double *coef;
} svector;

void
    ILLsvector_init (svector *s),
    ILLsvector_free (svector *s);

int
    ILLsvector_alloc (svector *s, int nzcnt),
    ILLsvector_copy (const svector *s_in, svector *s_out);

#endif /* __SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef __HEAP_H
#define __HEAP_H

typedef struct{
   int  *entry;
   int  *loc;
   int  hexist;
   int  maxsize;
   int  size;
   double *key;
} heap;

void
    ILLheap_insert (heap *h, int ix),
    ILLheap_modify (heap *h, int ix),
    ILLheap_delete (heap *h, int ix),
    ILLheap_init (heap *h),
    ILLheap_free (heap *h);

int
    ILLheap_findmin (heap *h),
    ILLheap_build (heap *h, int nelems, double *key);

#endif /* __HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef __MATRIX_H
#define __MATRIX_H

typedef struct ILLmatrix {
    double *matval;             /* The coefficients.                       */
    int    *matcnt;             /* Number of coefs in each col.            */
    int    *matind;             /* The row indices of the coefs.           */
    int    *matbeg;             /* The start of each col.                  */
    int     matcols;            /* Number of columns.                      */
    int     matrows;
    int     matcolsize;         /* Length of matbeg and matcnt.            */
    int     matsize;            /* Length of matind and matval.            */
    int     matfree;            /* Free space at end of matind.            */
                                /* Note: free elements marked by -1 in     */
                                /* matind; we keep at least 1 free at end. */
} ILLmatrix;

void ILLmatrix_init (ILLmatrix *A);
void ILLmatrix_free (ILLmatrix *A);
void ILLmatrix_prt (FILE *fd, ILLmatrix *A);

#endif /* __MATRIX_H */
