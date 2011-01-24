/* $RCSfile: stddefs.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef ILL_STDDEFS_H
#define ILL_STDDEFS_H

#define SWAP(x,y,temp) {temp = x; x = y; y = temp;}

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define ABS(x) ((x) >= 0 ? (x) : -(x))

#endif /* ILL_STDDEFS_H */
