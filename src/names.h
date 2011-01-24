/* RCSINFO $Id: names.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_NAMES_H
#define ILL_NAMES_H

extern void ILLfree_names (char **names, int count); 
extern int ILLgenerate_names(char prefix, int nnames, char ***names); 

#endif 
