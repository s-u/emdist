/* RCSINFO $Id: names.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#include "names.h" 
#include "util.h" 
#include "except.h" 
#include "symtab.h" 
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


void ILLfree_names (char **names, int count)
{
    int i;

    if (names) {
        for (i = 0; i < count; i++) {
            ILL_IFFREE (names[i], char);
        }
        ILL_IFFREE (names, char *);
    }
}

int ILLgenerate_names(char prefix, int nnames, char ***names) 
{
    int rval = 0;
    int i, len;
    char *buf = (char *) NULL;

    *names = (char **) NULL;
    if (nnames == 0) goto CLEANUP;

    ILL_SAFE_MALLOC (buf, ILL_namebufsize, char);
    ILL_SAFE_MALLOC (*names, nnames, char *);
    for (i = 0; i < nnames; i++) {
        (*names)[i] = (char *) NULL;
    }

    for (i = 0; i < nnames; i++) {
        sprintf (buf, "%c%d", prefix, i);
        len = strlen (buf) + 1;
        ILL_SAFE_MALLOC ((*names)[i], len, char);
        strcpy ((*names)[i], buf); 
    }

CLEANUP:

    if (rval) {
        if (*names) {
            ILLfree_names (*names, nnames);
            *names = (char **) NULL;
        }
    }

    ILL_IFFREE (buf, char);
    if (rval) {
         fprintf (stderr, "ILLsymboltab_generate_names failed\n");
    } 
    return rval;
}
