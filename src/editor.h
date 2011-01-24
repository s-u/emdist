/* RCSINFO $Id: editor.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef EDITOR_H
#define EDITOR_H

extern void ILLeditor_init(void) ; 
extern void ILLeditor(QSdata* p); 
extern int  ILLeditor_solve(QSdata* p, int salgo) ; 

#endif
