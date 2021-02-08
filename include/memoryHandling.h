#ifndef _memoryHandling_h
#define _memoryHandling_h

#include <stdio.h>
#include <stdlib.h>
#include "dataStructures.h"
#include "linearPotential.h"



Rarray rarrDef(int n); // Allocate real vector of n elements

Carray carrDef(int n); // Allocate complex vector of n elements

CMKLarray cmklDef(int n); // Allocate MKL's complex vector

Rmatrix rmatDef(int m, int n); // real matrix with m rows and n columns

Cmatrix cmatDef(int m, int n); // complex matrix with m rows and n columns



CCScmat ccscmatDef(int n, int max_nonzeros);
// Allocate CCS matrix structure with  n  rows
// with the maximum number of nonzero elements
// max_nonzeros

CCSrmat ccsrmatDef(int n, int max_nonzeros);



EqDataPkg equation_structure(int, int, double, int, double, double,
        double, double, double, double, double);
// Pack all equation parameters



void rmatFree(int m, Rmatrix M); // Release real matrix with m rows

void cmatFree(int m, Cmatrix M); // Release complex matrix with m rows

void ccscmatFree(CCScmat M); // Relase CCS matrix of complex numbers

void ccsrmatFree(CCSrmat M); // Relase CCS matrix of real numbers

void ReleaseEqDataPkg(EqDataPkg);

#endif
