#ifndef _calculus_h
#define _calculus_h

#include <math.h>

#include "memoryHandling.h"
#include "arrayOperations.h"

double complex Csimps1D(int, Carray, double);

double Rsimps1D(int, Rarray, double);

double complex Csimps2D(int, int, Carray, double, double);

double Rsimps2D(int, int, Rarray, double, double);

double Rsimps2D_sphere(int, int, Rarray, Rarray, double);

void renormalize(int, int, Carray, double, double, double);

void renormalize_spheric(EqDataPkg EQ, Carray S);

void renormalizeReal(int, int, Rarray, double, double, double);

void derivative_periodic(int, Carray, double, Carray);






#endif
