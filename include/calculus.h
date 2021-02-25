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

double complex Csimps2D_sphere(int, int, Rarray, Carray, double);

void renormalize(int, int, Carray, double, double, double);

void renormalize_spheric(EqDataPkg EQ, Carray S);

void renormalizeReal(int, int, Rarray, double, double, double);

void derivative_1dperiodic(int, Carray, double, Carray);

void twice_derivative_1dperiodic(int, Carray, double, Carray);

void sph_phi_derivative(int, int, Carray, double, Carray);

void sph_phi_twice_derivative(int, int, Carray, double, Carray);

void sph_theta_derivative(int, int, Carray, double, Carray);

void sph_theta_twice_derivative(int, int, Carray, double, Carray);






#endif
