#ifndef _arrayOperations_h
#define _arrayOperations_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <math.h>
#include "dataStructures.h"


void carrRandom(int n, int seed, Carray v);

void carrFill(int n, double complex z, Carray v);
void rarrFill(int n, double x, Rarray v);
// Fill the elements of an array with number z(complex) or x(real)



void rarrFillInc(int n, double x0, double dx, Rarray v);
// Fill starting from x0 increasing dx elementwise through the vector
// x[0] = x0, x[1] = x0 + dx, ... , x[n-1] = x0 + (n-1)*dx



void carrCopy(int n, Carray from, Carray to);
void carrCopy_from_real(int n, Rarray from, Carray to);
void rarrCopy(int n, Rarray from, Rarray to);
// Copy elements "from" array "to" array

void carrCutBounds(int n, int chunk, Carray inp, Carray out);



void MKL2Carray(int n, CMKLarray a, Carray b);
void Carray2MKL(int n, Carray b, CMKLarray a);
// Convert array datatypes



/*          ***********************************************

                     BASIC OPERATIONS ELEMENT-WISE

            ***********************************************          */



void carrRealPart(int n, Carray v, Rarray vreal);

void carrImagPart(int n, Carray v, Rarray vimag);

void carrConj(int n, Carray v, Carray v_conj);
// Element-wise conjugation

void carrAdd(int n, Carray v1, Carray v2, Carray v);
void rarrAdd(int n, Rarray v1, Rarray v2, Rarray v);
// Element-wise addition

void carrSub(int n, Carray v1, Carray v2, Carray v);
void rarrSub(int n, Rarray v1, Rarray v2, Rarray v);
// Element-wise subtraction. convention v1 - v2

void carrConjMultiply(int n, Carray v1, Carray v2, Carray v);
void carrMultiply(int n, Carray v1, Carray v2, Carray v);
void rarrMultiply(int n, Rarray v1, Rarray v2, Rarray v);
// Element-wise product v1 x v2

void carrScalarMultiply(int n, Carray v, double complex z, Carray ans);
void rarrScalarMultiply(int n, Rarray v, double z, Rarray ans);
// multiply all elements by a scalar

void carrScalarAdd(int n, Carray v, double complex z, Carray ans);
void rarrScalarAdd(int n, Rarray v, double z, Rarray ans);
// Element-wise addition by a scalar

void carrDiv(int n, Carray v1, Carray v2, Carray v);
void rarrDiv(int n, Rarray v1, Rarray v2, Rarray v);
// Element-Wise division v1 / v2

void carrUpdate(int n, Carray v1, double complex z, Carray v2, Carray v);
void rcarrUpdate(int n, Carray v1, double complex z, Rarray v2, Carray v);
void rarrUpdate(int n, Rarray v1, double z, Rarray v2, Rarray v);
// Update-like operation v[i] = v1[i] + z * v2[i];

void carrAbs(int n, Carray v, Rarray vabs);
void rarrAbs(int n, Rarray v, Rarray vabs);
// Take absolute(complex) value for each component

void carrAbs2(int n, Carray v, Rarray vabs);
void rarrAbs2(int n, Rarray v, Rarray vabs);
// Take absolute(complex) squared value for each component

void renormalizeVector(int n, Carray v, double norm);
// Renormalize vector such that || v || = norm



/*          ***********************************************

                               FUNCTIONALS

            ***********************************************          */



// Default scalar product. Convention <v1*, v2>
double complex carrDot(int n, Carray v1, Carray v2);

// Product and sum elements. Convention <v1, v2>
double complex unconj_carrDot(int n, Carray v1, Carray v2);

double rarrDot(int n, Rarray v1, Rarray v2);

// Vector modulus
double carrMod(int n, Carray v);
double rarrMod(int n, Rarray v);

// Vector Modulus squared
double carrMod2(int n, Carray v);
double rarrMod2(int n, Rarray v);

// Sum all elements of a vector
double complex carrReduction(int n, Carray v);
double rarrReduction(int n, Rarray v);



/*          ***********************************************

                     FUNCTION COMPUTED ELEMENT-WISE

            ***********************************************          */



// take exponential of each component of z * v.
void carrExp(int n, double complex z, Carray v, Carray ans);
void rcarrExp(int n, double complex z, Rarray v, Carray ans);

#endif
