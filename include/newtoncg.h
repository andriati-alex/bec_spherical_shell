#ifndef _newtoncg_h
#define _newtoncg_h

#include <mkl.h>
#include <mkl_dfti.h>
#include "inout.h"
#include "calculus.h"
#include "observables.h"
#include "arrayOperations.h"

void stationaryNewton(EqDataPkg EQ, Carray Sa, Carray Sb,
        double err_tol, int iter_tol, double mu_a, double mu_b);

void stationaryFixedNorm(EqDataPkg EQ, Carray Sa, Carray Sb,
        double err_tol, int iter_tol);

#endif
