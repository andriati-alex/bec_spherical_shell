#ifndef _newtoncg_h
#define _newtoncg_h

#include <mkl.h>
#include <mkl_dfti.h>
#include "inout.h"
#include "calculus.h"
#include "observables.h"
#include "arrayOperations.h"

void stationaryNewton(EqDataPkg, Carray, Carray, double, int);

#endif
