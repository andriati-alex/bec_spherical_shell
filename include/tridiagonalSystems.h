#ifndef _tridiagonalSystems_h
#define _tridiagonalSystems_h

#include "memoryHandling.h"
#include "arrayOperations.h"

void tridiag(int, Carray, Carray, Carray, Carray, Carray);
void lu_decomposition(int, Carray, Carray, Carray, Carray, Carray);
void tridiag_lu(int, Carray, Carray, Carray, Carray, Carray, Carray);


#endif
