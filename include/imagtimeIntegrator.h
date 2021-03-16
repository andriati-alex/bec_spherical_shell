#ifndef _imagtimeIntegrator_h
#define _imagtimeIntegrator_h

#include "tridiagonalSystems.h"
#include "arrayOperations.h"
#include "observables.h"
#include "inout.h"

int splitstep_spherical_shell(EqDataPkg, Carray, Carray);

int splitstep_spherical_shell_single(EqDataPkg, Carray);

int splitstep_theta_sphere(EqDataPkg, Carray, Carray, int, int);

#endif
