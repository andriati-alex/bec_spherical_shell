#ifndef _imagtimeIntegrator_h
#define _imagtimeIntegrator_h

#include "tridiagonalSystems.h"
#include "arrayOperations.h"
#include "observables.h"
#include "inout.h"

void _get_psi_theta(int, int, int, Carray, Carray);

void _set_psi_theta(int, int, int, Carray, Carray);

int splitstep_spherical_shell(EqDataPkg, Carray, Carray);

int splitstep_spherical_shell_single(EqDataPkg, Carray);

int splitstep_theta_sphere(EqDataPkg, Carray, Carray, int, int);

int splitstep_theta_sphere_equals(EqDataPkg, Carray, Carray, int, int);

#endif
