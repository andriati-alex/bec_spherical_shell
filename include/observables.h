#ifndef _observables_h
#define _observables_h

#include "calculus.h"

double angular_momentum_lz(int, int, double, Rarray, Carray);

double functionals_single(EqDataPkg, Carray, double *, double *);

double functionals(EqDataPkg, Carray, Carray, double *, double *, double *);

double functionals_theta(
        EqDataPkg, Carray, Carray, double *, double *, double *, int, int
);

double density_overlap(EqDataPkg, Rarray, Rarray);
double theta_density_overlap(EqDataPkg, Rarray, Rarray);

double avg_residue(EqDataPkg, Carray, Carray, double, double);

#endif
