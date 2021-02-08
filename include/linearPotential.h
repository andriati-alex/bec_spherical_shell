
#ifndef _linearPotential_h
#define _linearPotential_h

#include <string.h>
#include <stdio.h>
#include "arrayOperations.h"

void harmonic(int,int,Rarray,Rarray,Rarray,double,double);

void quartic(int,int,Rarray,Rarray,Rarray,double,double);

void quarticQuadratic(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wy, double a, double b);

void HarmonicMexicanHat(int,int,Rarray,Rarray,Rarray,double,double,double,
        double,double);

void QuarticMexicanHat(int,int,Rarray,Rarray,Rarray,double,double,double,
        double,double);

void GetPotential(char [],int,int,Rarray,Rarray,Rarray,double []);

#endif
