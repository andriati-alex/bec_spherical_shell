#ifndef _dataStructures_h
#define _dataStructures_h

#include <complex.h>
#include <mkl.h>



/* DATATYPES OF THE PACKAGE
 * *********************************************************
 *
 * Description
 * ===========
 *
 * Header with all datatypes used on this package
 *
 *
 *
 * Few words on Compressed-Column-Storage(CCS) of a Matrix
 * =======================================================
 *
 * CCS is a method to save memory when storing sparse  matrices
 * It has a good efficiency if the number of zeros in each  row
 * are almost the same,  and maximum efficiency when  they  are
 * identical. It works as follows, given a matrix expressed  in
 * the conventional form with 'M' rows and 'N' columns one need
 * to  shift  all nonzero elements to the left.  The  remaining
 * matrix will have a column 'm'  whose any other column  k > m
 * will contain only zeros, these ones are then cutted off. The
 * final step is to concatenate  the  remaining  columns  in  a
 * vector with the values and create another vector of integers
 * with the original column number. The variable  'm'  indicate
 * the maximum number of nonzero elements in  a  same  row. The
 * leading vectors with values and former columns  numbers  are
 * storaged in a column major format where given v[i] the former
 * row number were (i % M)
 *
 *
 *
 * The Equation Package structure
 * ==============================================
 *
 * It is created a structure to hold all equation's coefficients
 * and external potential (usually trap).
 *
 * *************************************************************/





#define PI 3.141592653589793



/* Shortcut for datatypes and pointers to them
 * =========================================== */

typedef double complex doublec;    // default complex number
typedef MKL_Complex16 * CMKLarray; // mkl complex vector
typedef double complex * Carray;   // default complex vector
typedef double * Rarray;           // real vector

typedef double ** Rmatrix;
typedef double complex ** Cmatrix;










/* Compressed Column storage of an n x n Real Matrix
 * ================================================= */

struct _CCScmat
{
    int m;     // max number of non-zero elements in a same row
	int * col;  // Column index of elemetns.
	Carray vec; // column oriented vector.
};

typedef struct _CCScmat * CCScmat;


/* Compressed Column storage of an n x n Real Matrix
 * ================================================= */

struct _CCSrmat
{
    // same fields from complex case above
    int m;
	int * col;
	Rarray vec;
};

typedef struct _CCSrmat * CCSrmat;










/* Structure to store all relevant equation information */

struct _EquationDataPkg
{
    int
        nt,         // total number of time steps to propagate
        nphi,       // number of grid points in phi direction
        ntheta;     // number of grid points in theta direction
    double
        dphi,       // grid spacing in phi-direction
        dtheta,     // grid spacing in theta-direction
        dt,         // time step spacing
        nabla_coef, // real number multiplying laplacian
        omega,      // rotation frequency
        frac_a,     // population fraction of species type A
        frac_b,     // population fraction of species type B
        ga,         // contact interaction strength
        gb,
        gab;
    Rarray
        phi,        // phi grid points
        theta;      // theta grid points
};

typedef struct _EquationDataPkg * EqDataPkg;

#endif
