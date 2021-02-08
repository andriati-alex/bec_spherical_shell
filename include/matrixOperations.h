#ifndef _matrixOperations_h
#define _matrixOperations_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "memoryHandling.h"
#include "arrayOperations.h"



/*          ***********************************************          */
/*                   SETUP VALUES IN MATRIX ENTRIES                  */
/*          ***********************************************          */



void cmatFill(int m, int n, double complex z, Cmatrix M);
/* Fill all matrix with a constant value
 * *************************************
 *
 *  m the number of lines
 *  n the number of columns
 *  z the value to fill the entire matrix
 *
 * *************************************/



void cmatFillDK(int n, int k, Carray z, Cmatrix M);
/* Fill k Diagonal of square matrix of n by n elements with an array z
 * *******************************************************************
 *
 * z must contain n - k elements
 * k = 0 is the main diagonal
 * k > 0 upper diagonals
 * k < 0 lower diagonals
 *
 * *******************************************************************/



void cmatFillTri(int n, Carray upper, Carray mid, Carray lower, Cmatrix M);
/*   Fill a matrix with just the tridiagonals entries, rest with zeros   */



void RowMajor(int m, int n, Cmatrix M, Carray v);
/* Store Matrix M(m x n) in a vector v using Row Major scheme */



void setValueCCScmat(int n, int i, int j, int col, doublec z, CCScmat M);
/*        Set a value in the i row and original column number col      */

void setValueCCSrmat(int n, int i, int j, int col, double x, CCSrmat M);
/*        Set a value in the i row and original column number col      */





/*          **********************************************          */
/*          MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION          */
/*          **********************************************          */





void cmatvec(int m, int n, Cmatrix M, Carray v, Carray ans);
/* General Matrix Vector multiplication: M . v = ans
 * *************************************************
 *
 *  m number of lines of M
 *  n number of columns of M and components of v
 *
 * *************************************************/



void cmatmat(int m, int n, int l, Cmatrix M, Cmatrix A, Cmatrix ans);
/* General Matrix Matrix multiplication: M . A = ans
 * *************************************************
 *
 *  M has m(lines) by n(columns)
 *  A has n(lines) by l(columns)
 *  ans has m(lines) by l(columns)
 *
 * *************************************************/



void CCScmatvec(int n, Carray vals, int * cols, int m, Carray vec, Carray ans);
/* Matrix(in CCS format) vector multiplication
 * 
 * Given CCSmat A the arguments taken are
 *
 * vals = A->vec
 * cols = A->col
 * m    = A->m
 *
 * *******************************************/



void CCSrmatvec(int n, Rarray vals, int * cols, int m, Rarray vec, Rarray ans);


#endif
