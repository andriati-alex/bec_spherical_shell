#ifndef _inout_h
#define _inout_h



#include <stdio.h>
#include <stdlib.h>
#include "dataStructures.h"



void sepline();



void cprint(double complex z);



void carr_print(int n, Carray v);
/** Print complex array of size n as column vector on screen **/



void rarr_print(int n, Rarray v);
/** Print real array of size n as column vector on screen **/



void cmat_print(int m, int n, Cmatrix M);
/** Print complex matrix M of m rows and n columns on screen **/



void rmat_print(int m, int n, Rmatrix M);
/** Print real matrix M of m rows and n columns on screen **/



void carr_txt(char fname [], int M, Carray v);
/** Record complex array v of M elements as column vector in file
  * called fname **/



void rarr_txt(char fname [], int M, Rarray v);
/** Record real array v of M elements as column vector in file
  * called fname **/



void carr_inline(FILE * f, int M, Carray v);
/** Record complex array v of M elements on a opened file in the
  * current buffer line. After done print linebreak **/



void rarr_inline(FILE * f, int M, Rarray v);
/** Record real array v of M elements on a opened file in the
  * current buffer line. After done print linebreak **/


#endif
