#include "tridiagonalSystems.h"



void tridiag(int n, Carray upper, Carray lower, Carray mid, Carray RHS,
     Carray ans)
{

/** Solve a tridiagonal system of equations. Output parameter : ans **/

    unsigned int
        i,
        k;

    double complex
        RHS1,
        RHS2;

    // Intermediate steps - L . U decomposition

    Carray
        u,
        l,
        z;



    u = carrDef(n);
    l = carrDef(n - 1);
    z = carrDef(n);

    if (cabs(mid[0]) == 0 )
    {
        // In this case there is a system reduction
        // where we solve for [x1  x3  x4 ... xn]
        // what is equivalent to adjust the two first 
        // equations and starts counters from 1

        RHS1 = RHS[1] - mid[1]   * RHS[0] / upper[0];
        RHS2 = RHS[2] - lower[1] * RHS[0] / upper[0];

        // u and z factor initizlization with the changed system
        u[1] = lower[0];
        z[1] = RHS1;

        // One iteration need to be performed out because
        // the change of values in lower and RHS

        l[1] = 0;
        u[2] = mid[2] - l[1] * upper[1];
        z[2] = RHS2 - l[1] * z[1];

        for (i = 2; i < n - 1; i++)
        {
            k = i + 1;
            l[i] = lower[i] / u[i];
            u[k] = mid[k] - l[i] * upper[i];
            z[k] = RHS[k] - l[i] * z[i];
        }

        ans[n-1] = z[n-1] / u[n-1];

        for (i = 2; i <= n - 1; i++)
        {
            k = n - i;
            ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
        }

        // Obtained order ans[0..n] = [nan  x1  x3  x4 .. xn]
        // Organize ans[0..n] = [x1  x2  x3  .. xn]
        ans[0] = ans[1];
        ans[1] = RHS[0] / upper[0];

        // Free local alilocated memory
        free(u);
        free(l);
        free(z);

        return;
    }
    
    u[0] = mid[0];
    z[0] = RHS[0];

    for (i = 0;  i < n - 1; i++)
    {
        k = i + 1;
        l[i] = lower[i] / u[i];
        u[k] = mid[k] - l[i] * upper[i];
        z[k] = RHS[k] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];

    for (i = 2; i <= n; i++)
    {
        k = n - i;
        ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
    }

    // Free local allocated memory
    free(u);
    free(l);
    free(z);
}


void lu_decomposition(
        int n, Carray upper, Carray lower, Carray mid, Carray l, Carray u)
{

/** Compute vectors of LU decomposition for tridigonal system **/

    unsigned int
        i;

    u[0] = mid[0];
    for (i = 0;  i < n - 1; i++)
    {
        l[i] = lower[i] / u[i];
        u[i + 1] = mid[i + 1] - l[i] * upper[i];
    }

}


void tridiag_lu(
        int n, Carray upper, Carray l, Carray u, Carray z,
        Carray RHS, Carray ans)
{

/** Improved routine to solve tridiagonal system given lu-decomposition **/

    int
        i;

    z[0] = RHS[0];

    for (i = 0;  i < n - 1; i++)
    {
        z[i + 1] = RHS[i + 1] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];
    for (i = n - 2; i >= 0; i--)
    {
        ans[i] = (z[i] - upper[i] * ans[i + 1]) / u[i];
    }
}
