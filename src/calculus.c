#include "calculus.h"


double complex Csimps1D(int n, Carray f, double h)
{

    int
        i;

    double complex
        sum;

    sum = 0;

    if (n < 3)
    {
        printf("\n\n\tERROR : less than 3 point to integrate by simps !\n\n");
        exit(EXIT_FAILURE);
    }

    if (n % 2 == 0)
    {

    //  Case the number of points is even then must integrate the last
    //  chunk using simpson's 3/8 rule to maintain accuracy

        for (i = 0; i < (n - 4); i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals
        sum = sum + (f[n-4] + 3 * (f[n-3] + f[n-2]) + f[n-1]) * 3 * h / 8;

    }

    else
    {

        for (i = 0; i < n - 2; i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals

    }

    return sum;

}





double Rsimps1D(int n, Rarray f, double h)
{

    int
        i;

    double
        sum;

    sum = 0;

    if (n % 2 == 0)
    {

    //  Case the number of points is even then must integrate the last
    //  chunk using simpson's 3/8 rule to maintain accuracy

        for (i = 0; i < (n - 4); i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals
        sum = sum + (f[n-4] + 3 * (f[n-3] + f[n-2]) + f[n-1]) * 3 * h / 8;

    }

    else
    {

        for (i = 0; i < n - 2; i = i + 2)
        {
            sum = sum + f[i] + 4 * f[i + 1] + f[i + 2];
        }
        sum = sum * h / 3; // End 3-point simpsons intervals

    }

    return sum;

}





double complex Csimps2D(int nx, int ny, Carray f, double hx, double hy)
{

/** INTEGRATION OF A FUNCTION OF 2 VARIABLES
  * ----------------------------------------
  *
  * Here the function must be stored in a vector, where given a
  * point in the grid (xi,yj) then
  *
  *     funtion(xi,yj) = f[i + nx * j]
  *
  * where i = 0, 1, ..., nx-1 and j = 0, 1, ..., ny-1 and the grid points
  * are linearly spaced xi = x0 + i*hx as well as yj = y0 + j*hy      **/

    unsigned int
        j;

    double complex
        result;

    Carray
        fy;

    fy = carrDef(ny);

#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < ny; j++)
    {
        // Integrate in x-direction and end up with a function of y
        fy[j] = Csimps1D(nx,&f[j*nx],hx);
    }

    result = Csimps1D(ny,fy,hy);

    free(fy);

    return result;

}





double Rsimps2D(int nx, int ny, Rarray f, double hx, double hy)
{

/** INTEGRATION OF A FUNCTION OF 2 VARIABLES
  * ----------------------------------------
  *
  * Here the function must be stored in a vector, where given point a
  * point in the grid (xi,yj) then
  *
  *     funtion(xi,yj) = f[i + nx * j]
  *
  * where i = 0, 1, ..., nx-1 and j = 0, 1, ..., ny-1 and the grid points
  * are linearly spaced xi = x0 + i*hx as well as yj = y0 + j*hy      **/

    unsigned int
        j;

    double
        result;

    Rarray
        fy;

    fy = rarrDef(ny);

#pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < ny; j++)
    {
        // Integrate in x-direction and end up with a function of y
        fy[j] = Rsimps1D(nx,&f[j*nx],hx);
    }

    result = Rsimps1D(ny,fy,hy);

    free(fy);

    return result;

}


double Rsimps2D_sphere(
        int nphi, int ntheta, Rarray theta, Rarray f, double dphi)
{

/** INTEGRATE FUNCTION OF 2 VARIABLES IN SPHERICAL COORDINATES **/

    unsigned int
        j;
    double
        dtheta,
        result;
    Rarray
        f_theta;

    dtheta = theta[1] - theta[0];
    f_theta = rarrDef(ntheta);

//  #pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < ntheta; j++)
    {
        // Integrate in phi and multiply by the jacobian term sin(theta)
        f_theta[j] = sin(theta[j]) * Rsimps1D(nphi,&f[j*nphi],dphi);
    }

    result = Rsimps1D(ntheta,f_theta,dtheta);

    free(f_theta);
    return result;
}


double complex Csimps2D_sphere(
        int nphi, int ntheta, Rarray theta, Carray f, double dphi)
{

/** INTEGRATE FUNCTION OF 2 VARIABLES IN SPHERICAL COORDINATES **/

    unsigned int
        j;
    double
        dtheta,
        result;
    Carray
        f_theta;

    dtheta = theta[1] - theta[0];
    f_theta = carrDef(ntheta);

//  #pragma omp parallel for private(j) schedule(static)
    for (j = 0; j < ntheta; j++)
    {
        // Integrate in phi and multiply by the jacobian term sin(theta)
        f_theta[j] = sin(theta[j]) * Csimps1D(nphi,&f[j*nphi],dphi);
    }

    result = Csimps1D(ntheta,f_theta,dtheta);

    free(f_theta);
    return result;
}


void renormalize_spheric(EqDataPkg EQ, Carray S)
{
    int
        i,
        grid_points;
    double
        old_norm;
    Rarray
        abs_square;

    grid_points = EQ->nphi * EQ->ntheta;
    abs_square = rarrDef(grid_points);
    carrAbs2(grid_points, S, abs_square);

    old_norm = Rsimps2D_sphere(EQ->nphi, EQ->ntheta, EQ->theta,
            abs_square, EQ->dphi
    );

    for (i = 0; i < grid_points; i++) S[i] = S[i] / sqrt(old_norm);
    free(abs_square);
}



void renormalize(int nx, int ny, Carray f, double hx, double hy, double norm)
{

    int
        i;

    double
        renorm;

    Rarray
        ModSquared;

    ModSquared = rarrDef(nx*ny);

    carrAbs2(nx*ny,f,ModSquared);

    renorm = norm * sqrt(1.0 / Rsimps2D(nx,ny,ModSquared,hx,hy));

    for (i = 0; i < nx*ny; i++) f[i] = f[i] * renorm;

    free(ModSquared);
}





void renormalizeReal(int nx, int ny, Rarray f, double hx, double hy,
                     double norm)
{

    int
        i;

    double
        renorm;

    Rarray
        f2;

    f2 = rarrDef(nx*ny);

    for (i = 0; i < nx*ny; i++) f2[i] = f[i] * f[i];

    renorm = norm * sqrt(1.0 / Rsimps2D(nx,ny,f2,hx,hy));

    for (i = 0; i < nx*ny; i++) f[i] = f[i] * renorm;

    free(f2);
}


void derivative_periodic(int Npts, Carray f, double dx, Carray dfdx)
{

/** Compute derivative of function `f` in `Npts` grid points with
    spacing `dx` considering  periodic boundary conditions,  that
    is f[n-1] = f[0],  with 4th-order Finite-Differences accuracy
    Output parameter : `dfdx`                                 **/

    int
        n,
        i;
    double
        r;

    n = Npts;            // make the life easier
    r = 1.0 / (12 * dx); // ratio for a fourth-order scheme

    // COMPUTE USING PERIODIC BOUNDARY CONDITIONS

    dfdx[0]   = ( f[n-3] - f[2] + 8 * (f[1] - f[n-2]) ) * r;
    dfdx[1]   = ( f[n-2] - f[3] + 8 * (f[2] - f[0]) )   * r;
    dfdx[n-2] = ( f[n-4] - f[1] + 8 * (f[0] - f[n-3]) ) * r;
    dfdx[n-1] = dfdx[0]; // assume last point as the boundary

    for (i = 2; i < n - 2; i++)
    {
        dfdx[i] = ( f[i-2] - f[i+2] + 8 * (f[i+1] - f[i-1]) ) * r;
    }

}
