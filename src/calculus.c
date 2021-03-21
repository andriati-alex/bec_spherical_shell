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

    #pragma omp parallel for private(j) schedule(static)
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
        dtheta;
    double complex
        result;
    Carray
        f_theta;

    dtheta = theta[1] - theta[0];
    f_theta = carrDef(ntheta);

    #pragma omp parallel for private(j) schedule(static)
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


void renormalize_theta_sphere(EqDataPkg EQ, Carray S)
{
    int
        i;
    double
        old_norm;
    Rarray
        abs_square_sin;

    abs_square_sin = rarrDef(EQ->ntheta);

    for (i = 0; i < EQ->ntheta; i++)
    {
        abs_square_sin[i] = sin(EQ->theta[i]) * (
                creal(S[i]) * creal(S[i]) + cimag(S[i]) * cimag(S[i])
        );
    }

    old_norm = Rsimps1D(EQ->ntheta, abs_square_sin, EQ->dtheta);

    for (i = 0; i < EQ->ntheta; i++) S[i] = S[i] / sqrt(old_norm);
    free(abs_square_sin);
}


void derivative_1dperiodic(int n, Carray f, double dx, Carray dfdx)
{

/** Compute derivative of function `f` in `n` grid  points  with
    spacing `dx` considering  periodic boundary conditions, that
    is f[n-1] = f[0],  with 4th-order Finite-Differences  method
    Output parameter : `dfdx`                                **/

    int
        i;
    double
        r;

    r = 1.0 / (12 * dx); // ratio for a fourth-order scheme

    dfdx[0]   = ( f[n-3] - f[2] + 8 * (f[1] - f[n-2]) ) * r;
    dfdx[1]   = ( f[n-2] - f[3] + 8 * (f[2] - f[0]) )   * r;
    dfdx[n-2] = ( f[n-4] - f[1] + 8 * (f[0] - f[n-3]) ) * r;
    dfdx[n-1] = dfdx[0]; // assume last point as the boundary

    for (i = 2; i < n - 2; i++)
    {
        dfdx[i] = ( f[i-2] - f[i+2] + 8 * (f[i+1] - f[i-1]) ) * r;
    }

}


void derivative_1dreflection(int n, Carray f, double dx, Carray dfdx, int azi)
{

    int
        i,
        sign;
    double
        r;

    sign = 1 - 2 * (abs(azi) % 2);

    r = 1.0 / (12 * dx); // ratio for a fourth-order scheme

    dfdx[0]   = ( sign * f[2] - f[2] + 8 * (f[1] - sign * f[1]) ) * r;
    dfdx[1]   = ( sign * f[1] - f[3] + 8 * (f[2] - f[0]) ) * r;
    dfdx[n-2] = ( f[n-4] - sign * f[n-2] + 8 * (f[n-1] - f[n-3]) ) * r;
    dfdx[n-1] = ( f[n-3] - sign * f[n-3] + 8 * (sign * f[n-2] - f[n-2]) ) * r;

    for (i = 2; i < n - 2; i++)
    {
        dfdx[i] = ( f[i-2] - f[i+2] + 8 * (f[i+1] - f[i-1]) ) * r;
    }

}


void twice_derivative_1dperiodic(int n, Carray f, double dx, Carray d2fdx)
{

/** Compute second derivative of function `f` in `n` grid points
    with spacing `dx`. Consider f[0] = f[n - 1] the boundary **/

    int
        i,
        l;

    double
        r;

    r = 1.0 / (12 * dx * dx);

    l = n - 2; // last index before boundary (at n - 1 => 0)
    d2fdx[0] = r * (- f[2] + 16 * f[1] - 30 * f[0] + 16 * f[l] - f[l - 1]);
    d2fdx[1] = r * (- f[3] + 16 * f[2] - 30 * f[1] + 16 * f[0] - f[l]);
    d2fdx[l] = r * (- f[1] + 16 * f[0] - 30 * f[l] + 16 * f[l - 1] - f[l - 2]);
    d2fdx[n - 1] = d2fdx[0];

    for (i = 2; i < n - 2; i++)
    {
        d2fdx[i] = r * (
                - f[i + 2] + 16 * f[i + 1] - 30 * f[i]
                + 16 * f[i - 1] - f[i - 2]
        );
    }

}


void twice_derivative_1dreflection(
        int n, Carray f, double dx, Carray d2fdx, int azi)
{

/** Compute second derivative of function `f` in `n` grid points
    with spacing `dx`. Consider f[0] = f[n - 1] the boundary **/

    int
        i,
        l,
        sign;
    double
        r;

    sign = 1 - 2 * (abs(azi) % 2);
    r = 1.0 / (12 * dx * dx);

    l = n - 1; // last index
    d2fdx[0] = r * (
            -f[2] + 16 * f[1] - 30 * f[0] + 16 * sign * f[1] - sign * f[2]
    );
    d2fdx[1] = r * (
            -f[3] + 16 * f[2] - 30 * f[1] + 16 * f[0] - sign * f[1]
    );
    d2fdx[l] = r * (
            -sign * f[l-2] + sign * 16 * f[l-1] -
            30 * f[l] + 16 * f[l-1] - f[l-2]);
    d2fdx[l - 1] = r * (
            -sign * f[l-1] + 16 * f[l] - 30 * f[l-1]
            + 16 * f[l-2] - f[l-3]
    );

    for (i = 2; i < n - 2; i++)
    {
        d2fdx[i] = r * (
                - f[i + 2] + 16 * f[i + 1] - 30 * f[i]
                + 16 * f[i - 1] - f[i - 2]
        );
    }

}


void sph_phi_derivative(
        int nphi, int ntht, Carray f, double dphi, Carray der)
{

/** Compute derivative of function with respect to phi in spherical
    coordinates. At the poles the phi-derivative must vanish    **/

    for (int i = 1; i < ntht - 1; i++)
    {
        derivative_1dperiodic(
                nphi, &f[i * nphi], dphi, &der[i * nphi]
        );
    }

    // at the poles the function must be constant with respect to phi
    carrFill(nphi, 0.0 + 0.0 * I, &der[0]);
    carrFill(nphi, 0.0 + 0.0 * I, &der[(ntht - 1) * nphi]);

}


void sph_phi_twice_derivative(
        int nphi, int ntht, Carray f, double dphi, Carray der)
{

/** Compute derivative of function with respect to phi in spherical
    coordinates. At the poles the phi-derivative must vanish    **/

    for (int i = 1; i < ntht - 1; i++)
    {
        twice_derivative_1dperiodic(
                nphi, &f[i * nphi], dphi, &der[i * nphi]
        );
    }

    // at the poles the function must be constant with respect to phi
    carrFill(nphi, 0.0 + 0.0 * I, &der[0]);
    carrFill(nphi, 0.0 + 0.0 * I, &der[(ntht - 1) * nphi]);

}


void sph_theta_derivative(
        int nphi, int ntht, Carray f, double dtht, Carray der)
{
    int
        j,
        i_phi,
        adv,
        bck,
        adv_twice,
        bck_twice,
        pi_rotation;
    double
        r;

    r = 1.0 / (12 * dtht); // ratio for a fourth-order scheme

    for (j = 2; j < ntht - 2; j++)
    {
        for (i_phi = 0; i_phi < nphi; i_phi++)
        {
            adv_twice = (j + 2) * nphi + i_phi;
            adv = (j + 1) * nphi + i_phi;
            bck = (j - 1) * nphi + i_phi;
            bck_twice = (j - 2) * nphi + i_phi;
            der[j * nphi + i_phi] = r * (
                    f[bck_twice] - f[adv_twice] + 8 * (f[adv] - f[bck])
            );
        }
    }

    j = 0;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        adv_twice = (j + 2) * nphi + i_phi;
        adv = (j + 1) * nphi + i_phi;
        bck = adv + pi_rotation;
        bck_twice = adv_twice + pi_rotation;
        der[j * nphi + i_phi] = r * (
                f[bck_twice] - f[adv_twice] + 8 * (f[adv] - f[bck])
        );
    }

    j = 1;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        adv_twice = (j + 2) * nphi + i_phi;
        adv = (j + 1) * nphi + i_phi;
        bck = (j - 1) * nphi + i_phi;
        bck_twice = j * nphi + i_phi + pi_rotation;
        der[j * nphi + i_phi] = r * (
                f[bck_twice] - f[adv_twice] + 8 * (f[adv] - f[bck])
        );
    }

    j = ntht - 1;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        bck = (j - 1) * nphi + i_phi;
        bck_twice = (j - 2) * nphi + i_phi;
        adv_twice = bck_twice + pi_rotation;
        adv = bck + pi_rotation;
        der[j * nphi + i_phi] = r * (
                f[bck_twice] - f[adv_twice] + 8 * (f[adv] - f[bck])
        );
    }

    j = ntht - 2;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        bck = (j - 1) * nphi + i_phi;
        bck_twice = (j - 2) * nphi + i_phi;
        adv_twice = j * nphi + i_phi + pi_rotation;
        adv = (j + 1) * nphi + i_phi;
        der[j * nphi + i_phi] = r * (
                f[bck_twice] - f[adv_twice] + 8 * (f[adv] - f[bck])
        );
    }

    for (j = 0; j < ntht; j++)
    {
        der[j * nphi + nphi - 1] = der[j * nphi];
    }

}


void sph_theta_twice_derivative(
        int nphi, int ntht, Carray f, double dtht, Carray der)
{
    int
        j,
        i_phi,
        pnt,
        adv,
        bck,
        adv_twice,
        bck_twice,
        pi_rotation;
    double
        r;

    r = 1.0 / (12 * dtht * dtht); // ratio for a fourth-order scheme

    for (j = 2; j < ntht - 2; j++)
    {

        for (i_phi = 0; i_phi < nphi; i_phi++)
        {
            pnt = j * nphi + i_phi;
            adv_twice = (j + 2) * nphi + i_phi;
            adv = (j + 1) * nphi + i_phi;
            bck = (j - 1) * nphi + i_phi;
            bck_twice = (j - 2) * nphi + i_phi;
            der[pnt] = r * (
                    - f[bck_twice] - f[adv_twice] + 16 * (f[adv] + f[bck])
                    - 30 * f[pnt]
            );
        }
    }

    j = 0;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pnt = j * nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        adv_twice = (j + 2) * nphi + i_phi;
        adv = (j + 1) * nphi + i_phi;
        bck = adv + pi_rotation;
        bck_twice = adv_twice + pi_rotation;
        der[pnt] = r * (
                - f[bck_twice] - f[adv_twice] + 16 * (f[adv] + f[bck])
                - 30 * f[pnt]
        );
    }

    j = 1;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pnt = j * nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        adv_twice = (j + 2) * nphi + i_phi;
        adv = (j + 1) * nphi + i_phi;
        bck = (j - 1) * nphi + i_phi;
        bck_twice = pnt + pi_rotation;
        der[pnt] = r * (
                - f[bck_twice] - f[adv_twice] + 16 * (f[adv] + f[bck])
                - 30 * f[pnt]
        );
    }

    j = ntht - 1;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pnt = j * nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        bck = (j - 1) * nphi + i_phi;
        bck_twice = (j - 2) * nphi + i_phi;
        adv_twice = bck_twice + pi_rotation;
        adv = bck + pi_rotation;
        der[pnt] = r * (
                - f[bck_twice] - f[adv_twice] + 16 * (f[adv] + f[bck])
                - 30 * f[pnt]
        );
    }

    j = ntht - 2;
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        pnt = j * nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        bck = (j - 1) * nphi + i_phi;
        bck_twice = (j - 2) * nphi + i_phi;
        adv_twice = pnt + pi_rotation;
        adv = (j + 1) * nphi + i_phi;
        der[pnt] = r * (
                - f[bck_twice] - f[adv_twice] + 16 * (f[adv] + f[bck])
                - 30 * f[pnt]
        );
    }

    for (j = 0; j < ntht; j++)
    {
        der[j * nphi + nphi - 1] = der[j * nphi];
    }

}



void laplace_app(
        EqDataPkg EQ,
        Carray state_a,
        Carray state_b,
        Carray laplace_a,
        Carray laplace_b)
{
    int
        j,
        i,
        nphi,
        ntht,
        grid_pt;
    double
        dphi,
        dtht,
        sin_tht,
        cos_tht;
    double complex
        first_order,
        second_order,
        avg_pole_a,
        avg_pole_b;
    Rarray
        tht;
    Carray
        der2phi,
        der2tht,
        der_tht;

    // unpack equation data
    tht = EQ->theta;
    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    dtht = tht[1] - tht[0];

    der2phi = carrDef(nphi * ntht);
    der2tht = carrDef(nphi * ntht);
    der_tht = carrDef(nphi * ntht);

    // Compute grid residue for species A

    sph_theta_twice_derivative(nphi, ntht, state_a, dtht, der2tht);
    sph_theta_derivative(nphi, ntht, state_a, dtht, der_tht);
    sph_phi_twice_derivative(nphi, ntht, state_a, dphi, der2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            laplace_a[grid_pt] = (
                    der2tht[grid_pt] +
                    der_tht[grid_pt] * cos_tht / sin_tht +
                    der2phi[grid_pt] / sin_tht / sin_tht
            );
        }
    }

    sph_theta_twice_derivative(nphi, ntht, state_b, dtht, der2tht);
    sph_theta_derivative(nphi, ntht, state_b, dtht, der_tht);
    sph_phi_twice_derivative(nphi, ntht, state_b, dphi, der2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            laplace_b[grid_pt] = (
                    der2tht[grid_pt] +
                    der_tht[grid_pt] * cos_tht / sin_tht +
                    der2phi[grid_pt] / sin_tht / sin_tht
            );
        }
    }

    // approximate residue at the poles using the average
    // over the small ring at the nearest \theta grid point

    avg_pole_a = 0.0;
    avg_pole_b = 0.0;
    for (i = 0; i < nphi; i++)
    {
        grid_pt = nphi + i;
        first_order = - 0.5 * (
                - 3 * laplace_a[grid_pt] +
                4 * laplace_a[grid_pt + nphi] -
                laplace_a[grid_pt + 2 * nphi]
        );
        second_order = 0.5 * (
                laplace_a[grid_pt + 2 * nphi] -
                2 * laplace_a[grid_pt + nphi] +
                laplace_a[grid_pt]
        );
        avg_pole_a += (laplace_a[grid_pt] + first_order + second_order) / nphi;
        // compute for b species
        first_order = - 0.5 * (
                - 3 * laplace_b[grid_pt] +
                4 * laplace_b[grid_pt + nphi] -
                laplace_b[grid_pt + 2 * nphi]
        );
        second_order = 0.5 * (
                laplace_b[grid_pt + 2 * nphi] -
                2 * laplace_b[grid_pt + nphi] +
                laplace_b[grid_pt]
        );
        avg_pole_b += (laplace_b[grid_pt] + first_order + second_order) / nphi;
    }
    carrFill(nphi, avg_pole_a, &laplace_a[0]);
    carrFill(nphi, avg_pole_b, &laplace_b[0]);

    avg_pole_a = 0.0;
    avg_pole_b = 0.0;
    for (i = 0; i < nphi; i++)
    {
        grid_pt = (ntht - 2) * nphi + i;
        first_order = 0.5 * (
                3 * laplace_a[grid_pt] -
                4 * laplace_a[grid_pt - nphi] +
                laplace_a[grid_pt - 2 * nphi]
        );
        second_order = 0.5 * (
                laplace_a[grid_pt - 2 * nphi] -
                2 * laplace_a[grid_pt - nphi] +
                laplace_a[grid_pt]
        );
        avg_pole_a += (laplace_a[grid_pt] + first_order + second_order) / nphi;
        // species B
        first_order = 0.5 * (
                3 * laplace_b[grid_pt] -
                4 * laplace_b[grid_pt - nphi] +
                laplace_b[grid_pt - 2 * nphi]
        );
        second_order = 0.5 * (
                laplace_b[grid_pt - 2 * nphi] -
                2 * laplace_b[grid_pt - nphi] +
                laplace_b[grid_pt]
        );
        avg_pole_b += (laplace_b[grid_pt] + first_order + second_order) / nphi;
    }
    carrFill(nphi, avg_pole_a, &laplace_a[(ntht - 1) * nphi]);
    carrFill(nphi, avg_pole_b, &laplace_b[(ntht - 1) * nphi]);

    free(der2phi);
    free(der2tht);
    free(der_tht);
}
