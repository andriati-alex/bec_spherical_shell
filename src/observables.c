#include "observables.h"


double angular_momentum_lz(
        int nphi, int ntheta, double dphi, Rarray theta, Carray state)
{
    int
        j;
    double
        lz;
    Carray
        integ,
        ds_dphi;

    ds_dphi = carrDef(nphi * ntheta);
    integ = carrDef(nphi * ntheta);

    for (j = 1; j < ntheta - 1; j++)
    {
        derivative_1dperiodic(nphi, &state[j*nphi], dphi, &ds_dphi[j*nphi]);
    }


    // at the poles the function is constant with respect to phi
    carrFill(nphi, 0.0 + 0.0 * I, &ds_dphi[0]);
    carrFill(nphi, 0.0 + 0.0 * I, &ds_dphi[(ntheta-1)*nphi]);

    for (j = 0; j < nphi * ntheta; j++)
    {
        integ[j] = conj(state[j]) * ds_dphi[j];
    }

    lz = creal(-I * Csimps2D_sphere(nphi, ntheta, theta, integ, dphi));
    free(ds_dphi);
    free(integ);

    return lz;
}


double avg_residue(
        EqDataPkg EQ, Carray state_a, Carray state_b,
        double mu_a, double mu_b)
{
    int
        j,
        i,
        nphi,
        ntht,
        grid_pt;
    double
        dphi,
        nabla_coef,
        ga,
        gb,
        gab,
        sum_grid_res;
    Rarray
        tht,
        abs_square_a,
        abs_square_b;
    Carray
        grid_res,
        laplace_a,
        laplace_b;

    // unpack equation data
    tht = EQ->theta;
    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;

    abs_square_a = rarrDef(nphi * ntht);
    abs_square_b = rarrDef(nphi * ntht);
    grid_res = carrDef(nphi * ntht);
    laplace_a = carrDef(nphi * ntht);
    laplace_b = carrDef(nphi * ntht);

    sum_grid_res = 0;

    carrAbs2(nphi * ntht, state_a, abs_square_a);
    carrAbs2(nphi * ntht, state_b, abs_square_b);

    // Compute grid residue for species A

    laplace_app(EQ, state_a, state_b, laplace_a, laplace_b);

    for (j = 0; j < ntht; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res[grid_pt] = (
                    nabla_coef * laplace_a[grid_pt]
                    + state_a[grid_pt] * (
                    ga * abs_square_a[grid_pt] + gab * abs_square_b[grid_pt]
                    )
                    - mu_a * state_a[grid_pt]
            );
        }
    }

    sum_grid_res = cabs(Csimps2D_sphere(nphi, ntht, tht, grid_res, dphi));

    for (j = 0; j < ntht; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res[grid_pt] = (
                    nabla_coef * laplace_b[grid_pt]
                    + state_b[grid_pt] * (
                    gb * abs_square_b[grid_pt] + gab * abs_square_a[grid_pt]
                    )
                    - mu_b * state_b[grid_pt]
            );
        }
    }

    sum_grid_res += cabs(Csimps2D_sphere(nphi, ntht, tht, grid_res, dphi));

    free(abs_square_a);
    free(abs_square_b);
    free(grid_res);
    free(laplace_a);
    free(laplace_b);

    return sum_grid_res;
}


void abs_grad_sphere_square(
        int nphi, int ntheta, double dphi,
        Rarray theta, Carray state, Rarray abs_grad_square)
{
    int
        j,
        i,
        grid_pt;
    double
        sin_tht,
        dtheta;
    double complex
        z1,
        z2;
    Carray
        grad_phi,
        grad_theta;

    dtheta = theta[1] - theta[0];

    grad_phi = carrDef(nphi * ntheta);
    grad_theta = carrDef(nphi * ntheta);

    sph_theta_derivative(nphi, ntheta, state, dtheta, grad_theta);
    sph_phi_derivative(nphi, ntheta, state, dphi, grad_phi);

    for (j = 1; j < ntheta - 1; j++)
    {
        sin_tht = sin(theta[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grad_phi[grid_pt] = grad_phi[grid_pt] / sin_tht;
        }
    }

    for (j = 0; j < ntheta; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            z1 = grad_phi[grid_pt];
            z2 = grad_theta[grid_pt];
            abs_grad_square[grid_pt] = (
                    creal(z1) * creal(z1) + cimag(z1) * cimag(z1)
                    + creal(z2) * creal(z2) + cimag(z2) * cimag(z2)
            );
        }
    }

    free(grad_phi);
    free(grad_theta);
}


double functionals_single(
        EqDataPkg EQ, Carray state, double * kin, double * mu)
{

    int
        j,
        i,
        nphi,
        ntheta,
        grid_points;
    double
        total_energy,
        g,
        lz,
        dphi,
        omega,
        nabla_coef;
    Rarray
        theta,
        Integ_energy,
        Integ_kin,
        Integ_mu,
        abs_grad_square,
        abs_square;

    nphi = EQ->nphi;
    dphi = EQ->dphi;
    ntheta = EQ->ntheta;
    theta = EQ->theta;
    g = EQ->ga;
    omega = EQ->omega;
    nabla_coef = EQ->nabla_coef;
    grid_points = nphi * ntheta;

    Integ_energy = rarrDef(grid_points);
    Integ_kin = rarrDef(grid_points);
    Integ_mu = rarrDef(grid_points);
    abs_square = rarrDef(grid_points);
    abs_grad_square = rarrDef(grid_points);

    lz = angular_momentum_lz(nphi, ntheta, dphi, theta, state);
    carrAbs2(grid_points, state, abs_square);
    abs_grad_sphere_square(nphi, ntheta, dphi, theta, state, abs_grad_square);

    // Setup function to integrate
    for (j = 0; j < ntheta; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            Integ_kin[j * nphi + i] = (
                    - nabla_coef * abs_grad_square[j * nphi + i]
            );
            Integ_mu[j * nphi + i] = (
                    - nabla_coef * abs_grad_square[j * nphi + i]
                    + g * (abs_square[j * nphi + i]
                    * abs_square[j * nphi + i])
            );
            Integ_energy[j * nphi + i] = (
                    - nabla_coef * abs_grad_square[j * nphi + i]
                    + 0.5 * g * (abs_square[j * nphi + i]
                    * abs_square[j * nphi + i])
            );
        }
    }

    total_energy = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_energy, dphi)
            - omega * lz
    );
    (* kin) = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_kin, dphi)
            - omega * lz
    );
    (* mu) = Rsimps2D_sphere(nphi, ntheta, theta, Integ_mu, dphi) - omega * lz;

    free(abs_grad_square);
    free(abs_square);
    free(Integ_mu);
    free(Integ_kin);
    free(Integ_energy);

    return total_energy;
}


double functionals(EqDataPkg EQ, Carray state_a, Carray state_b,
        double * kin, double * mu_a, double * mu_b)
{

/** Gross-Pitaesvkii functional **/

    int
        j,
        i,
        nphi,
        ntheta,
        grid_points;
    double
        total_energy,
        ga,
        gb,
        gab,
        lza,
        lzb,
        frac_a,
        frac_b,
        dphi,
        omega,
        nabla_coef;
    Rarray
        theta,
        Integ,
        Integ_kin,
        Integ_mu_a,
        Integ_mu_b,
        abs_grad_square_a,
        abs_grad_square_b,
        abs_square_a,
        abs_square_b;

    nphi = EQ->nphi;
    dphi = EQ->dphi;
    ntheta = EQ->ntheta;
    theta = EQ->theta;
    frac_a = EQ->frac_a;
    frac_b = EQ->frac_b;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;
    omega = EQ->omega;
    nabla_coef = EQ->nabla_coef;
    grid_points = nphi * ntheta;

    Integ = rarrDef(grid_points);
    Integ_kin = rarrDef(grid_points);
    Integ_mu_a = rarrDef(grid_points);
    Integ_mu_b = rarrDef(grid_points);
    abs_square_a = rarrDef(grid_points);
    abs_square_b = rarrDef(grid_points);
    abs_grad_square_a = rarrDef(grid_points);
    abs_grad_square_b = rarrDef(grid_points);

    lza = angular_momentum_lz(nphi, ntheta, dphi, theta, state_a);
    lzb = angular_momentum_lz(nphi, ntheta, dphi, theta, state_b);
    carrAbs2(grid_points, state_a, abs_square_a);
    carrAbs2(grid_points, state_b, abs_square_b);
    abs_grad_sphere_square(
            nphi, ntheta, dphi, theta, state_a, abs_grad_square_a
    );
    abs_grad_sphere_square(
            nphi, ntheta, dphi, theta, state_b, abs_grad_square_b
    );

    // Setup function to integrate
    for (j = 0; j < ntheta; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            Integ_kin[j * nphi + i] = (
                    - nabla_coef * frac_a * abs_grad_square_a[j * nphi + i]
                    - nabla_coef * frac_b * abs_grad_square_b[j * nphi + i]
            );
            Integ_mu_a[j * nphi + i] = (
                    - nabla_coef * abs_grad_square_a[j * nphi + i]
                    + ga * (abs_square_a[j * nphi + i]
                    * abs_square_a[j * nphi + i])
                    + sqrt(frac_b / frac_a) * gab * (abs_square_a[j * nphi + i]
                    * abs_square_b[j * nphi + i])
            );
            Integ_mu_b[j * nphi + i] = (
                    - nabla_coef * abs_grad_square_b[j * nphi + i]
                    + gb * (abs_square_b[j * nphi + i]
                    * abs_square_b[j * nphi + i])
                    + sqrt(frac_a / frac_b) * gab * (abs_square_a[j * nphi + i]
                    * abs_square_b[j * nphi + i])
            );
            Integ[j * nphi + i] = (
                    - nabla_coef * frac_a * abs_grad_square_a[j * nphi + i]
                    - nabla_coef * frac_b * abs_grad_square_b[j * nphi + i]
                    + 0.5 * frac_a * ga * (abs_square_a[j * nphi + i]
                    * abs_square_a[j * nphi + i])
                    + 0.5 * frac_b * gb * (abs_square_b[j * nphi + i]
                    * abs_square_b[j * nphi + i])
                    + sqrt(frac_a * frac_b) * gab * (abs_square_a[j * nphi + i]
                    * abs_square_b[j * nphi + i])
            );
        }
    }

    total_energy = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ, dphi)
            - omega * (lza + lzb)
    );
    (* kin) = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_kin, dphi)
            - omega * (lza + lzb)
    );
    (* mu_a) = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_mu_a, dphi)
    );
    (* mu_b) = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_mu_b, dphi)
    );

    free(abs_grad_square_a);
    free(abs_grad_square_b);
    free(abs_square_a);
    free(abs_square_b);
    free(Integ_mu_a);
    free(Integ_mu_b);
    free(Integ_kin);
    free(Integ);

    return total_energy;
}


double functionals_theta(EqDataPkg EQ, Carray state_a, Carray state_b,
        double * kin, double * mu_a, double * mu_b, int azi_a, int azi_b)
{

    int
        j,
        ntheta;
    double
        sin_th,
        dtheta,
        nabla_part_a,
        nabla_part_b,
        total_energy,
        ga,
        gb,
        gab,
        frac_a,
        frac_b,
        nabla_coef;
    Carray
        der_a,
        der_b;
    Rarray
        theta,
        Integ,
        Integ_kin,
        Integ_mu_a,
        Integ_mu_b,
        abs_square_der_a,
        abs_square_der_b,
        abs_square_a,
        abs_square_b;

    ntheta = EQ->ntheta;
    dtheta = EQ->dtheta;
    theta = EQ->theta;
    frac_a = EQ->frac_a;
    frac_b = EQ->frac_b;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;
    nabla_coef = EQ->nabla_coef;

    der_a = carrDef(ntheta);
    der_b = carrDef(ntheta);

    Integ = rarrDef(ntheta);
    Integ_kin = rarrDef(ntheta);
    Integ_mu_a = rarrDef(ntheta);
    Integ_mu_b = rarrDef(ntheta);
    abs_square_a = rarrDef(ntheta);
    abs_square_b = rarrDef(ntheta);
    abs_square_der_a = rarrDef(ntheta);
    abs_square_der_b = rarrDef(ntheta);

    derivative_1dreflection(ntheta, state_a, dtheta, der_a, azi_a);
    derivative_1dreflection(ntheta, state_b, dtheta, der_b, azi_b);
    carrAbs2(ntheta, der_a, abs_square_der_a);
    carrAbs2(ntheta, der_b, abs_square_der_b);

    carrAbs2(ntheta, state_a, abs_square_a);
    carrAbs2(ntheta, state_b, abs_square_b);

    // Setup function to integrate
    for (j = 1; j < ntheta - 1; j++)
    {
        sin_th = sin(theta[j]);
        nabla_part_a = -nabla_coef * (
                abs_square_der_a[j] +
                azi_a * azi_a * abs_square_a[j] / sin_th / sin_th
        );
        nabla_part_b = -nabla_coef * (
                abs_square_der_b[j] +
                azi_b * azi_b * abs_square_b[j] / sin_th / sin_th
        );
        Integ_kin[j] = sin_th * (frac_a * nabla_part_a + frac_b * nabla_part_b);
        Integ_mu_a[j] = sin_th * (
                nabla_part_a
                + ga * abs_square_a[j] * abs_square_a[j]
                + sqrt(frac_b / frac_a) * gab * (abs_square_a[j]
                * abs_square_b[j])
        );
        Integ_mu_b[j] = sin_th * (
                nabla_part_b
                + gb * abs_square_b[j] * abs_square_b[j]
                + sqrt(frac_a / frac_b) * gab * (abs_square_a[j]
                * abs_square_b[j])
        );
        Integ[j] = sin_th * (
                frac_a * nabla_part_a
                + frac_b * nabla_part_b
                + 0.5 * frac_a * ga * abs_square_a[j] * abs_square_a[j]
                + 0.5 * frac_b * gb * abs_square_b[j] * abs_square_b[j]
                + sqrt(frac_a * frac_b) * gab * (abs_square_a[j]
                * abs_square_b[j])
        );
    }
    Integ[0] = 0.0;
    Integ[ntheta - 1] = 0.0;
    Integ_kin[0] = 0.0;
    Integ_kin[ntheta - 1] = 0.0;
    Integ_mu_a[0] = 0.0;
    Integ_mu_a[ntheta - 1] = 0.0;
    Integ_mu_b[0] = 0.0;
    Integ_mu_b[ntheta - 1] = 0.0;

    total_energy = Rsimps1D(ntheta, Integ, dtheta);
    (* kin) = Rsimps1D(ntheta, Integ_kin, dtheta);
    (* mu_a) = Rsimps1D(ntheta, Integ_mu_a, dtheta);
    (* mu_b) = Rsimps1D(ntheta, Integ_mu_b, dtheta);

    free(abs_square_der_a);
    free(abs_square_der_b);
    free(abs_square_a);
    free(abs_square_b);
    free(Integ_mu_a);
    free(Integ_mu_b);
    free(Integ_kin);
    free(Integ);
    free(der_a);
    free(der_b);

    return total_energy;
}


double theta_density_overlap(EqDataPkg EQ, Rarray density_a, Rarray density_b)
{
    int
        i;
    double
        norm_a,
        norm_b,
        raw_overlap;
    Rarray
        integ_a,
        integ_b,
        integ_ab;

    integ_ab = rarrDef(EQ->ntheta);
    integ_a = rarrDef(EQ->ntheta);
    integ_b = rarrDef(EQ->ntheta);

    for (i = 0; i < EQ->ntheta; i++)
    {
        integ_ab[i] = sin(EQ->theta[i]) * density_a[i] * density_b[i];
        integ_a[i] = sin(EQ->theta[i]) * density_a[i] * density_a[i];
        integ_b[i] = sin(EQ->theta[i]) * density_b[i] * density_b[i];
    }

    raw_overlap = Rsimps1D(EQ->ntheta, integ_ab, EQ->dtheta);
    norm_a = Rsimps1D(EQ->ntheta, integ_a, EQ->dtheta);
    norm_b = Rsimps1D(EQ->ntheta, integ_b, EQ->dtheta);

    free(integ_ab);
    free(integ_a);
    free(integ_b);
    return raw_overlap * raw_overlap / norm_a / norm_b;
}


double density_overlap(EqDataPkg EQ, Rarray density_a, Rarray density_b)
{
    int
        i;
    double
        norm_a,
        norm_b,
        raw_overlap;
    Rarray
        integ_a,
        integ_b,
        integ_ab;

    integ_ab = rarrDef(EQ->ntheta * EQ->nphi);
    integ_a = rarrDef(EQ->ntheta * EQ->nphi);
    integ_b = rarrDef(EQ->ntheta * EQ->nphi);

    for (i = 0; i < EQ->ntheta * EQ->nphi; i++)
    {
        integ_ab[i] = density_a[i] * density_b[i];
        integ_a[i] = density_a[i] * density_a[i];
        integ_b[i] = density_b[i] * density_b[i];
    }

    raw_overlap = Rsimps2D_sphere(
            EQ->nphi, EQ->ntheta, EQ->theta, integ_ab, EQ->dphi
    );
    norm_a = Rsimps2D_sphere(
            EQ->nphi, EQ->ntheta, EQ->theta, integ_a, EQ->dphi
    );
    norm_b = Rsimps2D_sphere(
            EQ->nphi, EQ->ntheta, EQ->theta, integ_b, EQ->dphi
    );

    free(integ_ab);
    free(integ_a);
    free(integ_b);
    return raw_overlap * raw_overlap / norm_a / norm_b;
}
