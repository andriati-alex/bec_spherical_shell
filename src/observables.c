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
        dtht,
        nabla_coef,
        ga,
        gb,
        gab,
        sin_tht,
        cos_tht;
    double complex
        sum_grid_res;
    Rarray
        tht,
        abs_square_a,
        abs_square_b;
    Carray
        grid_res,
        nabla2phi,
        der2tht,
        der_tht;

    // unpack equation data
    tht = EQ->theta;
    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    dtht = tht[1] - tht[0];
    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;

    nabla2phi = carrDef(nphi * ntht);
    der2tht = carrDef(nphi * ntht);
    der_tht = carrDef(nphi * ntht);
    abs_square_a = rarrDef(nphi * ntht);
    abs_square_b = rarrDef(nphi * ntht);
    grid_res = carrDef(nphi * ntht);

    sum_grid_res = 0;

    carrAbs2(nphi * ntht, state_a, abs_square_a);
    carrAbs2(nphi * ntht, state_b, abs_square_b);

    // Compute grid residue for species A

    sph_theta_twice_derivative(nphi, ntht, state_a, dtht, der2tht);
    sph_theta_derivative(nphi, ntht, state_a, dtht, der_tht);
    sph_phi_twice_derivative(nphi, ntht, state_a, dphi, nabla2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            nabla2phi[grid_pt] = nabla2phi[grid_pt] / sin_tht / sin_tht;
        }
    }

    for (j = 0; j < ntht; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res[grid_pt] = (
                    nabla_coef * (
                    sin_tht * der2tht[grid_pt] +
                    cos_tht * der_tht[grid_pt] +
                    sin_tht * nabla2phi[grid_pt]
                    )
                    + sin_tht * state_a[grid_pt] * (
                    ga * abs_square_a[grid_pt] + gab * abs_square_b[grid_pt]
                    )
                    - sin_tht * mu_a * state_a[grid_pt]
            );
        }
    }

    sum_grid_res = Csimps2D(nphi, ntht, grid_res, dphi, dtht);

    // Compute grid residue for species B

    sph_theta_twice_derivative(nphi, ntht, state_b, dtht, der2tht);
    sph_theta_derivative(nphi, ntht, state_b, dtht, der_tht);
    sph_phi_twice_derivative(nphi, ntht, state_b, dphi, nabla2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            nabla2phi[grid_pt] = nabla2phi[grid_pt] / sin_tht / sin_tht;
        }
    }

    for (j = 0; j < ntht; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res[grid_pt] = (
                    nabla_coef * (
                    sin_tht * der2tht[grid_pt] +
                    cos_tht * der_tht[grid_pt] +
                    sin_tht * nabla2phi[grid_pt]
                    )
                    + sin_tht * state_b[grid_pt] * (
                    gb * abs_square_b[grid_pt] + gab * abs_square_a[grid_pt]
                    )
                    - sin_tht * mu_b * state_b[grid_pt]
            );
        }
    }

    sum_grid_res = sum_grid_res + Csimps2D(nphi, ntht, grid_res, dphi, dtht);

    free(nabla2phi);
    free(der2tht);
    free(der_tht);
    free(abs_square_a);
    free(abs_square_b);
    free(grid_res);

    return cabs(sum_grid_res);
}


void abs_grad_sphere_square(
        int nphi, int ntheta, double dphi,
        Rarray theta, Carray state, Rarray abs_grad_square)
{
    int
        j,
        i_phi,
        phi_point,
        pi_rotation;
    double
        dtheta;
    double complex
        z1,
        z2;
    Carray
        ds_dphi,
        grad_phi,
        grad_theta;

    dtheta = theta[1] - theta[0];

    ds_dphi = carrDef(nphi * ntheta);
    grad_phi = carrDef(nphi * ntheta);
    grad_theta = carrDef(nphi * ntheta);

    for (j = 1; j < ntheta - 1; j++)
    {
        // Compute gradient in phi axis
        derivative_1dperiodic(nphi, &state[j*nphi], dphi, &ds_dphi[j*nphi]);
        carrScalarMultiply(
                nphi, &ds_dphi[j*nphi], 1.0 / sin(theta[j]), &grad_phi[j*nphi]
        );
        // Compute gradient in theta axis
        for (i_phi = 0; i_phi < nphi; i_phi++)
        {
            grad_theta[j * nphi + i_phi] = (
                    (state[(j + 1) * nphi + i_phi]
                    - state[(j - 1) * nphi + i_phi]) / (2 * dtheta)
            );
        }
    }

    // compute derivatives at poles theta = 0 and theta = PI following:
    // 1. backward point in theta = 0 is the forward one rotated PI in phi
    // 2. forward point in theta = PI is the backward one rotated PI in phi
    // 3. Avoid exceeding boundary if PI rotation yield phi > 2 * PI
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        phi_point = nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        grad_theta[i_phi] = (
                (state[phi_point] - state[phi_point + pi_rotation])
                / ( 2 * dtheta)
        );
    }
    for (i_phi = 0; i_phi < nphi; i_phi++)
    {
        phi_point = (ntheta - 2) * nphi + i_phi;
        pi_rotation = nphi / 2 - (i_phi / (nphi / 2 + 1)) * (nphi - 1);
        grad_theta[(ntheta - 1) * nphi + i_phi] = (
                (state[phi_point + pi_rotation] - state[phi_point])
                / ( 2 * dtheta)
        );
    }

    // at the poles the function is constant with respect to phi
    carrFill(nphi, 0.0 + 0.0 * I, &grad_phi[0]);
    carrFill(nphi, 0.0 + 0.0 * I, &grad_phi[(ntheta-1)*nphi]);

    for (j = 0; j < ntheta; j++)
    {
        for (i_phi = 0; i_phi < nphi; i_phi++)
        {
            z1 = grad_phi[j * nphi + i_phi];
            z2 = grad_theta[j * nphi + i_phi];
            abs_grad_square[j * nphi + i_phi] = (
                    creal(z1) * creal(z1) + cimag(z1) * cimag(z1)
                    + creal(z2) * creal(z2) + cimag(z2) * cimag(z2)
            );
        }
    }

    free(ds_dphi);
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
            - omega * lza
    );
    (* mu_b) = (
            Rsimps2D_sphere(nphi, ntheta, theta, Integ_mu_b, dphi)
            - omega * lzb
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
