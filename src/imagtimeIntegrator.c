#include "imagtimeIntegrator.h"


void _explicit_theta(EqDataPkg EQ, double azi_num, Carray psi, Carray cn_psi)
{
    int
        i,
        size;
    double
        d,
        nabla_coef,
        sin_th,
        cot_th,
        azi_num_sq;
    Rarray
        theta;
    complex
        idt,
        second_der,
        first_der;

    idt = - I * EQ->dt;
    d = EQ->dtheta;
    size = EQ->ntheta;
    theta = EQ->theta;
    nabla_coef = EQ->nabla_coef;
    azi_num_sq = azi_num * azi_num;

    for (i = 2; i < size - 2; i ++)
    {
        sin_th = sin(theta[i]);
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * idt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );
    }

    // Work near the boundary. `azi_num` = 0 implies Neumann boundary with
    // the outsize domain point as the one preeceding \pi or next of 0.
    // `azi_num` != 0 implies zero at the sphere poles, \pi and 0.
    if (azi_num_sq == 0)
    {
        // same as before but with `azi_num` == 0
        i = 1;
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * idt * nabla_coef * (second_der + cot_th * first_der)
        );

        i = size - 2;
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * idt * nabla_coef * (second_der + cot_th * first_der)
        );

        // at the poles `first_der` == 0
        i = 0;
        second_der = (psi[i+1] - 2 * psi[i] + psi[i+1]) / d / d;
        cn_psi[i] = (I * psi[i] + 0.5 * idt * nabla_coef * second_der);

        i = size - 1;
        second_der = (psi[i-1] - 2 * psi[i] + psi[i-1]) / d / d;
        cn_psi[i] = (I * psi[i] + 0.5 * idt * nabla_coef * second_der);
    }
    else
    {
        
        // same as before but with `psi` == 0 at the boundaries
        i = 1;
        sin_th = sin(theta[i]);
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i]) / d / d;
        first_der = 0.5 * (psi[i+1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * idt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );

        i = size - 2;
        sin_th = sin(theta[i]);
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (- 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (- psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * idt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );

        cn_psi[0] = 0.0 + 0.0 * I;
        cn_psi[size-1] = 0.0 + 0.0 * I;
    }
}


void _get_psi_theta(int ntheta, int nphi, int phi_fft_index,
        Carray psi, Carray psi_th)
{
    if (phi_fft_index > nphi - 2)
    {
        printf("\n\nInvalid access separating psi as function of theta\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < ntheta; j++)
    {
        psi_th[j] = psi[j * nphi + phi_fft_index];
    }
}


void _set_psi_theta(int ntheta, int nphi, int phi_fft_index,
        Carray psi_th, Carray psi)
{
    if (phi_fft_index > nphi - 2)
    {
        printf("\n\nInvalid access separating psi as function of theta\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < ntheta; j++)
    {
        psi[j * nphi + phi_fft_index] = psi_th[j];
    }
}


int splitstep_spherical_shell(EqDataPkg EQ, Carray Sa, Carray Sb)
{

    int
        m,
        l,
        j,
        k,
        nphi,
        ntheta,
        grid_points,
        display_info_stride;
    double
        mu_a,
        mu_b,
        energy,
        kin_energy,
        sin_sq,
        azi_sq,
        dphi,
        dtheta,
        nabla_coef,
        frac_a,
        frac_b,
        ga,
        gb,
        gab,
        Idt,
        norm,
        den_overlap;
    double complex
        dt;
    Cmatrix
        mid;
    Carray
        inter_evol_op,
        phi_fft,
        upper,
        lower,
        upper_m0,
        lower_m0,
        psi_th,
        cn_psi_th,
        cn_solution;
    Rarray
        azi,
        theta,
        abs_square_a,
        abs_square_b,
        inter_pot;

    dt  = - I * EQ->dt;
    Idt = - EQ->dt;

    if (EQ->nt > 1000) display_info_stride = (EQ->nt / 1000);
    else               display_info_stride = 1;

    norm = 1.0;
    nphi = EQ->nphi;
    ntheta = EQ->ntheta;
    grid_points = nphi * ntheta;
    dphi = EQ->dphi;
    dtheta = EQ->dtheta;
    nabla_coef = EQ->nabla_coef;
    frac_a = EQ->frac_a;
    frac_b = EQ->frac_b;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;
    theta = EQ->theta;

    abs_square_a = rarrDef(grid_points);    // |Psi_a|^2
    abs_square_b = rarrDef(grid_points);    // |Psi_b|^2
    inter_pot = rarrDef(grid_points);       // interaction potential
    inter_evol_op = carrDef(grid_points);   // exp[-i * dt/2 * inter_pot]
    phi_fft = carrDef(grid_points);         // FFT with respect to phi
    azi = rarrDef(nphi - 1);                // FFT frequencies/azimuthal number
    psi_th = carrDef(ntheta);               // fixed phi value
    cn_psi_th = carrDef(ntheta);            // Crank-Nicolson explict part
    cn_solution = carrDef(ntheta);          // solution of tridiagonal system

    // Arrays in tridiagonal system from Crank-Nicolson method for theta
    // As for m = 0 we have different boundary conditions, the upper and
    // lower arrays of tridiagonal matrix change
    upper = carrDef(ntheta - 3);
    lower = carrDef(ntheta - 3);
    upper_m0 = carrDef(ntheta - 1);
    lower_m0 = carrDef(ntheta - 1);
    // the diagonal change for all azimuthal numbers
    mid = cmatDef(nphi - 1, ntheta);

    // setup descriptor - MKL implementation of FFT
    // `l`enght of FFT vector must ignore phi = 2 * PI (boundary)
    l = nphi - 1;
    DFTI_DESCRIPTOR_HANDLE desc;
    m = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, l);
    m = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(l));
    m = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(l));
    m = DftiCommitDescriptor(desc);

    // Fourier Frequencies = azimutal quantum numbers = integers
    // the first frequency/azimuthal number is aways zero
    for (m = 0; m < l; m++)
    {
        if (m <= (l - 1) / 2) { azi[m] = m; }
        else                  { azi[m] = m - l; }
    }

    // Crank-Nicolson implicit scheme for theta in phi-momentum space
    // note there are `nphi - 1` azimuthal numbers. m = 0 must be set
    // separately as the boundary condition change.
    for (m = 1; m < nphi - 1; m++)
    {
        azi_sq = azi[m] * azi[m];
        for (j = 0; j < ntheta - 2; j++)
        {
            sin_sq = sin(theta[j + 1]) * sin(theta[j + 1]);
            mid[m][j] = (I + nabla_coef * dt *
                    (1.0 / dtheta / dtheta + 0.5 * azi_sq / sin_sq)
            );
        }
    }
    for (j = 0; j < ntheta; j++)
    {
        mid[0][j] = I + nabla_coef * dt / dtheta / dtheta;
    }

    for (j = 1; j < ntheta - 2; j++)
    {
        upper[j - 1] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                - 0.25 * cos(theta[j]) / sin(theta[j]) / dtheta
        );
        lower[j - 1] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                + 0.25 * cos(theta[j + 1]) / sin(theta[j + 1]) / dtheta
        );
    }

    for (j = 1; j < ntheta - 1; j++)
    {
        upper_m0[j] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                - 0.25 * cos(theta[j]) / sin(theta[j]) / dtheta
        );
        lower_m0[j - 1] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                + 0.25 * cos(theta[j]) / sin(theta[j]) / dtheta
        );
    }
    upper_m0[0] = 2 * (-0.5 * nabla_coef * dt / dtheta / dtheta);
    lower_m0[ntheta - 2] = 2 * (-0.5 * nabla_coef * dt / dtheta / dtheta);

    printf("\n\nProgrs  Energy       Kinect");
    printf("       mu_a         mu_b         overlap");
    sepline();

    // Start time evolution
    for (k = 0; k < EQ->nt; k++)
    {
        carrAbs2(grid_points, Sa, abs_square_a);
        carrAbs2(grid_points, Sb, abs_square_b);

        // EVOLUTION OF SPECIES A

        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (ga * abs_square_a[j]
                    + sqrt(frac_b / frac_a) * gab * abs_square_b[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, Sa, phi_fft);

        // go to momentum space in phi axis. For a fixed theta there is
        // a stride corresponding to phi varying from 0 to 2*PI  though
        // the point 2*PI is ignore due to periodic boundary condition.
        for (j = 0; j < ntheta; j++)
        {
            m = DftiComputeForward(desc, &phi_fft[j*nphi]);
        }

        for (j = 0; j < nphi - 1; j++)
        {
            _get_psi_theta(ntheta, nphi, j, phi_fft, psi_th);
            _explicit_theta(EQ, azi[j], psi_th, cn_psi_th);
            if (j == 0)
            {
                triDiag(ntheta, upper_m0, lower_m0, mid[0],
                        cn_psi_th, cn_solution
                );
            }
            else
            {
                triDiag(ntheta - 2, upper, lower, mid[j],
                        &cn_psi_th[1], &cn_solution[1]
                );
                cn_solution[0] = 0.0 + 0.0 * I;
                cn_solution[ntheta-1] = 0.0 + 0.0 * I;
            }
            _set_psi_theta(ntheta, nphi, j, cn_solution, phi_fft);
        }

        // return to spatial coordinates in phi axis
        for (j = 0; j < ntheta; j++)
        {
            m = DftiComputeBackward(desc, &phi_fft[j*nphi]);
        }

        // back to positional space in phi-axis
        carrCopy(grid_points, phi_fft, Sa);

        // set periodic boundary at 2 * Pi
        for (j = 0; j < ntheta; j++)
        {
            Sa[j * nphi + nphi - 1] = Sa[j * nphi];
        }

        // EVOLUTION OF SPECIES B

        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (gb * abs_square_b[j]
                    + sqrt(frac_a / frac_b) * gab * abs_square_a[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, Sb, phi_fft);

        // go to momentum space in phi axis. For a fixed theta there is
        // a stride corresponding to phi varying from 0 to 2*PI
        for (j = 0; j < ntheta; j++)
        {
            m = DftiComputeForward(desc, &phi_fft[j*nphi]);
        }

        for (j = 0; j < nphi - 1; j++)
        {
            _get_psi_theta(ntheta, nphi, j, phi_fft, psi_th);
            _explicit_theta(EQ, azi[j], psi_th, cn_psi_th);
            if (j == 0)
            {
                triDiag(ntheta, upper_m0, lower_m0, mid[0],
                        cn_psi_th, cn_solution
                );
            }
            else
            {
                triDiag(ntheta - 2, upper, lower, mid[j],
                        &cn_psi_th[1], &cn_solution[1]
                );
                cn_solution[0] = 0.0 + 0.0 * I;
                cn_solution[ntheta-1] = 0.0 + 0.0 * I;
            }
            _set_psi_theta(ntheta, nphi, j, cn_solution, phi_fft);
        }

        // return to spatial coordinates in phi axis
        for (j = 0; j < ntheta; j++)
        {
            m = DftiComputeBackward(desc, &phi_fft[j*nphi]);
        }

        // set periodic boundary at 2 * Pi
        carrCopy(grid_points, phi_fft, Sb);
        for (j = 0; j < ntheta; j++)
        {
            Sb[j * nphi + nphi - 1] = Sb[j * nphi];
        }

        // ANOTHER HALF STEP OF INTERACTION EVOLUTION OPERATOR

        carrAbs2(grid_points, Sa, abs_square_a);
        carrAbs2(grid_points, Sb, abs_square_b);

        // species A
        carrCopy(grid_points, Sa, phi_fft);
        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (ga * abs_square_a[j]
                    + sqrt(frac_b / frac_a) * gab * abs_square_b[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, phi_fft, Sa);

        // species B
        carrCopy(grid_points, Sb, phi_fft);
        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (gb * abs_square_b[j]
                    + sqrt(frac_a / frac_b) * gab * abs_square_a[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, phi_fft, Sb);

        // Renormalize
        renormalize_spheric(EQ, Sa);
        renormalize_spheric(EQ, Sb);

        if ( (k + 1) % (display_info_stride) == 0 )
        {
            carrAbs2(grid_points, Sa, abs_square_a);
            carrAbs2(grid_points, Sb, abs_square_b);
            energy = functionals(EQ, Sa, Sb, &kin_energy, &mu_a, &mu_b);
            den_overlap = density_overlap(EQ, abs_square_a, abs_square_b);
            printf("%5.1lf%%",(100.0 * k) / EQ->nt);
            printf("  %11.8lf  %11.8lf  %11.8lf  %11.8lf  %8.6lf\n",
                    energy, kin_energy, mu_a, mu_b, den_overlap);
        }

    }

    free(inter_evol_op);
    free(phi_fft);
    free(abs_square_a);
    free(abs_square_b);
    free(azi);
    free(upper);
    free(lower);
    free(upper_m0);
    free(lower_m0);
    free(inter_pot);
    free(psi_th);
    free(cn_psi_th);
    free(cn_solution);
    cmatFree(nphi - 1, mid);

    m = DftiFreeDescriptor(&desc);

    return EQ->nt + 1;
}
