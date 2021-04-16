#include "realtimeIntegrator.h"


void _explicit_theta_dyn(
        EqDataPkg EQ, double azi_num, Carray psi, Carray cn_psi
    )
{

/** Apply Crank-Nicolson finite-difference operator in theta direction
    ------------------------------------------------------------------
    After using a Fourier transform in phi axis, the equation in theta
    direction turns out to be dependent on the azimuthal number  which
    also influences the boundary condition at the poles.

    Input Parameters
    ----------------
    `EQ` : physical parameters of GP-equation and grid information
    `azi_num` : azimuthal number from phi-Fourier transformation
    `psi` : sizeof(`EQ->ntheta`) function of theta for a fixed `azi_num`

    Output Parameter
    ----------------
    `cn_psi` : sizeof(`psi`) result of CN operator acting on `psi`

**/

    int
        i,
        size;
    double
        d,
        dt,
        nabla_coef,
        sin_th,
        cot_th,
        azi_num_sq;
    Rarray
        theta;
    double complex
        second_der,
        first_der;

    dt = EQ->dt;
    d = EQ->dtheta;
    size = EQ->ntheta;
    theta = EQ->theta;
    nabla_coef = EQ->nabla_coef;
    azi_num_sq = azi_num * azi_num;

    // set the core points not influenced by the boundaries
    for (i = 2; i < size - 2; i ++)
    {
        sin_th = sin(theta[i]);
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * dt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );
    }

    // Work near the boundary. `azi_num` = 0 implies Neumann boundary condition
    // i.e., the derivative must be zero at theta = 0,\pi. Thus, the forward
    // point at \theta = \pi is just the backward one rotate by \pi in \phi,
    // which turns out to be the simply backward point for `azi_num` = 0.
    // `azi_num` != 0 implies zero at \theta = \pi and 0.
    if (azi_num_sq == 0)
    {
        // same as before but with `azi_num` == 0
        i = 1;
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * dt * nabla_coef * (second_der + cot_th * first_der)
        );

        i = size - 2;
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (psi[i+1] - 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (psi[i+1] - psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * dt * nabla_coef * (second_der + cot_th * first_der)
        );

        // at the poles `first_der` == 0
        i = 0;
        second_der = (psi[i+1] - 2 * psi[i] + psi[i+1]) / d / d;
        cn_psi[i] = (I * psi[i] + 0.5 * dt * nabla_coef * second_der);

        i = size - 1;
        second_der = (psi[i-1] - 2 * psi[i] + psi[i-1]) / d / d;
        cn_psi[i] = (I * psi[i] + 0.5 * dt * nabla_coef * second_der);
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
                + 0.5 * dt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );

        i = size - 2;
        sin_th = sin(theta[i]);
        cot_th = cos(theta[i]) / sin(theta[i]);
        second_der = (- 2 * psi[i] + psi[i-1]) / d / d;
        first_der = 0.5 * (- psi[i-1]) / d;
        cn_psi[i] = (I * psi[i]
                + 0.5 * dt * nabla_coef * (second_der + cot_th * first_der
                - psi[i] * azi_num_sq / sin_th / sin_th)
        );

        cn_psi[0] = 0.0 + 0.0 * I;
        cn_psi[size-1] = 0.0 + 0.0 * I;
    }
}


void _propagate_linear_dyn(
        EqDataPkg EQ, Carray upper_m0, Carray upper,
        Cmatrix l_decomp, Cmatrix u_decomp, Rarray azi,
        DFTI_DESCRIPTOR_HANDLE desc, Carray phi_fft)
{
    int
        j,
        nphi,
        ntheta;

    Carray
        psi_th,
        cn_psi_th,
        cn_solution,
        aux_workspace;

    nphi = EQ->nphi;
    ntheta = EQ->ntheta;

#pragma omp parallel private(j,psi_th,cn_psi_th,cn_solution,aux_workspace)
    {
        psi_th = carrDef(ntheta);
        cn_psi_th = carrDef(ntheta);
        cn_solution = carrDef(ntheta);
        aux_workspace = carrDef(ntheta);

        #pragma omp for
        for (j = 0; j < ntheta; j++)
        {
            DftiComputeForward(desc, &phi_fft[j*nphi]);
        }

        // Solve for azimuthal number == 0

        #pragma omp single
        {
            _get_psi_theta(ntheta, nphi, 0, phi_fft, psi_th);
            _explicit_theta_dyn(EQ, azi[0], psi_th, cn_psi_th);
            tridiag_lu(
                    ntheta, upper_m0, l_decomp[0], u_decomp[0],
                    aux_workspace, cn_psi_th, cn_solution
            );
            _set_psi_theta(ntheta, nphi, 0, cn_solution, phi_fft);
        }

        // Solve for azimuthal number != 0

        #pragma omp for
        for (j = 1; j < nphi - 1; j++)
        {
            _get_psi_theta(ntheta, nphi, j, phi_fft, psi_th);
            _explicit_theta_dyn(EQ, azi[j], psi_th, cn_psi_th);
            tridiag_lu(
                    ntheta - 2, upper, l_decomp[j], u_decomp[j],
                    aux_workspace, &cn_psi_th[1], &cn_solution[1]);
            cn_solution[0] = 0.0 + 0.0 * I;
            cn_solution[ntheta-1] = 0.0 + 0.0 * I;
            _set_psi_theta(ntheta, nphi, j, cn_solution, phi_fft);
        }

        // return to spatial coordinates in phi axis
        #pragma omp for
        for (j = 0; j < ntheta; j++)
        {
            DftiComputeBackward(desc, &phi_fft[j*nphi]);
        }

        free(psi_th);
        free(cn_psi_th);
        free(cn_solution);
        free(aux_workspace);
    }

    // set periodic boundary
    for (j = 0; j < ntheta; j++)
    {
        phi_fft[j * nphi + nphi - 1] = phi_fft[j * nphi];
    }
}


int real_time_evolution(
        EqDataPkg EQ, Carray Sa, Carray Sb, char prefix [], int n_rec)
{

    char
        fname_a[100],
        fname_b[100],
        fname_obs[100];
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
        norm_a,
        norm_b,
        dt,
        lza,
        lzb,
        mu_a,
        mu_b,
        energy,
        kin_energy,
        sin_sq,
        azi_sq,
        dphi,
        dtheta,
        nabla_coef,
        omega,
        frac_a,
        frac_b,
        ga,
        gb,
        gab,
        den_overlap;
    double complex
        Idt;
    Cmatrix
        l_decomp,
        u_decomp,
        mid;
    Carray
        inter_evol_op,
        phi_fft,
        upper,
        lower,
        upper_m0,
        lower_m0;
    Rarray
        azi,
        theta,
        abs_square_a,
        abs_square_b,
        inter_pot;
    FILE
        * file_ptr_a,
        * file_ptr_b,
        * file_obs,
        * file_overlap;

    dt = EQ->dt;
    Idt = - I * dt;
    nphi = EQ->nphi;
    ntheta = EQ->ntheta;
    grid_points = nphi * ntheta;

    strcpy(fname_a, "output/");
    strcat(fname_a, prefix);
    strcat(fname_a, "_obs_realtime.dat");
    file_overlap = fopen(fname_a, "w");
    if (file_overlap == NULL)
    {
        printf("\n\nImpossible to open file %s\n\n", fname_a);
        exit(EXIT_FAILURE);
    }

    strcpy(fname_a, "output/");
    strcat(fname_a, prefix);
    strcat(fname_a, "_speciesA_realtime.dat");
    file_ptr_a = fopen(fname_a, "w");
    if (file_ptr_a == NULL)
    {
        printf("\n\nImpossible to open file %s\n\n", fname_a);
        exit(EXIT_FAILURE);
    }
    carr_inline(file_ptr_a, grid_points, Sa);
    fclose(file_ptr_a);

    strcpy(fname_b, "output/");
    strcat(fname_b, prefix);
    strcat(fname_b, "_speciesB_realtime.dat");
    file_ptr_b = fopen(fname_b, "w");
    if (file_ptr_b == NULL)
    {
        printf("\n\nImpossible to open file %s\n\n", fname_b);
        exit(EXIT_FAILURE);
    }
    carr_inline(file_ptr_b, grid_points, Sb);
    fclose(file_ptr_b);

    strcpy(fname_obs, "output/");
    strcat(fname_obs, prefix);
    strcat(fname_obs, "_log_realtime.dat");

    if (EQ->nt > 1000) display_info_stride = (EQ->nt / 1000);
    else               display_info_stride = 1;

    dphi = EQ->dphi;
    dtheta = EQ->dtheta;
    nabla_coef = EQ->nabla_coef;
    omega = EQ->omega;
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

    // Arrays in tridiagonal system from Crank-Nicolson method for theta
    // As for m = 0 we have different boundary conditions, the upper and
    // lower arrays of tridiagonal matrix change in size and last elements
    upper = carrDef(ntheta - 3);
    lower = carrDef(ntheta - 3);
    upper_m0 = carrDef(ntheta - 1);
    lower_m0 = carrDef(ntheta - 1);

    // the diagonal change for all azimuthal numbers
    mid = cmatDef(nphi - 1, ntheta);

    // Corresponding LU decomposition
    l_decomp = cmatDef(nphi - 1, ntheta);
    u_decomp = cmatDef(nphi - 1, ntheta);

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
        if (m <= (l - 1) / 2) azi[m] = m;
        else                  azi[m] = m - l;
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
                    + dt * omega * azi[m]
            );
        }
    }
    for (j = 0; j < ntheta; j++)
    {
        mid[0][j] = I + nabla_coef * dt / dtheta / dtheta;
    }

    // set off-diagonals for general m != 0
    for (j = 1; j < ntheta - 2; j++)
    {
        upper[j - 1] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                - 0.25 * cos(theta[j]) / sin(theta[j]) / dtheta
        );
        lower[j - 1] = nabla_coef * dt * (-0.5 / dtheta / dtheta
                + 0.25 * cos(theta[j + 1]) / sin(theta[j + 1]) / dtheta
        );
    }

    // set off-diagonals for m = 0
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

    // compute LU-decomposition for all azimuthal numbers
    lu_decomposition(
            ntheta, upper_m0, lower_m0, mid[0], l_decomp[0], u_decomp[0]
    );
    for (j = 1; j < nphi - 1; j++)
    {
        lu_decomposition(
                ntheta - 2, upper, lower, mid[j], l_decomp[j], u_decomp[j]
        );
    }

    printf("\n\nProgrs  Energy       Kinect");
    printf("       mu_a         mu_b         overlap");
    sepline();
    carrAbs2(grid_points, Sa, abs_square_a);
    carrAbs2(grid_points, Sb, abs_square_b);
    norm_a = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_a, dphi);
    norm_b = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_b, dphi);
    energy = functionals(EQ, Sa, Sb, &kin_energy, &mu_a, &mu_b);
    den_overlap = density_overlap(EQ, abs_square_a, abs_square_b);
    printf("%5.1lf%%", 0.0);
    lza = angular_momentum_lz(nphi, ntheta, dphi, theta, Sa);
    lzb = angular_momentum_lz(nphi, ntheta, dphi, theta, Sb);
    printf("  %11.8lf  %11.8lf  %11.8lf  %11.8lf  %8.6lf",
            energy, kin_energy, norm_a, norm_b, den_overlap);
    printf("  %11.8lf  %11.8lf\n", lza, lzb);
    file_obs = fopen(fname_obs, "w");
    if (file_obs == NULL)
    {
        printf("\n\nImpossible to open file %s\n\n", fname_obs);
        exit(EXIT_FAILURE);
    }
    fprintf(file_obs,
            "%.6lf %.6lf %.6lf %.6lf\n",
            0.0, energy, norm_a, norm_b
    );
    fclose(file_obs);
    fprintf(file_overlap,
            "%.10lf %.10lf %.10lf %.10lf\n",
            0.0, den_overlap, lza, lzb
    );

    // Start time evolution
    for (k = 0; k < EQ->nt; k++)
    {
        carrAbs2(grid_points, Sa, abs_square_a);
        carrAbs2(grid_points, Sb, abs_square_b);

        // A) interaction part half time step

        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (ga * abs_square_a[j]
                    + sqrt(frac_b / frac_a) * gab * abs_square_b[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, Sa, phi_fft);

        // A) linear part entire time step

        _propagate_linear_dyn(
                EQ, upper_m0, upper, l_decomp, u_decomp, azi, desc, phi_fft
        );

        // A) transfer data
        carrCopy(grid_points, phi_fft, Sa);

        // B) interaction part half time step

        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (gb * abs_square_b[j]
                    + sqrt(frac_a / frac_b) * gab * abs_square_a[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, Sb, phi_fft);

        // B) linear part entire time step

        _propagate_linear_dyn(
                EQ, upper_m0, upper, l_decomp, u_decomp, azi, desc, phi_fft
        );

        // B) transfer data

        carrCopy(grid_points, phi_fft, Sb);

        // ANOTHER HALF STEP OF INTERACTION EVOLUTION OPERATOR
        // Note that the interaction potential will be updated

        carrAbs2(grid_points, Sa, abs_square_a);
        carrAbs2(grid_points, Sb, abs_square_b);

        // species A - phi_fft as workspace array

        carrCopy(grid_points, Sa, phi_fft);
        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (ga * abs_square_a[j]
                    + sqrt(frac_b / frac_a) * gab * abs_square_b[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, phi_fft, Sa);

        // species B - phi_fft as workspace array

        carrCopy(grid_points, Sb, phi_fft);
        for (j = 0; j < grid_points; j++)
        {
            inter_pot[j] = (gb * abs_square_b[j]
                    + sqrt(frac_a / frac_b) * gab * abs_square_a[j]
            );
        }
        rcarrExp(grid_points, 0.5 * Idt, inter_pot, inter_evol_op);
        carrMultiply(grid_points, inter_evol_op, phi_fft, Sb);

        if ( (k + 1) % (display_info_stride) == 0 )
        {
            carrAbs2(grid_points, Sa, abs_square_a);
            carrAbs2(grid_points, Sb, abs_square_b);
            norm_a = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_a, dphi);
            norm_b = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_b, dphi);
            energy = functionals(EQ, Sa, Sb, &kin_energy, &mu_a, &mu_b);
            den_overlap = density_overlap(EQ, abs_square_a, abs_square_b);
            printf("%5.1lf%%",(100.0 * k) / EQ->nt);
            lza = angular_momentum_lz(nphi, ntheta, dphi, theta, Sa);
            lzb = angular_momentum_lz(nphi, ntheta, dphi, theta, Sb);
            printf("  %11.8lf  %11.8lf  %11.8lf  %11.8lf  %8.6lf",
                    energy, kin_energy, norm_a, norm_b, den_overlap);
            printf("  %11.8lf  %11.8lf\n", lza, lzb);
        }

        if ((k + 1) % n_rec == 0)
        {
            carrAbs2(grid_points, Sa, abs_square_a);
            carrAbs2(grid_points, Sb, abs_square_b);
            norm_a = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_a, dphi);
            norm_b = Rsimps2D_sphere(nphi, ntheta, theta, abs_square_b, dphi);
            energy = functionals(EQ, Sa, Sb, &kin_energy, &mu_a, &mu_b);
            file_obs = fopen(fname_obs, "a");
            if (file_obs == NULL)
            {
                printf("\n\nImpossible to open file %s\n\n", fname_obs);
                exit(EXIT_FAILURE);
            }
            fprintf(file_obs,
                    "%.6lf %.6lf %.6lf %.6lf\n",
                    0.0, energy, norm_a, norm_b
            );
            fclose(file_obs);
            file_ptr_a = fopen(fname_a, "a");
            if (file_ptr_a == NULL)
            {
                printf("\n\nImpossible to open file %s\n\n", fname_a);
                exit(EXIT_FAILURE);
            }
            carr_inline(file_ptr_a, grid_points, Sa);
            fclose(file_ptr_a);
            file_ptr_b = fopen(fname_b, "a");
            if (file_ptr_b == NULL)
            {
                printf("\n\nImpossible to open file %s\n\n", fname_b);
                exit(EXIT_FAILURE);
            }
            carr_inline(file_ptr_b, grid_points, Sb);
            fclose(file_ptr_b);
        }

        carrAbs2(grid_points, Sa, abs_square_a);
        carrAbs2(grid_points, Sb, abs_square_b);
        den_overlap = density_overlap(EQ, abs_square_a, abs_square_b);
        lza = angular_momentum_lz(nphi, ntheta, dphi, theta, Sa);
        lzb = angular_momentum_lz(nphi, ntheta, dphi, theta, Sb);
        fprintf(file_overlap, "%.10lf %.10lf %.10lf %.10lf\n",
                (k + 1) * dt, den_overlap, lza, lzb
        );
        if (fabs(kin_energy) > 1.0 / dtheta)
        {
            printf("WARNING : possible collapse. Abort...");
            break;
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
    cmatFree(nphi - 1, mid);
    cmatFree(nphi - 1, l_decomp);
    cmatFree(nphi - 1, u_decomp);
    fclose(file_overlap);

    m = DftiFreeDescriptor(&desc);

    return k;
}
