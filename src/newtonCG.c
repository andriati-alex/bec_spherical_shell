#include "newtoncg.h"


void set_zero_poles_cg(TwoSpeciesState S)
{
    int
        i,
        j,
        grid_pt;

    j = 0;
    for (i = 0; i < S->nphi; i++)
    {
        grid_pt = j * S->nphi + i;
        S->speca[grid_pt] = 0.0;
        S->specb[grid_pt] = 0.0;
        S->speca_re[grid_pt] = 0.0;
        S->speca_im[grid_pt] = 0.0;
        S->specb_re[grid_pt] = 0.0;
        S->specb_im[grid_pt] = 0.0;
    }
    j = S->ntht - 1;
    for (i = 0; i < S->nphi; i++)
    {
        grid_pt = j * S->nphi + i;
        S->speca[grid_pt] = 0.0;
        S->specb[grid_pt] = 0.0;
        S->speca_re[grid_pt] = 0.0;
        S->speca_im[grid_pt] = 0.0;
        S->specb_re[grid_pt] = 0.0;
        S->specb_im[grid_pt] = 0.0;
    }
    for (j = 0; j < S->ntht; j++)
    {
        grid_pt = j * S->nphi + (S->nphi - 1);
        S->speca[grid_pt] = 0.0;
        S->specb[grid_pt] = 0.0;
        S->speca_re[grid_pt] = 0.0;
        S->speca_im[grid_pt] = 0.0;
        S->specb_re[grid_pt] = 0.0;
        S->specb_im[grid_pt] = 0.0;
    }
}



void grid_residue(
        EqDataPkg EQ,
        Carray state_a,
        Carray state_b,
        Carray grid_res_a,
        Carray grid_res_b,
        double mu_a,
        double mu_b)
{
    int
        j,
        i,
        nphi,
        ntht,
        grid_pt;
    double
        nabla_coef,
        ga,
        gb,
        gab;
    Rarray
        abs_square_a,
        abs_square_b;
    Carray
        laplace_a,
        laplace_b;

    // unpack equation data
    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;

    laplace_a = carrDef(nphi * ntht);
    laplace_b = carrDef(nphi * ntht);
    abs_square_a = rarrDef(nphi * ntht);
    abs_square_b = rarrDef(nphi * ntht);

    carrAbs2(nphi * ntht, state_a, abs_square_a);
    carrAbs2(nphi * ntht, state_b, abs_square_b);

    laplace_app(EQ, state_a, state_b, laplace_a, laplace_b);

    // Compute grid residue for species A

    for (j = 0; j < ntht; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res_a[grid_pt] = (
                    nabla_coef * laplace_a[grid_pt] +
                    state_a[grid_pt] * (
                    ga * abs_square_a[grid_pt] + gab * abs_square_b[grid_pt]
                    )
                    - mu_a * state_a[grid_pt]
            );
        }
    }

    // Compute grid residue for species B

    for (j = 0; j < ntht; j++)
    {
        for (i = 0; i < nphi; i++)
        {
            grid_pt = j * nphi + i;
            grid_res_b[grid_pt] = (
                    nabla_coef * laplace_b[grid_pt] +
                    state_b[grid_pt] * (
                    gb * abs_square_b[grid_pt] + gab * abs_square_a[grid_pt]
                    )
                    - mu_b * state_b[grid_pt]
            );
        }
    }

    free(laplace_a);
    free(laplace_b);
    free(abs_square_a);
    free(abs_square_b);
}



double maxNorm(int N, Carray S)
{

    int
        i;
    double
        maxRes;

    maxRes = 0;

    for (i = 0; i < N; i++)
    {
        if (sqrt(cabs(S[i])) > maxRes)
        {
            maxRes = sqrt(cabs(S[i]));
        }
    }

    return maxRes;
}



void phi_twice_der_cg(
        int nphi, int ntht, double dphi, Carray state, Carray der)
{
    int
        i,
        j,
        grid_pt;

    for (j = 1; j < ntht - 1; j++)
    {
        for (i = 1; i < nphi - 2; i++)
        {
            grid_pt = j * nphi + i;
            der[grid_pt] = (
                    state[grid_pt + 1] -
                    2 * state[grid_pt] +
                    state[grid_pt - 1]
            ) / dphi / dphi;
        }
        grid_pt = j * nphi + 0;
        der[grid_pt] = (
                state[grid_pt + 1] -
                2 * state[grid_pt] +
                state[grid_pt + nphi - 2]
        ) / dphi / dphi;
        grid_pt = j * nphi + nphi - 2;
        der[grid_pt] = (
                state[grid_pt - (nphi - 2)] -
                2 * state[grid_pt] +
                state[grid_pt - 1]
        ) / dphi / dphi;
    }
}



void tht_twice_der_cg(
        int nphi, int ntht, double dtht, Carray state, Carray der)
{
    int
        i,
        j,
        grid_pt;

    for (j = 2; j < ntht - 2; j++)
    {
        for (i = 0; i < nphi - 1; i++)
        {
            grid_pt = j * nphi + i;
            der[grid_pt] = (
                    state[grid_pt + nphi] -
                    2 * state[grid_pt] +
                    state[grid_pt - nphi]
            ) / dtht / dtht;
        }
    }

    j = 1;
    for (i = 0; i < nphi - 1; i++)
    {
        grid_pt = j * nphi + i;
        der[grid_pt] = (
                state[grid_pt + nphi] - 2 * state[grid_pt]
        ) / dtht / dtht;
    }

    j = ntht - 2;
    for (i = 0; i < nphi - 1; i++)
    {
        grid_pt = j * nphi + i;
        der[grid_pt] = (
                -2 * state[grid_pt] + state[grid_pt - nphi]
        ) / dtht / dtht;
    }
}



void tht_der_cg(
        int nphi, int ntht, double dtht, Carray state, Carray der)
{
    int
        i,
        j,
        grid_pt;

    for (j = 2; j < ntht - 2; j++)
    {
        for (i = 0; i < nphi - 1; i++)
        {
            grid_pt = j * nphi + i;
            der[grid_pt] = 0.5 * (
                    state[grid_pt + nphi] - state[grid_pt - nphi]
            ) / dtht;
        }
    }

    j = 1;
    for (i = 0; i < nphi - 1; i++)
    {
        grid_pt = j * nphi + i;
        der[grid_pt] = 0.5 * state[grid_pt + nphi] / dtht;
    }

    j = ntht - 2;
    for (i = 0; i < nphi - 1; i++)
    {
        grid_pt = j * nphi + i;
        der[grid_pt] = - 0.5 * state[grid_pt - nphi] / dtht;
    }
}



void laplace_cg(
        EqDataPkg EQ,
        Carray state,
        Carray laplace)
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

    tht_twice_der_cg(nphi, ntht, dtht, state, der2tht);
    tht_der_cg(nphi, ntht, dtht, state, der_tht);
    phi_twice_der_cg(nphi, ntht, dphi, state, der2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi - 1; i++)
        {
            grid_pt = j * nphi + i;
            laplace[grid_pt] = (
                    der2tht[grid_pt] +
                    der_tht[grid_pt] * cos_tht / sin_tht +
                    der2phi[grid_pt] / sin_tht / sin_tht
            );
        }
    }

    free(der2phi);
    free(der2tht);
    free(der_tht);
}



void laplace_cg_dagger(
        EqDataPkg EQ,
        Carray state,
        Carray laplace)
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

    tht_twice_der_cg(nphi, ntht, dtht, state, der2tht);
    tht_der_cg(nphi, ntht, dtht, state, der_tht);
    phi_twice_der_cg(nphi, ntht, dphi, state, der2phi);

    for (j = 1; j < ntht - 1; j++)
    {
        sin_tht = sin(tht[j]);
        cos_tht = cos(tht[j]);
        for (i = 0; i < nphi - 1; i++)
        {
            grid_pt = j * nphi + i;
            laplace[grid_pt] = (
                    der2tht[grid_pt] -
                    der_tht[grid_pt] * cos_tht / sin_tht +
                    der2phi[grid_pt] / sin_tht / sin_tht
            );
        }
    }

    free(der2phi);
    free(der2tht);
    free(der_tht);
}



void linearized_op(
        EqDataPkg EQ, TwoSpeciesState newton, TwoSpeciesState conj_grad,
        TwoSpeciesState lin_op, double mu_a, double mu_b, int dagger)
{

    int
        i,
        j,
        nphi,
        ntht,
        grid_pt;
    double
        nabla_coef,
        ga,
        gb,
        gab,
        newton_a_re,
        newton_a_im,
        newton_b_re,
        newton_b_im,
        newton_a_re_sq,
        newton_a_im_sq,
        newton_b_re_sq,
        newton_b_im_sq,
        conj_grad_a_re,
        conj_grad_a_im,
        conj_grad_b_re,
        conj_grad_b_im;
    Carray
        laplace_a,
        laplace_b;

    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;
    nphi = EQ->nphi;
    ntht = EQ->ntheta;

    laplace_a = carrDef(nphi * ntht);
    laplace_b = carrDef(nphi * ntht);

    if (dagger)
    {
        laplace_cg_dagger(EQ, conj_grad->speca, laplace_a);
        laplace_cg_dagger(EQ, conj_grad->specb, laplace_b);
    }
    else
    {
        laplace_cg(EQ, conj_grad->speca, laplace_a);
        laplace_cg(EQ, conj_grad->specb, laplace_b);
    }

    set_zero_poles_cg(lin_op);

    for (j = 1; j < ntht - 1; j++)
    {
        for (i = 0; i < nphi - 1; i++)
        {
            grid_pt = j * ntht + i;
            conj_grad_a_re = conj_grad->speca_re[grid_pt];
            conj_grad_a_im = conj_grad->speca_im[grid_pt];
            conj_grad_b_re = conj_grad->specb_re[grid_pt];
            conj_grad_b_im = conj_grad->specb_im[grid_pt];
            newton_a_re = newton->speca_re[grid_pt];
            newton_a_im = newton->speca_im[grid_pt];
            newton_b_re = newton->specb_re[grid_pt];
            newton_b_im = newton->specb_im[grid_pt];
            newton_a_re_sq = newton_a_re * newton_a_re;
            newton_b_re_sq = newton_b_re * newton_b_re;
            newton_a_im_sq = newton_a_im * newton_a_im;
            newton_b_im_sq = newton_b_im * newton_b_im;
            lin_op->speca_re[grid_pt] = (
                    nabla_coef * creal(laplace_a[grid_pt]) +
                    conj_grad_a_re * (
                            - mu_a
                            + ga * (3 * newton_a_re_sq + newton_a_im_sq)
                            + gab * (newton_b_re_sq + newton_b_im_sq)
                    ) +
                    2 * ga * newton_a_im * newton_a_re * conj_grad_a_im +
                    2 * gab * newton_b_re * newton_a_re * conj_grad_b_re +
                    2 * gab * newton_b_im * newton_a_re * conj_grad_b_im
            );
            lin_op->speca_im[grid_pt] = (
                    nabla_coef * cimag(laplace_a[grid_pt]) +
                    conj_grad_a_im * (
                            - mu_a
                            + ga * (newton_a_re_sq + 3 * newton_a_im_sq)
                            + gab * (newton_b_re_sq + newton_b_im_sq)
                    ) +
                    2 * ga * newton_a_re * newton_a_im * conj_grad_a_re +
                    2 * gab * newton_b_re * newton_a_im * conj_grad_b_re +
                    2 * gab * newton_b_im * newton_a_im * conj_grad_b_im
            );
            lin_op->specb_re[grid_pt] = (
                    nabla_coef * creal(laplace_b[grid_pt]) +
                    conj_grad_b_re * (
                            - mu_b
                            + gb * (3 * newton_b_re_sq + newton_b_im_sq)
                            + gab * (newton_a_re_sq + newton_a_im_sq)
                    ) +
                    2 * gb * newton_b_im * newton_b_re * conj_grad_b_im +
                    2 * gab * newton_a_re * newton_b_re * conj_grad_a_re +
                    2 * gab * newton_a_im * newton_b_re * conj_grad_a_im
            );
            lin_op->specb_im[grid_pt] = (
                    nabla_coef * cimag(laplace_b[grid_pt]) +
                    conj_grad_b_im * (
                            - mu_b
                            + gb * (newton_b_re_sq + 3 * newton_b_im_sq)
                            + gab * (newton_a_re_sq + newton_a_im_sq)
                    ) +
                    2 * gb * newton_b_re * newton_b_im * conj_grad_b_re +
                    2 * gab * newton_a_re * newton_b_im * conj_grad_a_re +
                    2 * gab * newton_a_im * newton_b_im * conj_grad_a_im
            );
        }
    }

    set_states_from_parts(lin_op);
    free(laplace_a);
    free(laplace_b);
}



int conjgrad(
        EqDataPkg EQ, double mu_a, double mu_b, TwoSpeciesState newton,
        TwoSpeciesState newton_res, TwoSpeciesState conj_grad, double tol)
{

/** CONJUGATE GRADIENT ITERATIVE METHOD FOR 2D-BEC
  * **********************************************
  *
  * Since the Newton method is dereived from a complex equation the system
  * can be splited in two parts, here we chose in real and imaginaty parts
  *
  * The differential operator that is discretized in the given gird assume
  * a 2x2 block form if organized in a matrix. The final vector space size
  * is 2 * nx * ny, where nx and ny are the number of points in each direc
  * tion in the discretized domain.
  *
  * Therefore any time the vector is updated, we call twice  the underlying
  * routine for the first nx*ny elements corresponding to real part and and
  * the later nx*ny related to imaginary part
***************************************************************************/

    int
        N,
        l,
        nphi,
        ntht,
        maxiter;
    double
        a,
        beta,
        err;
    TwoSpeciesState
        res,
        dir,
        prev_res,
        lin_op;

    l = 0;
    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    N = nphi * ntht;
    maxiter = 3 * N;
    res = alloc_two_species_struct(nphi, ntht);
    // res_tilde = alloc_two_species_struct(nphi, ntht);
    dir = alloc_two_species_struct(nphi, ntht);
    // dir_tilde = alloc_two_species_struct(nphi, ntht);
    prev_res = alloc_two_species_struct(nphi, ntht);
    // prev_res_tilde = alloc_two_species_struct(nphi, ntht);
    lin_op = alloc_two_species_struct(nphi, ntht);
    // lin_op_dag = alloc_two_species_struct(nphi, ntht);

    // linearizedOp(EQ,phir,phii,mu,fr,fi,lfr,lfi);
    linearized_op(EQ, newton, conj_grad, lin_op, mu_a, mu_b, 0);

    // compute residue
    // rarrSub(nx*ny,Fr,lfr,realRes);
    // rarrSub(nx*ny,Fi,lfi,imagRes);
    rarrSub(nphi * ntht, newton_res->speca_re, lin_op->speca_re, res->speca_re);
    rarrSub(nphi * ntht, newton_res->speca_im, lin_op->speca_im, res->speca_im);
    rarrSub(nphi * ntht, newton_res->specb_re, lin_op->specb_re, res->specb_re);
    rarrSub(nphi * ntht, newton_res->specb_im, lin_op->specb_im, res->specb_im);
    set_zero_poles_cg(res);
    set_states_from_parts(res);

    // pkg_states(res->speca, res->specb, res_tilde);
    // pkg_states(res->speca, res->specb, dir_tilde);

    pkg_states(res->speca, res->specb, dir);

    err = (
            rarrMod(nphi * ntht, res->speca_re) +
            rarrMod(nphi * ntht, res->speca_im) +
            rarrMod(nphi * ntht, res->specb_re) +
            rarrMod(nphi * ntht, res->specb_im)
    );

    err = maxNorm(N, res->speca) + maxNorm(N, res->specb);
    printf("\t\t%.8lf", err);

    while (err > tol)
    {
        set_states_from_parts(res);
        set_states_from_parts(dir);
        // set_states_from_parts(res_tilde);
        // set_states_from_parts(dir_tilde);

        // hold the result of operator acting on direction
        linearized_op(EQ, newton, dir, lin_op, mu_a, mu_b, 0);
        // linearized_op(EQ, newton, dir_tilde, lin_op_dag, mu_a, mu_b, 1);

        // scalar to update solution and residue
        a = (
                (
                    rarrDot(N, res->speca_re, res->speca_re) +
                    rarrDot(N, res->speca_im, res->speca_im) +
                    rarrDot(N, res->specb_re, res->specb_re) +
                    rarrDot(N, res->specb_im, res->specb_im)
                ) /
                (
                    rarrDot(N, dir->speca_re, lin_op->speca_re) +
                    rarrDot(N, dir->speca_im, lin_op->speca_im) +
                    rarrDot(N, dir->specb_re, lin_op->specb_re) +
                    rarrDot(N, dir->specb_im, lin_op->specb_im)
                )
        );

        // record residue from this iteration to safely update
        pkg_states(res->speca, res->specb, prev_res);

        //pkg_states(res_tilde->speca, res_tilde->specb, prev_res_tilde);

        // update residue
        rarrUpdate(
                N, prev_res->speca_re, - a, lin_op->speca_re, res->speca_re
        );
        rarrUpdate(
                N, prev_res->speca_im, - a, lin_op->speca_im, res->speca_im
        );
        rarrUpdate(
                N, prev_res->specb_re, - a, lin_op->specb_re, res->specb_re
        );
        rarrUpdate(
                N, prev_res->specb_im, - a, lin_op->specb_im, res->specb_im
        );
        set_states_from_parts(res);

        rarrUpdate(
                N, conj_grad->speca_re, a, dir->speca_re, conj_grad->speca_re
        );
        rarrUpdate(
                N, conj_grad->speca_im, a, dir->speca_im, conj_grad->speca_im
        );
        rarrUpdate(
                N, conj_grad->specb_re, a, dir->specb_re, conj_grad->specb_re
        );
        rarrUpdate(
                N, conj_grad->specb_im, a, dir->specb_im, conj_grad->specb_im
        );
        set_states_from_parts(conj_grad);

        // scalar to update direction
        beta = (
                (
                    rarrDot(N, res->speca_re, res->speca_re) +
                    rarrDot(N, res->speca_im, res->speca_im) +
                    rarrDot(N, res->specb_re, res->specb_re) +
                    rarrDot(N, res->specb_im, res->specb_im)
                ) /
                (
                    rarrDot(N, prev_res->speca_re, prev_res->speca_re) +
                    rarrDot(N, prev_res->speca_im, prev_res->speca_im) +
                    rarrDot(N, prev_res->specb_re, prev_res->specb_re) +
                    rarrDot(N, prev_res->specb_im, prev_res->specb_im)
                )
        );

        rarrUpdate(N, res->speca_re, beta, dir->speca_re, dir->speca_re);
        rarrUpdate(N, res->speca_im, beta, dir->speca_im, dir->speca_im);
        rarrUpdate(N, res->specb_re, beta, dir->specb_re, dir->specb_re);
        rarrUpdate(N, res->specb_im, beta, dir->specb_im, dir->specb_im);

        err = (
                rarrMod(nphi * ntht, res->speca_re) +
                rarrMod(nphi * ntht, res->speca_im) +
                rarrMod(nphi * ntht, res->specb_re) +
                rarrMod(nphi * ntht, res->specb_im)
        );

        printf("\n\t%.6lf", err);

        l = l + 1; // Update iteration counter

        if (l == maxiter)
        {
            printf("\n\nWARNING : exit before achieve desired residual ");
            printf("value in Conjugate Gradient method due to max number ");
            printf("of iterations given =  %d\n\n", maxiter);
            break;
        }
    }

    // Free function local memory
    release_two_species_state(res);
    release_two_species_state(dir);
    release_two_species_state(lin_op);
    release_two_species_state(prev_res);
    // release_two_species_state(res_tilde);
    // release_two_species_state(dir_tilde);
    // release_two_species_state(lin_op_dag);
    // release_two_species_state(prev_res_tilde);

    return l;
}



void stationaryNewton(EqDataPkg EQ, Carray Sa, Carray Sb,
        double err_tol, int iter_tol)
{

    int
        i,
        N,
        nphi,
        ntht,
        Niter,
        CGiter;
    double
        E,
        kin,
        mu_a,
        mu_b,
        dphi,
        norm_a,
        norm_b,
        error_a,
        error_b,
        CGtol;
    Carray
        grid_res_a,
        grid_res_b;
    Rarray
        abs_square_a,
        abs_square_b;
    TwoSpeciesState
        conj_grad,
        newton_step,
        newton_res;

    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    N = nphi * ntht;

    newton_res = alloc_two_species_struct(nphi, ntht);
    newton_step = alloc_two_species_struct(nphi, ntht);
    conj_grad = alloc_two_species_struct(nphi, ntht);
    grid_res_a = carrDef(N);
    grid_res_b = carrDef(N);
    abs_square_a = rarrDef(N);
    abs_square_b = rarrDef(N);

    carrAbs2(N, Sa, abs_square_a);
    carrAbs2(N, Sb, abs_square_b);
    norm_a = sqrt(Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_a, dphi));
    norm_b = sqrt(Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_b, dphi));
    E = functionals(EQ, Sa, Sb, &kin, &mu_a, &mu_b);

    // stationaryOp(EQ,mu,fr,fi,Fcg_real,Fcg_imag);
    grid_residue(EQ, Sa, Sb, grid_res_a, grid_res_b, mu_a, mu_b);
    error_a = maxNorm(N, grid_res_a);
    error_b = maxNorm(N, grid_res_b);

    // right hand side of linear operator passed to iteratice CG method
    for (i = 0; i < N; i++)
    {
        grid_res_a[i] = - grid_res_a[i];
        grid_res_b[i] = - grid_res_b[i];
    }
    pkg_states(grid_res_a, grid_res_b, newton_res);

    printf("\nNewton It.      Error    CG It.     ");
    printf("Energy     mu        Norm");
    sepline();

    // NEWTON LOOP
    Niter = 0;
    CGiter = 0;
    while (error_a + error_b > err_tol)
    {
        printf("\n%5d       %10.5lf   %6d   ", Niter, error_a + error_b, CGiter);
        printf("%9.5lf  %9.5lf  %9.6lf   ", E, mu_a, mu_b);
        printf("%9.5lf  %9.6lf\n", norm_a, norm_b);

        // Initial guess for conjugate-gradient method
        carrFill(N, 0.0, conj_grad->speca);
        carrFill(N, 0.0, conj_grad->specb);
        set_states_real_imag(conj_grad);

        pkg_states(Sa, Sb, newton_step);

        // Solve with conj. grad. according to current newton error
        if (1E-4 * (error_a + error_b) > 1E-2)
        {
            CGtol = 1E-2;
        }
        else
        {
            CGtol = 1E-4 * (error_a + error_b);
        }
        // set_zero_poles_cg(newton_res);
        CGiter = conjgrad(
                EQ, mu_a, mu_b, newton_step, newton_res, conj_grad, CGtol
        );

        // update solution from newton iteration
        carrAdd(N, Sa, conj_grad->speca, Sa);
        carrAdd(N, Sb, conj_grad->specb, Sb);

        grid_residue(
                EQ, Sa, Sb, grid_res_a, grid_res_b, mu_a, mu_b
        );
        error_a = maxNorm(N, grid_res_a);
        error_b = maxNorm(N, grid_res_b);

        for (i = 0; i < N; i++)
        {
            grid_res_a[i] = - grid_res_a[i];
            grid_res_b[i] = - grid_res_b[i];
        }
        pkg_states(grid_res_a, grid_res_b, newton_res);

        Niter = Niter + 1;

        if (Niter > iter_tol) break;
    }

    printf("\n\n%5d       %10.5lf   %6d   ",Niter, error_a + error_b, CGiter);
    printf("%9.5lf  %9.5lf  %9.6lf", E, mu_a, mu_b);
    sepline();

    if (Niter > iter_tol)
    {
        printf("\nWARNING : Achieved maximum iterations allowed");
        printf(" by the user in job-ncg.conf file. Exiting\n");
    }

    release_two_species_state(newton_step);
    release_two_species_state(newton_res);
    release_two_species_state(conj_grad);
}
