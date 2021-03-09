#include "newtoncg.h"


void appr_poles(int nphi, int ntht, Carray state)
{
    int
        i,
        grid_pt;
    double complex
        avg_pole,
        first_order,
        second_order;

    avg_pole = 0.0;
    for (i = 0; i < nphi; i++)
    {
        grid_pt = nphi + i;
        first_order = - 0.5 * (
                - 3 * state[grid_pt] +
                4 * state[grid_pt + nphi] -
                state[grid_pt + 2 * nphi]
        );
        second_order = 0.5 * (
                state[grid_pt + 2 * nphi] -
                2 * state[grid_pt + nphi] +
                state[grid_pt]
        );
        avg_pole += (state[grid_pt] + first_order + second_order) / nphi;
    }
    carrFill(nphi, avg_pole, &state[0]);

    avg_pole = 0.0;
    for (i = 0; i < nphi; i++)
    {
        grid_pt = (ntht - 2) * nphi + i;
        first_order = 0.5 * (
                3 * state[grid_pt] -
                4 * state[grid_pt - nphi] +
                state[grid_pt - 2 * nphi]
        );
        second_order = 0.5 * (
                state[grid_pt - 2 * nphi] -
                2 * state[grid_pt - nphi] +
                state[grid_pt]
        );
        avg_pole += (state[grid_pt] + first_order + second_order) / nphi;
    }
    carrFill(nphi, avg_pole, &state[(ntht - 1) * nphi]);
}


double inner_cg(TwoSpeciesState s1, TwoSpeciesState s2, TwoSpeciesState aux,
        Rarray theta, double dphi)
{

/** Inner product of two state of two component system. Use real and imag
    time separately and integrate over the sphere                    **/

    int
        nphi,
        ntht,
        N;

    nphi = s1->nphi;
    ntht = s1->ntht;
    N = nphi * ntht;

    rarrMultiply(N, s1->speca_re, s2->speca_re, aux->speca_re);
    rarrMultiply(N, s1->speca_im, s2->speca_im, aux->speca_im);
    rarrMultiply(N, s1->specb_re, s2->specb_re, aux->specb_re);
    rarrMultiply(N, s1->specb_im, s2->specb_im, aux->specb_im);
    return (
        Rsimps2D_sphere(nphi, ntht, theta, aux->speca_re, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->speca_im, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->specb_re, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->specb_im, dphi)
    );
}


void update_state_cg(
        TwoSpeciesState s1,
        double scalar,
        TwoSpeciesState s2,
        TwoSpeciesState out)
{

    int
        nphi,
        ntht,
        N;

    nphi = s1->nphi;
    ntht = s1->ntht;
    N = nphi * ntht;

    rarrUpdate(N, s1->speca_re, scalar, s2->speca_re, out->speca_re);
    rarrUpdate(N, s1->speca_im, scalar, s2->speca_im, out->speca_im);
    rarrUpdate(N, s1->specb_re, scalar, s2->specb_re, out->specb_re);
    rarrUpdate(N, s1->specb_im, scalar, s2->specb_im, out->specb_im);

    set_states_from_parts(out);
}


double self_inner_cg(TwoSpeciesState s, TwoSpeciesState aux,
        Rarray theta, double dphi)
{
    int
        nphi,
        ntht,
        N;

    nphi = s->nphi;
    ntht = s->ntht;
    N = nphi * ntht;

    rarrMultiply(N, s->speca_re, s->speca_re, aux->speca_re);
    rarrMultiply(N, s->speca_im, s->speca_im, aux->speca_im);
    rarrMultiply(N, s->specb_re, s->specb_re, aux->specb_re);
    rarrMultiply(N, s->specb_im, s->specb_im, aux->specb_im);
    return (
        Rsimps2D_sphere(nphi, ntht, theta, aux->speca_re, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->speca_im, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->specb_re, dphi) +
        Rsimps2D_sphere(nphi, ntht, theta, aux->specb_im, dphi)
    );
}


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
    for (j = 1; j < S->ntht - 1; j++)
    {
        grid_pt = j * S->nphi + (S->nphi - 1);
        S->speca[grid_pt] = S->speca[j * S->nphi];
        S->specb[grid_pt] = S->specb[j * S->nphi];
        S->speca_re[grid_pt] = S->speca_re[j * S->nphi];
        S->speca_im[grid_pt] = S->speca_im[j * S->nphi];
        S->specb_re[grid_pt] = S->specb_re[j * S->nphi];
        S->specb_im[grid_pt] = S->specb_im[j * S->nphi];
    }
}


double poles_std(int nphi, Carray f_pole)
{
    int
        i;
    double complex
        avg;
    double
        std;

    avg = 0.0;
    std = 0.0;
    for (i = 0; i < nphi; i++) avg += f_pole[i] / nphi;
    for (i = 0; i < nphi; i++) std += cabs(avg - f_pole[i]) / nphi;
    return std;
}


double constant_poles(TwoSpeciesState S)
{
    int nphi = S->nphi;
    int j;
    j = S->ntht - 1;
    return (
            poles_std(nphi, S->speca) +
            poles_std(nphi, S->specb) +
            poles_std(nphi, &S->speca[j * nphi]) +
            poles_std(nphi, &S->specb[j * nphi])
    );
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
        TwoSpeciesState lin_op, double mu_a, double mu_b)
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

    laplace_app(EQ, conj_grad->speca, conj_grad->specb, laplace_a, laplace_b);

    for (j = 0; j < ntht; j++)
    {
        for (i = 0; i < nphi; i++)
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
    // set_zero_poles_cg(lin_op);

    free(laplace_a);
    free(laplace_b);
}


int conjgrad_stab(
        EqDataPkg EQ, double mu_a, double mu_b, TwoSpeciesState newton,
        TwoSpeciesState newton_res, TwoSpeciesState conj_grad, double tol,
        unsigned int init_random)
{
    int
        nphi,
        ntht,
        it_counter;
    Rarray
        tht;
    double
        dphi,
        err,
        alpha_scalar,
        omega_scalar,
        beta_scalar;
    TwoSpeciesState
        res,
        dir,
        s_aux,
        res_star,
        prev_res,
        inner,
        lin_op,
        lin_op_s;

    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    tht = EQ->theta;
    res = alloc_two_species_struct(nphi, ntht);
    prev_res = alloc_two_species_struct(nphi, ntht);
    s_aux = alloc_two_species_struct(nphi, ntht);
    dir = alloc_two_species_struct(nphi, ntht);
    lin_op = alloc_two_species_struct(nphi, ntht);
    lin_op_s = alloc_two_species_struct(nphi, ntht);
    inner = alloc_two_species_struct(nphi, ntht);
    res_star = alloc_two_species_struct(nphi, ntht);

    linearized_op(EQ, newton, conj_grad, lin_op, mu_a, mu_b);
    update_state_cg(newton_res, -1.0, lin_op, res);
    pkg_states(res->speca, res->specb, dir);
    pkg_states(res->speca, res->specb, res_star);

    err = sqrt(self_inner_cg(res, inner, tht, dphi));

    if (init_random)
    {
        carrRandom(nphi * ntht, init_random, res_star->speca);
        carrRandom(nphi * ntht, init_random, res_star->specb);
        set_zero_poles_cg(res_star);
        appr_poles(nphi, ntht, res_star->speca);
        appr_poles(nphi, ntht, res_star->specb);
        set_states_real_imag(res_star);
    }

    it_counter = 0;
    printf("\n\titer %5d error : %.7lf", it_counter, err);

    while (err > tol)
    {
        linearized_op(EQ, newton, dir, lin_op, mu_a, mu_b);
        if (inner_cg(res_star, lin_op, inner, tht, dphi) == 0)
        {
            printf("\n\nZero Division in alpha scalar!\n\n");
            exit(EXIT_FAILURE);
        }
        alpha_scalar = (
                inner_cg(res_star, res, inner, tht, dphi) /
                inner_cg(res_star, lin_op, inner, tht, dphi)
        );
        update_state_cg(res, -alpha_scalar, lin_op, s_aux);
        linearized_op(EQ, newton, s_aux, lin_op_s, mu_a, mu_b);
        if (self_inner_cg(lin_op_s, inner, tht, dphi) == 0)
        {
            printf("\n\nZero Division in omega scalar!\n\n");
            exit(EXIT_FAILURE);
        }
        omega_scalar = (
                inner_cg(s_aux, lin_op_s, inner, tht, dphi) /
                self_inner_cg(lin_op_s, inner, tht, dphi)
        );
        // update solution - use `prev_res` as step workspace variable
        update_state_cg(conj_grad, alpha_scalar, dir, prev_res);
        update_state_cg(prev_res, omega_scalar, s_aux, conj_grad);
        if (omega_scalar == 0)
        {
            err = sqrt(self_inner_cg(res, inner, tht, dphi));
            printf("\n\nOmega == 0 found. Exiting with error %.8lf\n", err);
            return it_counter;
        }
        // residual update
        pkg_states(res->speca, res->specb, prev_res);
        update_state_cg(s_aux, -omega_scalar, lin_op_s, res);
        if (inner_cg(res_star, prev_res, inner, tht, dphi) == 0)
        {
            printf("\n\nZero Division in beta scalar!\n\n");
            return it_counter;
            // exit(EXIT_FAILURE);
        }
        beta_scalar = (
                inner_cg(res_star, res, inner, tht, dphi) /
                inner_cg(res_star, prev_res, inner, tht, dphi)
        ) * alpha_scalar / omega_scalar;
        // update direction - use `prev_res` as step workspace variable
        update_state_cg(res, beta_scalar, dir, prev_res);
        update_state_cg(
                prev_res, -beta_scalar * omega_scalar, lin_op, dir
        );
        // update solution error
        err = sqrt(self_inner_cg(res, inner, tht, dphi));

        it_counter++;
        if (it_counter % 20 == 0)
        {
            printf("\n\titer %5d error : %.7lf", it_counter, err);
            linearized_op(EQ, newton, conj_grad, lin_op, mu_a, mu_b);
            update_state_cg(newton_res, -1.0, lin_op, prev_res);
            err = sqrt(self_inner_cg(prev_res, inner, tht, dphi));
            printf("\n\titer %5d error : %.7lf", it_counter, err);
        }
    }

    release_two_species_state(res);
    release_two_species_state(prev_res);
    release_two_species_state(dir);
    release_two_species_state(res_star);
    release_two_species_state(inner);
    release_two_species_state(lin_op);
    release_two_species_state(lin_op_s);
    release_two_species_state(s_aux);
    return it_counter;
}


void stationaryNewton(EqDataPkg EQ, Carray Sa, Carray Sb,
        double err_tol, int iter_tol, double mu_a, double mu_b)
{

    int
        i,
        N,
        nphi,
        ntht,
        Niter,
        CGiter;
    double
        dphi,
        cg_error,
        error_newton;
    Carray
        grid_res_a,
        grid_res_b;
    Rarray
        abs_square_a,
        abs_square_b;
    TwoSpeciesState
        workspace,
        conj_grad,
        newton_step,
        newton_res;

    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    N = nphi * ntht;

    workspace = alloc_two_species_struct(nphi, ntht);
    newton_res = alloc_two_species_struct(nphi, ntht);
    newton_step = alloc_two_species_struct(nphi, ntht);
    conj_grad = alloc_two_species_struct(nphi, ntht);
    grid_res_a = carrDef(N);
    grid_res_b = carrDef(N);
    abs_square_a = rarrDef(N);
    abs_square_b = rarrDef(N);

    grid_residue(EQ, Sa, Sb, grid_res_a, grid_res_b, mu_a, mu_b);

    // right hand side of linear operator passed to iteratice CG method
    for (i = 0; i < N; i++)
    {
        grid_res_a[i] = - grid_res_a[i];
        grid_res_b[i] = - grid_res_b[i];
    }
    pkg_states(grid_res_a, grid_res_b, newton_res);
    error_newton = sqrt(self_inner_cg(newton_res, workspace, EQ->theta, dphi));

    printf("\nNewton It.      Error    CG It.     ");
    printf("Energy     mu        Norm");
    sepline();

    // NEWTON LOOP
    Niter = 0;
    CGiter = 0;
    while (error_newton > err_tol)
    {
        printf("\n%5d       %10.5lf   %6d   ", Niter, error_newton, CGiter);
        printf("%9.5lf  %9.6lf   ", mu_a, mu_b);

        // Initial guess for conjugate-gradient method
        carrFill(N, 0.0, conj_grad->speca);
        carrFill(N, 0.0, conj_grad->specb);
        set_states_real_imag(conj_grad);

        pkg_states(Sa, Sb, newton_step);

        cg_error = 0.25 * error_newton;
        CGiter = conjgrad_stab(
                EQ, mu_a, mu_b, newton_step, newton_res, conj_grad, cg_error, 0
        );

        // update solution from newton iteration
        carrAdd(N, Sa, conj_grad->speca, Sa);
        carrAdd(N, Sb, conj_grad->specb, Sb);
        appr_poles(nphi, ntht, Sa);
        appr_poles(nphi, ntht, Sb);

        grid_residue(
                EQ, Sa, Sb, grid_res_a, grid_res_b, mu_a, mu_b
        );

        for (i = 0; i < N; i++)
        {
            grid_res_a[i] = - grid_res_a[i];
            grid_res_b[i] = - grid_res_b[i];
        }
        pkg_states(grid_res_a, grid_res_b, newton_res);
        error_newton = sqrt(
                self_inner_cg(newton_res, workspace, EQ->theta, dphi)
        );

        Niter = Niter + 1;

        if (Niter > iter_tol) break;
    }

    printf("\n\n%5d       %10.5lf   %6d   ",Niter, error_newton, CGiter);
    printf("%9.5lf  %9.6lf", mu_a, mu_b);
    sepline();

    if (Niter > iter_tol)
    {
        printf("\nWARNING : Achieved maximum iterations allowed");
        printf(" by the user in job-ncg.conf file. Exiting\n");
    }

    release_two_species_state(newton_step);
    release_two_species_state(newton_res);
    release_two_species_state(conj_grad);
    release_two_species_state(workspace);
    free(abs_square_a);
    free(abs_square_b);
    free(grid_res_a);
    free(grid_res_b);
}


void stationaryFixedNorm(EqDataPkg EQ, Carray Sa, Carray Sb,
        double err_tol, int iter_tol)
{

    int
        N,
        nphi,
        ntht;
    double
        E,
        kin,
        aux,
        mu_a,
        mu_b,
        dphi,
        norm_a,
        norm_b,
        old_mu_a,
        old_mu_b,
        old_norm_a,
        old_norm_b,
        der_a,
        der_b;
    Rarray
        abs_square_a,
        abs_square_b;

    nphi = EQ->nphi;
    ntht = EQ->ntheta;
    dphi = EQ->dphi;
    N = nphi * ntht;

    abs_square_a = rarrDef(N);
    abs_square_b = rarrDef(N);

    E = functionals(EQ, Sa, Sb, &kin, &old_mu_a, &old_mu_b);
    stationaryNewton(EQ, Sa, Sb, err_tol, iter_tol, old_mu_a, old_mu_b);
    carrAbs2(N, Sa, abs_square_a);
    carrAbs2(N, Sb, abs_square_b);
    old_norm_a = sqrt(
            Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_a, dphi)
    );
    old_norm_b = sqrt(
            Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_b, dphi)
    );

    renormalize_spheric(EQ, Sa);
    renormalize_spheric(EQ, Sb);

    while (fabs(old_norm_a - 1) > err_tol || fabs(old_norm_b - 1) > err_tol)
    {
        printf(
                "\n\nstep result : %.6lf %.6lf %.6lf %.6lf %.6lf\n\n",
                E, old_mu_a, old_mu_b, old_norm_a, old_norm_b
        );

        E = functionals(EQ, Sa, Sb, &kin, &mu_a, &mu_b);
        stationaryNewton(EQ, Sa, Sb, 0.5 * err_tol, iter_tol, mu_a, mu_b);
        carrAbs2(N, Sa, abs_square_a);
        carrAbs2(N, Sb, abs_square_b);
        norm_a = sqrt(
                Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_a, dphi)
        );
        norm_b = sqrt(
                Rsimps2D_sphere(nphi, ntht, EQ->theta, abs_square_b, dphi)
        );

        der_a = (norm_a - old_norm_a) / (mu_a - old_mu_a);
        der_b = (norm_b - old_norm_b) / (mu_b - old_mu_b);
        aux = mu_a;
        mu_a = old_mu_a - (1 - norm_a) / der_a;
        old_mu_a = aux;
        aux = mu_b;
        mu_b = old_mu_b - (1 - norm_b) / der_b;
        old_mu_b = aux;
        old_norm_a = norm_a;
        old_norm_b = norm_b;
    }

    printf(
            "\n\nstep result : %.6lf %.6lf %.6lf %.6lf %.6lf\n\n",
            E, old_mu_a, old_mu_b, old_norm_a, old_norm_b
    );

    free(abs_square_a);
    free(abs_square_b);

}
