#include "newtoncg.h"


double inner_cg(TwoSpeciesState s1, TwoSpeciesState s2, TwoSpeciesState aux,
        Rarray theta)
{

/** Inner product of two state of two component system. Use real and imag
    time separately and integrate over the sphere                    **/

    return (
        rarrDot(s1->ntht, s1->speca_re, s2->speca_re) +
        rarrDot(s1->ntht, s1->specb_re, s2->specb_re)
    );

}

void update_state_cg(
        TwoSpeciesState s1,
        double scalar,
        TwoSpeciesState s2,
        TwoSpeciesState out)
{
    rarrUpdate(s1->ntht, s1->speca_re, scalar, s2->speca_re, out->speca_re);
    rarrUpdate(s1->ntht, s1->specb_re, scalar, s2->specb_re, out->specb_re);
    rarrFill(out->ntht, 0.0, out->speca_im);
    rarrFill(out->ntht, 0.0, out->specb_im);
    set_states_from_parts(out);
}

double self_inner_cg(TwoSpeciesState s, TwoSpeciesState aux,
        Rarray theta)
{
    return (
        rarrDot(s->ntht, s->speca_re, s->speca_re) +
        rarrDot(s->ntht, s->specb_re, s->specb_re)
    );
}

void twice_derivative_2order(int n, Rarray f, double dx, Rarray d2f)
{
    d2f[0] = 2 * (f[1] - f[0]) / dx / dx;
    d2f[n - 1] = 2 * (f[n - 2] - f[n - 1]) / dx /dx;
    for (int i = 1; i < n - 1; i++)
    {
        d2f[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / dx / dx;
    }
}

void derivative_2order(int n, Rarray f, double dx, Rarray d2f)
{
    d2f[0] = 0;
    d2f[n - 1] = 0;
    for (int i = 1; i < n - 1; i++)
    {
        d2f[i] = (f[i + 1] - f[i - 1]) / 2 / dx;
    }
}

void laplace_2order(
        int n,
        Rarray f,
        int vort,
        double dx,
        Rarray tht,
        Rarray nabla,
        Rarray work1,
        Rarray work2)
{
    int
        vort_sq;
    double
        cot_tht,
        sin_tht_sq;

    vort_sq = vort * vort;
    derivative_2order(n, f, dx, work1);
    twice_derivative_2order(n, f, dx, work2);

    nabla[0] = work2[0];
    nabla[n - 1] = work2[n - 1];
    for (int i = 1; i < n - 1; i++)
    {
        cot_tht = cos(tht[i]) / sin(tht[i]);
        sin_tht_sq = sin(tht[i]) * sin(tht[i]);
        nabla[i] = work2[i] + cot_tht * work1[i] - vort_sq * f[i] / sin_tht_sq;
    }
}

void grid_residue(
        EqDataPkg EQ,
        Rarray state_a,
        Rarray state_b,
        int vort_a,
        int vort_b,
        Rarray grid_res_a,
        Rarray grid_res_b,
        double mu_a,
        double mu_b)
{
    int
        j,
        ntht;
    double
        nabla_coef,
        dtht,
        ga,
        gb,
        gab;
    Rarray
        tht_vals,
        abs_square_a,
        abs_square_b;
    Rarray
        work1,
        work2,
        nabla;

    // unpack equation data
    ntht = EQ->ntheta;
    dtht = EQ->dtheta;
    tht_vals = EQ->theta;
    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;

    nabla = rarrDef(ntht);
    work1 = rarrDef(ntht);
    work2 = rarrDef(ntht);
    abs_square_a = rarrDef(ntht);
    abs_square_b = rarrDef(ntht);

    rarrAbs2(ntht, state_a, abs_square_a);
    rarrAbs2(ntht, state_b, abs_square_b);

    // Compute grid residue for species A
    laplace_2order(
            ntht, state_a, vort_a, dtht, tht_vals, nabla, work1, work2
    );
    for (j = 0; j < ntht; j++)
    {
        grid_res_a[j] = (
                nabla_coef * nabla[j] +
                state_a[j] * (ga * abs_square_a[j] + gab * abs_square_b[j])
                - mu_a * state_a[j]
        );
    }

    // Compute grid residue for species B
    laplace_2order(
            ntht, state_b, vort_b, dtht, tht_vals, nabla, work1, work2
    );
    for (j = 0; j < ntht; j++)
    {
        grid_res_b[j] = (
                nabla_coef * nabla[j] +
                state_b[j] * (gb * abs_square_b[j] + gab * abs_square_a[j])
                - mu_b * state_b[j]
        );
    }
    if (vort_a != 0)
    {
        grid_res_a[0] = 0;
        grid_res_a[ntht - 1] = 0;
    }
    if (vort_b != 0)
    {
        grid_res_b[0] = 0;
        grid_res_b[ntht - 1] = 0;
    }

    free(nabla);
    free(work1);
    free(work2);
    free(abs_square_a);
    free(abs_square_b);
}

void linearized_op(
        EqDataPkg EQ,
        TwoSpeciesState newton,
        TwoSpeciesState conj_grad,
        TwoSpeciesState lin_op,
        double mu_a,
        double mu_b,
        int vort_a,
        int vort_b
        )
{

    int
        j,
        ntht;
    double
        nabla_coef,
        ga,
        gb,
        gab,
        dtht,
        fa,
        fb,
        fa_sq,
        fb_sq,
        cg_fa,
        cg_fb;
    Rarray
        tht_vals,
        work1,
        work2,
        nabla_a,
        nabla_b,
        conj_grad_a,
        conj_grad_b;

    conj_grad_a = conj_grad->speca_re;
    conj_grad_b = conj_grad->specb_re;

    // unpack Equation data
    nabla_coef = EQ->nabla_coef;
    ga = EQ->ga;
    gb = EQ->gb;
    gab = EQ->gab;
    dtht = EQ->dtheta;
    ntht = EQ->ntheta;
    tht_vals = EQ->theta;

    work1 = rarrDef(ntht);
    work2 = rarrDef(ntht);
    nabla_a = rarrDef(ntht);
    nabla_b = rarrDef(ntht);
    laplace_2order(
            ntht, conj_grad_a, vort_a, dtht, tht_vals, nabla_a, work1, work2
    );
    laplace_2order(
            ntht, conj_grad_b, vort_b, dtht, tht_vals, nabla_b, work1, work2
    );

    for (j = 0; j < ntht; j++)
    {
        fa = newton->speca_re[j];
        fb = newton->specb_re[j];
        cg_fa = conj_grad_a[j];
        cg_fb = conj_grad_b[j];
        fa_sq = fa * fa;
        fb_sq = fb * fb;
        lin_op->speca_re[j] = (
                nabla_coef * nabla_a[j] +
                cg_fa * (- mu_a + ga * 3 * fa_sq + gab * fb_sq) +
                cg_fb * (2 * gab * fb * fa)
        );
        lin_op->specb_re[j] = (
                nabla_coef * nabla_b[j] +
                cg_fb * (- mu_b + gb * 3 * fb_sq + gab * fa_sq) +
                cg_fa * (2 * gab * fb * fa)
        );
    }
    if (vort_a != 0)
    {
        lin_op->speca_re[0] = 0;
        lin_op->speca_re[ntht - 1] = 0;
    }
    if (vort_b != 0)
    {
        lin_op->specb_re[0] = 0;
        lin_op->specb_re[ntht - 1] = 0;
    }

    free(work1);
    free(work2);
    free(nabla_a);
    free(nabla_b);
}

int conjgrad_stab(
        EqDataPkg EQ,
        double mu_a,
        double mu_b,
        int vort_a,
        int vort_b,
        TwoSpeciesState newton,
        TwoSpeciesState newton_res,
        TwoSpeciesState conj_grad,
        double tol
        )
{
    int
        ntht,
        it_counter;
    Rarray
        tht;
    double
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

    ntht = EQ->ntheta;
    tht = EQ->theta;
    res = alloc_two_species_struct(1, ntht);
    prev_res = alloc_two_species_struct(1, ntht);
    s_aux = alloc_two_species_struct(1, ntht);
    dir = alloc_two_species_struct(1, ntht);
    lin_op = alloc_two_species_struct(1, ntht);
    lin_op_s = alloc_two_species_struct(1, ntht);
    inner = alloc_two_species_struct(1, ntht);
    res_star = alloc_two_species_struct(1, ntht);

    rarrFill(ntht, 0.0, res->speca_im);
    rarrFill(ntht, 0.0, res->specb_im);
    rarrFill(ntht, 0.0, lin_op->speca_im);
    rarrFill(ntht, 0.0, lin_op->specb_im);
    rarrFill(ntht, 0.0, res_star->speca_im);
    rarrFill(ntht, 0.0, res_star->specb_im);
    rarrFill(ntht, 0.0, dir->speca_im);
    rarrFill(ntht, 0.0, dir->specb_im);
    rarrFill(ntht, 0.0, s_aux->speca_im);
    rarrFill(ntht, 0.0, s_aux->specb_im);
    rarrFill(ntht, 0.0, inner->speca_im);
    rarrFill(ntht, 0.0, inner->specb_im);
    rarrFill(ntht, 0.0, lin_op_s->speca_im);
    rarrFill(ntht, 0.0, lin_op_s->specb_im);

    linearized_op(EQ, newton, conj_grad, lin_op, mu_a, mu_b, vort_a, vort_b);
    update_state_cg(newton_res, -1.0, lin_op, res);

    set_states_from_parts(lin_op);
    set_states_from_parts(res);

    pkg_real_states(res->speca_re, res->specb_re, dir);
    pkg_real_states(res->speca_re, res->specb_re, res_star);

    err = sqrt(self_inner_cg(res, inner, tht));

    it_counter = 0;
    // printf("\n\titer %5d error : %.7lf", it_counter, err);

    while (err > tol)
    {
        linearized_op(EQ, newton, dir, lin_op, mu_a, mu_b, vort_a, vort_b);
        if (inner_cg(res_star, lin_op, inner, tht) == 0)
        {
            printf("\n\nZero Division in alpha scalar!\n\n");
            exit(EXIT_FAILURE);
        }
        alpha_scalar = (
                inner_cg(res_star, res, inner, tht) /
                inner_cg(res_star, lin_op, inner, tht)
        );
        update_state_cg(res, -alpha_scalar, lin_op, s_aux);
        linearized_op(EQ, newton, s_aux, lin_op_s, mu_a, mu_b, vort_a, vort_b);
        if (self_inner_cg(lin_op_s, inner, tht) == 0)
        {
            printf("\n\nZero Division in omega scalar!\n\n");
            exit(EXIT_FAILURE);
        }
        omega_scalar = (
                inner_cg(s_aux, lin_op_s, inner, tht) /
                self_inner_cg(lin_op_s, inner, tht)
        );
        // update solution - use `prev_res` as step workspace variable
        update_state_cg(conj_grad, alpha_scalar, dir, prev_res);
        update_state_cg(prev_res, omega_scalar, s_aux, conj_grad);
        if (omega_scalar == 0)
        {
            err = sqrt(self_inner_cg(res, inner, tht));
            printf("\n\nOmega == 0 found. Exiting with error %.8lf\n", err);
            exit(EXIT_FAILURE);
        }
        // residual update
        pkg_states(res->speca, res->specb, prev_res);
        update_state_cg(s_aux, -omega_scalar, lin_op_s, res);
        if (inner_cg(res_star, prev_res, inner, tht) == 0)
        {
            printf("\n\nZero Division in beta scalar!\n\n");
            exit(EXIT_FAILURE);
        }
        beta_scalar = (
                inner_cg(res_star, res, inner, tht) /
                inner_cg(res_star, prev_res, inner, tht)
        ) * alpha_scalar / omega_scalar;
        // update direction - use `prev_res` as step workspace variable
        update_state_cg(res, beta_scalar, dir, prev_res);
        update_state_cg(
                prev_res, -beta_scalar * omega_scalar, lin_op, dir
        );
        // update solution error
        err = sqrt(self_inner_cg(res, inner, tht));

        it_counter++;
        /*
        if (it_counter % 100 == 0)
        {
            printf("\n\titer %5d error : %.7lf", it_counter, err);
            linearized_op(
                EQ, newton, conj_grad, lin_op, mu_a, mu_b, vort_a, vort_b
            );
            update_state_cg(newton_res, -1.0, lin_op, prev_res);
            err = sqrt(self_inner_cg(prev_res, inner, tht));
            printf("\n\titer %5d error : %.7lf", it_counter, err);
        }
        */
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


void stationaryNewton(
        EqDataPkg EQ,
        Rarray Sa,
        Rarray Sb,
        double err_tol,
        int iter_tol,
        double mu_a,
        double mu_b,
        int vort_a,
        int vort_b)
{

    int
        i,
        N,
        Niter,
        CGiter;
    double
        dtht,
        norm1,
        norm2,
        cg_error,
        error_newton;
    Rarray
        sin_th,
        grid_res_a,
        grid_res_b,
        abs_square_a,
        abs_square_b;
    TwoSpeciesState
        workspace,
        conj_grad,
        newton_step,
        newton_res;

    N = EQ->ntheta;
    dtht = EQ->dtheta;

    sin_th = rarrDef(N);
    for (i = 0; i < N; i++)
    {
        sin_th[i] = sin(EQ->theta[i]);
    }
    printf("\nNewton method for mu_a = %.5lf and mu_b = %.5lf", mu_a, mu_b);

    workspace = alloc_two_species_struct(1, N);
    newton_res = alloc_two_species_struct(1, N);
    newton_step = alloc_two_species_struct(1, N);
    conj_grad = alloc_two_species_struct(1, N);
    grid_res_a = rarrDef(N);
    grid_res_b = rarrDef(N);
    abs_square_a = rarrDef(N);
    abs_square_b = rarrDef(N);

    grid_residue(EQ, Sa, Sb, vort_a, vort_b, grid_res_a, grid_res_b, mu_a, mu_b);

    rarrAbs2(N, Sa, abs_square_a);
    rarrAbs2(N, Sb, abs_square_b);
    norm1 = Rsimps1D_jac(N, abs_square_a, dtht, sin_th);
    norm2 = Rsimps1D_jac(N, abs_square_b, dtht, sin_th);

    // right hand side of linear operator passed to iteratice CG method
    for (i = 0; i < N; i++)
    {
        grid_res_a[i] = - grid_res_a[i];
        grid_res_b[i] = - grid_res_b[i];
    }
    pkg_real_states(grid_res_a, grid_res_b, newton_res);
    error_newton = sqrt(self_inner_cg(newton_res, workspace, EQ->theta));

    printf("\nIt    Error      CGit   norm1      norm2");
    sepline();

    // NEWTON LOOP
    Niter = 0;
    CGiter = 0;
    while (error_newton > err_tol)
    {
        printf("%2d   %9.6lf  %5d  ", Niter, error_newton, CGiter);
        printf("%9.6lf  %9.6lf\n", norm1, norm2);

        // Initial guess for conjugate-gradient method
        carrFill(N, 0.0, conj_grad->speca);
        carrFill(N, 0.0, conj_grad->specb);
        set_states_real_imag(conj_grad);

        pkg_real_states(Sa, Sb, newton_step);

        cg_error = 0.1 * error_newton;
        CGiter = conjgrad_stab(
                EQ,
                mu_a,
                mu_b,
                vort_a,
                vort_b,
                newton_step,
                newton_res,
                conj_grad,
                cg_error
        );

        // update solution from newton iteration
        rarrAdd(N, Sa, conj_grad->speca_re, Sa);
        rarrAdd(N, Sb, conj_grad->specb_re, Sb);

        grid_residue(
                EQ, Sa, Sb, vort_a, vort_b, grid_res_a, grid_res_b, mu_a, mu_b
        );

        for (i = 0; i < N; i++)
        {
            grid_res_a[i] = - grid_res_a[i];
            grid_res_b[i] = - grid_res_b[i];
        }
        pkg_real_states(grid_res_a, grid_res_b, newton_res);
        error_newton = sqrt(
                self_inner_cg(newton_res, workspace, EQ->theta)
        );

        Niter = Niter + 1;

        if (Niter > iter_tol) break;

        rarrAbs2(N, Sa, abs_square_a);
        rarrAbs2(N, Sb, abs_square_b);
        norm1 = Rsimps1D_jac(N, abs_square_a, dtht, sin_th);
        norm2 = Rsimps1D_jac(N, abs_square_b, dtht, sin_th);
    }

    printf("%2d   %9.6lf  %5d  ", Niter, error_newton, CGiter);
    printf("%9.6lf  %9.6lf", norm1, norm2);
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
