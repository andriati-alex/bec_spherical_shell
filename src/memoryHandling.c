#include "memoryHandling.h"



Rarray rarrDef(int n)
{
    double * ptr;

    ptr = (double * ) malloc( n * sizeof(double) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for double\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Carray carrDef(int n)
{
    double complex * ptr;

    ptr = (double complex * ) malloc( n * sizeof(double complex) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



CMKLarray cmklDef(int n)
{
    MKL_Complex16 * ptr;

    ptr = (MKL_Complex16 *) malloc( n * sizeof(MKL_Complex16) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex(mkl)\n\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}



Rmatrix rmatDef(int m, int n)
{

/** Real matrix of m rows and n columns **/

    int i;

    double ** ptr;

    ptr = (double ** ) malloc( m * sizeof(double *) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for real matrix\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = rarrDef(n);

    return ptr;
}



Cmatrix cmatDef(int m, int n)
{

/** Complex matrix of m rows and n columns **/

    int i;

    double complex ** ptr;

    ptr = (double complex ** ) malloc( m * sizeof(double complex *) );

    if (ptr == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for complex matrix\n\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; i++) ptr[i] = carrDef(n);

    return ptr;
}



CCScmat ccscmatDef(int n, int max_nonzeros)
{

/** Return empty CCS representation of matrix of  n  rows
  * and the maximum number of nonzeros elements in a same
  * row given in 'max_nonzeros'                       **/

    CCScmat
        M;

    M = (struct _CCScmat *) malloc(sizeof(struct _CCScmat));

    if (M == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for CCS matrix\n\n");
        exit(EXIT_FAILURE);
    }

    M->m = max_nonzeros;
    M->vec = carrDef(max_nonzeros * n);
    M->col = (int *) malloc( max_nonzeros * n * sizeof(int) );

    if (M->col == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integers\n\n");
        exit(EXIT_FAILURE);
    }

    return M;
}



CCSrmat ccsrmatDef(int n, int max_nonzeros)
{

    CCSrmat
        M;

    M = (struct _CCSrmat *) malloc(sizeof(struct _CCSrmat));

    if (M == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for CCS matrix\n\n");
        exit(EXIT_FAILURE);
    }

    M->m = max_nonzeros;
    M->vec = rarrDef(max_nonzeros * n);
    M->col = (int *) malloc( max_nonzeros * n * sizeof(int) );

    if (M->col == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for integers\n\n");
        exit(EXIT_FAILURE);
    }

    return M;
}



EqDataPkg equation_structure(
        int phi_grid_size, int theta_grid_size, double time_step_size,
        int time_steps, double nabla_coef, double rotation_freq,
        double frac_a, double frac_b, double ga, double gb, double gab)
{

/** Return pointer to a basic data structure with all needed information
    to solve the Gross-Pitaevskii equation                           **/

    EqDataPkg Equation = (EqDataPkg) malloc(sizeof(struct _EquationDataPkg));

    if (Equation == NULL)
    {
        printf("\n\n\nMEMORY ERROR : malloc fail for EqData structure\n\n");
        exit(EXIT_FAILURE);
    }

    Equation->phi = rarrDef(phi_grid_size);
    Equation->theta = rarrDef(theta_grid_size);

    Equation->nphi = phi_grid_size;
    Equation->ntheta = theta_grid_size;
    Equation->dphi = 2 * PI / (phi_grid_size - 1);
    Equation->dtheta = PI / (theta_grid_size - 1);
    Equation->nt = time_steps;
    Equation->dt = time_step_size;
    Equation->frac_a = frac_a;
    Equation->frac_b = frac_b;
    Equation->ga = ga;
    Equation->gb = gb;
    Equation->gab = gab;
    Equation->nabla_coef = nabla_coef;
    Equation->omega = rotation_freq;

    rarrFillInc(Equation->nphi, 0, Equation->dphi, Equation->phi);
    rarrFillInc(Equation->ntheta, 0, Equation->dtheta, Equation->theta);

    return Equation;
}


TwoSpeciesState alloc_two_species_struct(int nphi, int ntheta)
{
    TwoSpeciesState S = (TwoSpeciesState) malloc(
            sizeof(struct _TwoSpeciesState)
    );

    if (S == NULL)
    {
        printf("\n\n\nMEMORY ERROR : malloc fail for TwoSpeciesState\n\n");
        exit(EXIT_FAILURE);
    }

    S->nphi = nphi;
    S->ntht = ntheta;
    S->speca = carrDef(nphi * ntheta);
    S->specb = carrDef(nphi * ntheta);
    S->speca_re = rarrDef(nphi * ntheta);
    S->speca_im = rarrDef(nphi * ntheta);
    S->specb_re = rarrDef(nphi * ntheta);
    S->specb_im = rarrDef(nphi * ntheta);

    return S;
}


void pkg_states(Carray Sa, Carray Sb, TwoSpeciesState S)
{
    carrCopy(S->nphi * S->ntht, Sa, S->speca);
    carrCopy(S->nphi * S->ntht, Sb, S->specb);
    carrRealPart(S->nphi * S->ntht, S->speca, S->speca_re);
    carrImagPart(S->nphi * S->ntht, S->speca, S->speca_im);
    carrRealPart(S->nphi * S->ntht, S->specb, S->specb_re);
    carrImagPart(S->nphi * S->ntht, S->specb, S->specb_im);
}


void unpkg_states(Carray Sa, Carray Sb, TwoSpeciesState S)
{
    carrCopy(S->nphi * S->ntht, S->speca, Sa);
    carrCopy(S->nphi * S->ntht, S->specb, Sb);
}


void set_states_from_parts(TwoSpeciesState S)
{
    int
        grid_pt;
    for (int j = 0; j < S->ntht; j++)
    {
        for (int i = 0; i < S->nphi; i++)
        {
            grid_pt = j * S->nphi + i;
            S->speca[grid_pt] =  (
                    S->speca_re[grid_pt] + I * S->speca_im[grid_pt]
            );
            S->specb[grid_pt] =  (
                    S->specb_re[grid_pt] + I * S->specb_im[grid_pt]
            );
        }
    }
}



void set_states_real_imag(TwoSpeciesState S)
{
    int
        grid_pt;
    for (int j = 0; j < S->ntht; j++)
    {
        for (int i = 0; i < S->nphi; i++)
        {
            grid_pt = j * S->nphi + i;
            S->speca_re[grid_pt] = creal(S->speca[grid_pt]);
            S->speca_im[grid_pt] = cimag(S->speca[grid_pt]);
            S->specb_re[grid_pt] = creal(S->specb[grid_pt]);
            S->specb_im[grid_pt] = cimag(S->specb[grid_pt]);
        }
    }
}





/* ========================================================================
   ========================================================================
 
                               MEMORY RELEASE

   ========================================================================
   ======================================================================== */





void rmatFree(int m, Rmatrix M)
{

/** Release a real matrix of m rows **/

    int i;

    for (i = 0; i < m; i++) free(M[i]);

    free(M);
}



void cmatFree(int m, Cmatrix M)
{

/** Release a complex matrix of m rows **/

    int i;

    for (i = 0; i < m; i++) free(M[i]);

    free(M);
}



void ccscmatFree(CCScmat M)
{

/** Release Compressed-Column Storaged matrix **/

    free(M->col);
    free(M->vec);
    free(M);
}



void ccsrmatFree(CCSrmat M)
{

/** Release Compressed-Column Storaged matrix **/

    free(M->col);
    free(M->vec);
    free(M);
}



void ReleaseEqDataPkg(EqDataPkg EQ)
{
    free(EQ->phi);
    free(EQ->theta);
    free(EQ);
}


void release_two_species_state(TwoSpeciesState S)
{
    free(S->speca);
    free(S->specb);
    free(S->speca_re);
    free(S->speca_im);
    free(S->specb_re);
    free(S->specb_im);
    free(S);
}
