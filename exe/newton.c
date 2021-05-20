#include <string.h>
#include <stdio.h>
#include <math.h>
#include "newtoncg.h"


void TimePrint(double t)
{

    // format and print time in days / hours / minutes

    int
        tt = (int) t,
        hours = 0,
        mins  = 0;

    if ( tt / 3600  > 0 )
    { hours = tt / 3600;  tt = tt % 3600;  }

    if ( tt / 60    > 0 )
    { mins  = tt / 60;    tt = tt % 60;    }

    printf(" %d hour(s) %d minute(s)",hours,mins);
}


void ReachNewLine(FILE * f)
{
    /* Read until get new line or 'end of file' in a opened file. */
    char
        sentinel;
    while (1)
    {
        fscanf(f, "%c", &sentinel);
        if (sentinel == '\n' || sentinel == EOF) return;
    }
}


void save_equation_setup(
        char prefix [],
        EqDataPkg EQ,
        double mu_a,
        double mu_b,
        int azi_a,
        int azi_b,
        int line)
{
    char
        fname[100];
    FILE
        * txt_file_ptr;

    strcpy(fname, "output/");
    strcat(fname, prefix);
    strcat(fname, "_equation_newton.dat");

    if (line > 1) txt_file_ptr = fopen(fname, "a");
    else          txt_file_ptr = fopen(fname, "w");

    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    // Equation parameters
    fprintf(
            txt_file_ptr,
            "%d %.2lf %.2lf %.10lf %.10lf %d %d\n",
            EQ->ntheta, EQ->nabla_coef, EQ->omega, mu_a, mu_b, azi_a, azi_b
    );

    fclose(txt_file_ptr);
}


void save_obs_2species(
        char prefix [], EqDataPkg EQ, Rarray Sa, Rarray Sb, int azi_a,
        int azi_b, int line)
{
    int
        ntheta;
    double
        ga,
        gb,
        gab,
        norm_a,
        norm_b,
        frac_a,
        frac_b,
        energy,
        kin_energy,
        den_overlap;
    Rarray
        sin_th,
        abs_square_a,
        abs_square_b;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    ntheta = EQ->ntheta;

    sin_th = rarrDef(ntheta);
    for (int i = 0; i < ntheta; i++)
    {
        sin_th[i] = sin(EQ->theta[i]);
    }

    abs_square_a = rarrDef(ntheta);
    abs_square_b = rarrDef(ntheta);

    rarrAbs2(ntheta, Sa, abs_square_a);
    rarrAbs2(ntheta, Sb, abs_square_b);

    norm_a = Rsimps1D_jac(ntheta, abs_square_a, EQ->dtheta, sin_th);
    norm_b = Rsimps1D_jac(ntheta, abs_square_b, EQ->dtheta, sin_th);
    frac_a = norm_a / (norm_a + norm_b);
    frac_b = norm_b / (norm_a + norm_b);
    ga = EQ->ga * norm_a * (2 * PI);
    gb = EQ->gb * norm_b * (2 * PI);
    gab = EQ->gab * sqrt(norm_b * norm_b) * (2 * PI);

    energy = functionals_newton(EQ, Sa, Sb, &kin_energy, azi_a, azi_b);
    den_overlap = theta_density_overlap(EQ, abs_square_a, abs_square_b);

    strcpy(fname, "output/");
    strcat(fname, prefix);
    strcat(fname, "_obs_newton.dat");

    if (line > 1) txt_file_ptr = fopen(fname, "a");
    else          txt_file_ptr = fopen(fname, "w");

    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    fprintf(txt_file_ptr,
            "%.10lf %.10lf %.10lf %.10lf %.10lf ",
            energy, kin_energy, den_overlap, norm_a, norm_b
    );
    fprintf(txt_file_ptr,
            "%.10lf %.10lf %.10lf %.10lf %.10lf\n",
            frac_a, frac_b, ga, gb, gab
    );

    free(sin_th);
    free(abs_square_b);
    free(abs_square_a);
    fclose(txt_file_ptr);
}


EqDataPkg setup_equation(
        char prefix [], int num_species, int line, int * azi_a, int * azi_b,
        double * mu_a, double * mu_b)
{

/** Read line by line of _domain file and _eq to setup the equation **/

    int
        scanf_returned,
        theta_grid_size;
    double
        nabla_coef,
        rotation_freq,
        frac_a,
        frac_b,
        ga,
        gb,
        gab;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    // Try openning file with equation parameters
    strcpy(fname, "input/");
    strcat(fname, prefix);
    strcat(fname, "_eq.dat");
    printf("\nLooking for %s", fname);
    txt_file_ptr = fopen(fname, "r");
    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ................ Found !");
    }

    // advance to the requested line
    for (int i = 1; i < line; i++) ReachNewLine(txt_file_ptr);

    scanf_returned = fscanf(
            txt_file_ptr, "%lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf",
            &nabla_coef, &rotation_freq, &frac_a, &frac_b, &ga, &gb, &gab,
            azi_a, azi_b, mu_a, mu_b
    );
    ga = ga / (2 * PI);
    gb = gb / (2 * PI);
    gab = gab / (2 * PI);

    fclose(txt_file_ptr);

    if (scanf_returned != 11)
    {
        printf("\n\nWrong number of parameters or bad format in %s", fname);
        exit(EXIT_FAILURE);
    }

    // Try openning file to configure grid domain
    strcpy(fname, "input/");
    strcat(fname, prefix);
    strcat(fname, "_domain.dat");
    printf("\nLooking for %s", fname);
    txt_file_ptr = fopen(fname, "r");
    if (txt_file_ptr == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ............ Found !");
    }

    // advance to the requested line
    for (int i = 1; i < line; i++) ReachNewLine(txt_file_ptr);

    scanf_returned = fscanf(txt_file_ptr, "%d", &theta_grid_size);

    fclose(txt_file_ptr);

    if (scanf_returned != 1)
    {
        printf("\n\nWrong number of parameters or bad format in %s", fname);
        exit(EXIT_FAILURE);
    }

    return equation_structure(
            1, theta_grid_size, 0, 0,
            nabla_coef, rotation_freq, frac_a, frac_b, ga, gb, gab
    );

}


void setup_initial_condition(
        EqDataPkg EQ, char prefix [], Rarray state_a, Rarray state_b)
{
    int
        i,
        ntheta,
        fscanf_val;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    ntheta = EQ->ntheta;

    // open file with initial condition set on the grid - species A
    strcpy(fname, "input/");
    strcat(fname, prefix);
    strcat(fname, "_speciesA_init.dat");
    printf("\nLooking for %s", fname);
    txt_file_ptr = fopen(fname, "r");
    if (txt_file_ptr == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ..... Found !");
    }

    for (i = 0; i < ntheta; i++)
    {
        fscanf_val = fscanf(txt_file_ptr, "%lf", &state_a[i]);
        if (fscanf_val != 1)
        {
            printf("\n\nProblem reading initial condition data: ");
            printf("%d data points read successfully\n\n", i);
            exit(EXIT_FAILURE);
        }
    }
    fclose(txt_file_ptr);

    // open file with initial condition set on the grid - species B
    strcpy(fname, "input/");
    strcat(fname, prefix);
    strcat(fname, "_speciesB_init.dat");
    printf("\nLooking for %s", fname);
    txt_file_ptr = fopen(fname, "r");
    if (txt_file_ptr == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf(" ..... Found !");
    }

    for (i = 0; i < ntheta; i++)
    {
        fscanf_val = fscanf(txt_file_ptr, "%lf", &state_b[i]);
        if (fscanf_val != 1)
        {
            printf("\n\nProblem reading initial condition data: ");
            printf("%d data points read successfully\n\n", i);
            exit(EXIT_FAILURE);
        }
    }
    fclose(txt_file_ptr);

}





int main(int argc, char * argv[])
{

    /*  DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER ARCHITECTURE
     *  --------------------------------------------------------------- */

    // omp_set_num_threads(omp_get_max_threads() / 2);
    omp_set_num_threads(1);
    mkl_set_num_threads(1);

    int
        i,
        njobs,
        num_species,
        azi_a,
        azi_b;
    double
        mu_a,
        mu_b,
        start,      // start trigger to measure time
        time_used;  // Time used in calling evolution routine
    char
        c,
        line_str[20],
        fname[100],     // name of files to open
        infname[100],   // input file name prefix
        outfname[100];  // output file name prefix
    FILE
        * txt_file_ptr;  // grid and time information
    Rarray
        Sa,
        Sb;
    EqDataPkg
        EQ;

    txt_file_ptr = fopen("job.conf", "r");
    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n", "job.conf");
        exit(EXIT_FAILURE);
    }

    i = 1;
    while ((c  = getc(txt_file_ptr)) != EOF)
    {
        // jump comment lines marked started with #
        if (c == '#') { ReachNewLine(txt_file_ptr); continue; }
        else          { fseek(txt_file_ptr, -1, SEEK_CUR);    }
        switch (i)
        {
            case 1:
                fscanf(txt_file_ptr, "%d", &num_species);
                i = i + 1;
                break;
            case 2:
                fscanf(txt_file_ptr, "%s", infname);
                i = i + 1;
                break;
            case 3:
                fscanf(txt_file_ptr, "%s", outfname);
                break;
        }
        ReachNewLine(txt_file_ptr);
    }

    fclose(txt_file_ptr);

    if (i != 3)
    {
        printf("\nWrong number of parameter read from job.conf file\n\n");
        exit(EXIT_FAILURE);
    }

    if (num_species != 2)
    {
        printf("\nInvalid number of species %d in job.conf\n\n", num_species);
        exit(EXIT_FAILURE);
    }

    printf("\n\nNewton with conjugate gradient method");

    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_eq.dat");
    njobs = NumberOfLines(fname);

    printf("\nNumber of jobs requested : %d\n", njobs);

    for (i = 1; i <= njobs; i++)
    {
        // line number to read
        sprintf(line_str, "%d", i);

        printf("\n\n");
        sepline();
        printf("Initiating job %d ...", i);

        // PACK DOMAIN AND PARAMETERS INFORMATION IN A STRUCTURE
        EQ = setup_equation(infname, num_species, i, &azi_a, &azi_b, &mu_a, &mu_b);

        save_equation_setup(outfname, EQ, mu_a, mu_b, azi_a, azi_b, i);

        if (i == 1)
        {
            // configure initial condition
            Sa = rarrDef(EQ->ntheta);
            Sb = rarrDef(EQ->ntheta);
            setup_initial_condition(EQ, infname, Sa, Sb);
        }

        printf("\n\n");
        printf("\t\t*********************************************\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*            START NEWTON METHOD            *\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*********************************************\n");

        start = omp_get_wtime();

        stationaryNewton(EQ, Sa, Sb, 5E-4, 99, mu_a, mu_b, azi_a, azi_b);

        time_used = (double) (omp_get_wtime() - start);
        printf("\nTime elapsed");
        printf(" : %.0lf sec = ",time_used);
        TimePrint(time_used);

        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_speciesA_job");
        strcat(fname, line_str);
        strcat(fname, "_newton.dat");
        rarr_txt(fname, EQ->ntheta, Sa);
        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_speciesB_job");
        strcat(fname, line_str);
        strcat(fname, "_newton.dat");
        rarr_txt(fname, EQ->ntheta, Sb);
        save_obs_2species(outfname, EQ, Sa, Sb, azi_a, azi_b, i);

        // free(Sa);
        // free(Sb);
        ReleaseEqDataPkg(EQ);
    }

    printf("\n\n");
    return 0;
}
