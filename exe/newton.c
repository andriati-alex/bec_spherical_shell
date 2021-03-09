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
        char prefix [], EqDataPkg EQ, int num_species, int line)
{
    char
        fname[100];
    FILE
        * txt_file_ptr;

    strcpy(fname, "output/");
    strcat(fname, prefix);
    if (num_species == 1) strcat(fname, "_1species_equation_imagtime.dat");
    else                  strcat(fname, "_2species_equation_imagtime.dat");

    if (line > 1) txt_file_ptr = fopen(fname, "a");
    else          txt_file_ptr = fopen(fname, "w");

    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    // Grid domain
    fprintf(txt_file_ptr,
            "%d %d %lf %d ",
            EQ->nphi, EQ->ntheta, EQ->dt, EQ->nt
    );

    // Equation parameters
    if (num_species == 2)
    {
        fprintf(
                txt_file_ptr,
                "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",
                EQ->nabla_coef, EQ->omega, EQ->frac_a, EQ->frac_b,
                EQ->ga, EQ->gb, EQ->gab
        );
    }
    else
    {
        fprintf(txt_file_ptr,
                "%.15lf %.15lf %.15lf\n",
                EQ->nabla_coef, EQ->omega, EQ->ga
        );
    }
    fclose(txt_file_ptr);
}


void save_obs_2species(
        char prefix [], EqDataPkg EQ, Carray Sa, Carray Sb, int line)
{
    int
        nphi,
        ntheta;
    double
        dphi,
        lza,
        lzb,
        mu_a,
        mu_b,
        energy,
        kin_energy,
        den_overlap;
    Rarray
        theta,
        abs_square_a,
        abs_square_b;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    nphi = EQ->nphi;
    ntheta = EQ->ntheta;
    dphi = EQ->dphi;
    theta = EQ->theta;

    abs_square_a = rarrDef(nphi * ntheta);
    abs_square_b = rarrDef(nphi * ntheta);
    carrAbs2(nphi * ntheta, Sa, abs_square_a);
    carrAbs2(nphi * ntheta, Sb, abs_square_b);

    energy = functionals(EQ, Sa, Sb, &kin_energy, &mu_a, &mu_b);
    den_overlap = density_overlap(EQ, abs_square_a, abs_square_b);
    lza = angular_momentum_lz(nphi, ntheta, dphi, theta, Sa);
    lzb = angular_momentum_lz(nphi, ntheta, dphi, theta, Sb);

    strcpy(fname, "output/");
    strcat(fname, prefix);
    strcat(fname, "_2species_obs_imagtime.dat");

    if (line > 1) txt_file_ptr = fopen(fname, "a");
    else          txt_file_ptr = fopen(fname, "w");

    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    fprintf(txt_file_ptr,
            "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
            energy, kin_energy, mu_a, mu_b, den_overlap, lza, lzb
    );

    free(abs_square_b);
    free(abs_square_a);
    fclose(txt_file_ptr);
}


void save_obs_1species(char prefix [], EqDataPkg EQ, Carray S, int line)
{
    int
        nphi,
        ntheta;
    double
        dphi,
        lz,
        mu,
        energy,
        kin_energy;
    Rarray
        theta,
        abs_square;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    nphi = EQ->nphi;
    ntheta = EQ->ntheta;
    dphi = EQ->dphi;
    theta = EQ->theta;

    abs_square = rarrDef(nphi * ntheta);
    carrAbs2(nphi * ntheta, S, abs_square);

    energy = functionals_single(EQ, S, &kin_energy, &mu);
    lz = angular_momentum_lz(nphi, ntheta, dphi, theta, S);

    strcpy(fname, "output/");
    strcat(fname, prefix);
    strcat(fname, "_1species_obs_imagtime.dat");

    if (line > 1) txt_file_ptr = fopen(fname, "a");
    else          txt_file_ptr = fopen(fname, "w");

    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }

    fprintf(txt_file_ptr,
            "%.10lf %.10lf %.10lf %.10lf\n",
            energy, kin_energy, mu, lz
    );

    free(abs_square);
    fclose(txt_file_ptr);
}


EqDataPkg setup_equation(char prefix [], int num_species, int line)
{

/** Read line by line of _domain file and _eq to setup the equation **/

    int
        scanf_returned,
        theta_grid_size,
        phi_grid_size,
        time_steps;
    double
        step_size,
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
            txt_file_ptr, "%lf %lf %lf %lf %lf %lf %lf",
            &nabla_coef, &rotation_freq, &frac_a, &frac_b, &ga, &gb, &gab
    );

    fclose(txt_file_ptr);

    if (scanf_returned != 7)
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

    scanf_returned = fscanf(
            txt_file_ptr,
            "%d %d %lf %d",
            &phi_grid_size, &theta_grid_size, &step_size, &time_steps
    );

    fclose(txt_file_ptr);

    if (scanf_returned != 4)
    {
        printf("\n\nWrong number of parameters or bad format in %s", fname);
        exit(EXIT_FAILURE);
    }

    if (num_species == 1)
    {
        frac_b = 0.0;
        frac_a = 1.0;
    }

    return equation_structure(
            phi_grid_size, theta_grid_size, step_size, time_steps,
            nabla_coef, rotation_freq, frac_a, frac_b, ga, gb, gab
    );

}


void setup_initial_condition(
        EqDataPkg EQ, char prefix [], Carray state_a, Carray state_b)
{
    int
        i,
        nphi,
        ntheta,
        fscanf_val;
    double
        real,
        imag;
    Rarray
        abs2;
    char
        fname[100];
    FILE
        * txt_file_ptr;

    nphi = EQ->nphi;
    ntheta = EQ->ntheta;

    abs2 = rarrDef(nphi * ntheta);

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

    for (i = 0; i < nphi * ntheta; i++)
    {
        fscanf_val = fscanf(txt_file_ptr, " (%lf%lfj)", &real, &imag);
        if (fscanf_val != 2)
        {
            printf("\n\nProblem reading initial condition data: ");
            printf("%d data points read successfully\n\n", i);
            exit(EXIT_FAILURE);
        }
        state_a[i] = real + I * imag;
        abs2[i] = real*real + imag*imag;
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

    for (i = 0; i < nphi * ntheta; i++)
    {
        fscanf_val = fscanf(txt_file_ptr, " (%lf%lfj)", &real, &imag);
        if (fscanf_val != 2)
        {
            printf("\n\nProblem reading initial condition data: ");
            printf("%d data points read successfully\n\n", i);
            exit(EXIT_FAILURE);
        }
        state_b[i] = real + I * imag;
        abs2[i] = real*real + imag*imag;
    }
    fclose(txt_file_ptr);

    free(abs2);
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
        num_species;
    double
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
    Carray
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

    if (num_species < 1 || num_species > 2)
    {
        printf("\nInvalid number of species %d in job.conf\n\n", num_species);
        exit(EXIT_FAILURE);
    }

    printf("\n\nImaginary time propagation for %d species", num_species);

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
        printf("\nInitiating job %d ...", i);

        printf("\n\n");
        printf("\t\t**********************************************\n");
        printf("\t\t*                                            *\n");
        printf("\t\t*           CONFIGURING FROM FILES           *\n");
        printf("\t\t*                                            *\n");
        printf("\t\t**********************************************\n");

        // PACK DOMAIN AND PARAMETERS INFORMATION IN A STRUCTURE
        EQ = setup_equation(infname, num_species, i);

        save_equation_setup(outfname, EQ, num_species, i);

        // configure initial condition
        Sa = carrDef(EQ->nphi * EQ->ntheta);
        Sb = carrDef(EQ->nphi * EQ->ntheta);
        setup_initial_condition(EQ, infname, Sa, Sb);

        printf("\n\n\n\n");
        printf("\t\t*********************************************\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*            GRID SPECIFICATIONS            *\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*********************************************\n");

        printf("\nphi = [ %.2lf , %.2lf , ... , %.2lf , %.2lf ]\n",
                EQ->phi[0], EQ->phi[1],
                EQ->phi[EQ->nphi - 2], EQ->phi[EQ->nphi - 1]
        );
        printf("%d grid points with spacing = %.3lf\n",EQ->nphi, EQ->dphi);

        printf("\ntheta = [ %.2lf , %.2lf , ... , %.2lf , %.2lf ]\n",
                EQ->theta[0], EQ->theta[1],
                EQ->theta[EQ->ntheta-2], EQ->theta[EQ->ntheta-1]
        );
        printf("%d grid points with spacing = %.3lf\n",
                EQ->ntheta, EQ->dtheta
        );
        printf("\n%d time steps of size %.6lf, final time = %.2lf\n",
                EQ->nt, EQ->dt, EQ->nt*EQ->dt
        );

        printf("\n\n\n\n");
        printf("\t\t*********************************************\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*         START INTEGRATION ROUTINE         *\n");
        printf("\t\t*                                           *\n");
        printf("\t\t*********************************************\n");

        start = omp_get_wtime();

        stationaryNewton(EQ, Sa, Sb, 5E-3, 30);

        time_used = (double) (omp_get_wtime() - start);
        printf("\n\nTime elapsed");
        printf(" : %.0lf sec = ",time_used);
        TimePrint(time_used);

        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_speciesA_job");
        strcat(fname, line_str);
        strcat(fname, "_imagtime.dat");
        carr_txt(fname, EQ->nphi * EQ->ntheta, Sa);
        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_speciesB_job");
        strcat(fname, line_str);
        strcat(fname, "_imagtime.dat");
        carr_txt(fname, EQ->nphi * EQ->ntheta, Sb);

        free(Sa);
        free(Sb);
        ReleaseEqDataPkg(EQ);
    }

    printf("\n\n");
    return 0;
}
