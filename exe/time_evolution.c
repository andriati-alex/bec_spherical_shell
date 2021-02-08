#include <string.h>
#include <stdio.h>
#include <math.h>
#include "imagtimeIntegrator.h"


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


void save_equation_setup(FILE * f, EqDataPkg EQ)
{
    // Grid domain
    fprintf(f, "%d %d %lf %d ", EQ->nphi, EQ->ntheta, EQ->dt, EQ->nt);
    // Equation parameters
    fprintf(
            f,
            "%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",
            EQ->nabla_coef, EQ->omega, EQ->frac_a, EQ->frac_b,
            EQ->ga, EQ->gb, EQ->gab
    );
}


EqDataPkg setup_equation(char prefix [])
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

    int
        N,
        i,
        j,
        k;
    double
        start,      // start trigger to measure time
        time_used;  // Time used in calling evolution routine
    char
        c,
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
                fscanf(txt_file_ptr, "%s", infname);
                i = i + 1;
                break;
            case 2:
                fscanf(txt_file_ptr, "%s", outfname);
                break;
        }
        ReachNewLine(txt_file_ptr);
    }
    if (i != 2)
    {
        printf("\n\nWrong number of parameter read from job.conf file\n\n");
        exit(EXIT_FAILURE);
    }

    fclose(txt_file_ptr);


    printf("\n\n");
    printf("\t\t**********************************************\n");
    printf("\t\t*                                            *\n");
    printf("\t\t*           CONFIGURING FROM FILES           *\n");
    printf("\t\t*                                            *\n");
    printf("\t\t**********************************************\n");

    // PACK DOMAIN AND PARAMETERS INFORMATION IN A STRUCTURE
    EQ = setup_equation(infname);

    // configure initial condition
    Sa = carrDef(EQ->nphi * EQ->ntheta);
    Sb = carrDef(EQ->nphi * EQ->ntheta);
    setup_initial_condition(EQ, infname, Sa, Sb);

    // save equation parameters in output file
    strcpy(fname, "output/");
    strcat(fname, outfname);
    strcat(fname, "_equation_imagtime.dat");
    txt_file_ptr = fopen(fname, "w");
    if (txt_file_ptr == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n\n", fname);
        exit(EXIT_FAILURE);
    }
    save_equation_setup(txt_file_ptr, EQ);
    fclose(txt_file_ptr);

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
    printf("%d grid points with spacing = %.3lf\n", EQ->ntheta, EQ->dtheta);
    printf("\nFinal time %.2lf in steps of %.6lf\n", EQ->nt*EQ->dt, EQ->dt);

    /*  ===============================================================
     
                             CALL INTEGRATION ROUTINE
     
        ===============================================================  */

    N = splitstep_spherical_shell(EQ, Sa, Sb);

    // Record data
    strcpy(fname, "output/");
    strcat(fname, outfname);
    strcat(fname, "_final_speciesA_imagtime.dat");
    carr_txt(fname, EQ->nphi * EQ->ntheta, Sa);
    strcpy(fname, "output/");
    strcat(fname, outfname);
    strcat(fname, "_final_speciesB_imagtime.dat");
    carr_txt(fname, EQ->nphi * EQ->ntheta, Sb);

    free(Sa);
    free(Sb);
    ReleaseEqDataPkg(EQ);

    printf("\n\n");
    return 0;
}
