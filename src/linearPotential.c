#include "linearPotential.h"



void harmonic(int nx, int ny, Rarray x, Rarray y, Rarray V, double wx, double wy)
{
    int
        i,
        j;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            V[i + j*nx] = 0.5 * (wx*wx*x[i]*x[i] + wy*wy*y[j]*y[j]);
        }
    }
}



void quartic(int nx, int ny, Rarray x, Rarray y, Rarray V, double wx, double wy)
{
    int
        i,
        j;

    double
        x4,
        y4;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x4 = x[i] * x[i] * x[i] * x[i];
            y4 = y[j] * y[j] * y[j] * y[j];
            V[i + j*nx] = 0.5 * (wx*wx*x4 + wy*wy*y4);
        }
    }
}



void quarticXQuadratic(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wy, double a, double b)
{
    int
        i,
        j;

    double
        a4,
        a2,
        b2,
        x2,
        x4,
        y2;

    a2 = a * a;
    b2 = b * b;
    a4 = a2 * a2;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x4 = x[i] * x[i] * x[i] * x[i];
            x2 = x[i] * x[i];
            y2 = y[j] * y[j];
            V[i + j*nx] = 0.5*(-a2*x2 + wy*wy*y2) + 0.25*b2*x4 + 0.25*a4/b2;
        }
    }
}



void quarticRQuadratic(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wy, double a, double b)
{
    int
        i,
        j;

    double
        a4,
        a2,
        b2,
        x2,
        r4,
        y2;

    a2 = a * a;
    b2 = b * b;
    a4 = a2 * a2;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x2 = x[i] * x[i];
            y2 = y[j] * y[j];
            r4 = (x2 + y2) * (x2 + y2);
            V[i + j*nx] = 0.5*(-a2*x2 + wy*wy*y2) + 0.25*b2*r4 + 0.25*a4/b2;
        }
    }
}



void HarmonicMexicanHat(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wx, double wy, double height, double width, double asym)
{
    int
        i,
        j;

    double
        x2,
        y2,
        y2asym;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x2 = x[i] * x[i];
            y2 = y[j] * y[j];
            y2asym = asym * asym * y2;
            V[i + j*nx] = 0.5 * (wx*wx*x2 + wy*wy*y2) + \
                          height * exp(-(x2 + y2asym) / width / width);
        }
    }
}



void QuarticMexicanHat(int nx, int ny, Rarray x, Rarray y, Rarray V,
     double wx, double wy, double height, double width, double asym)
{
    int
        i,
        j;

    double
        x2,
        y2,
        y2asym;

    for (i = 0; i < nx; i ++)
    {
        for (j = 0; j < ny; j++)
        {
            x2 = x[i] * x[i];
            y2 = y[j] * y[j];
            y2asym = asym * asym * y2;
            V[i + j*nx] = 0.5 * (wx*wx*x2*x2 + wy*wy*y2*y2) + \
                          height * exp(-(x2 + y2asym) / width / width);
        }
    }
}



void GetPotential(char name [], int nx, int ny, Rarray x, Rarray y, Rarray V,
                  double p [])
{

    if (strcmp(name, "zero") == 0)
    {
        rarrFill(nx*ny,0,V);
        return;
    }

    if (strcmp(name, "harmonic") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0)
        {
            printf("\n\nHarmonic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        harmonic(nx,ny,x,y,V,p[0],p[1]);
        return;
    }

    if (strcmp(name, "quartic") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0)
        {
            printf("\n\nQuartic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        quartic(nx,ny,x,y,V,p[0],p[1]);
        return;
    }

    if (strcmp(name, "quarticXQuadratic") == 0)
    {
        if (p[0] <= 0 || p[2] <= 0)
        {
            printf("\n\nQuartic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        quarticXQuadratic(nx,ny,x,y,V,p[0],p[1],p[2]);
        return;
    }

    if (strcmp(name, "quarticRQuadratic") == 0)
    {
        if (p[2] <= 0)
        {
            printf("\n\nQuartic Trap parameters must be positive\n");
            exit(EXIT_FAILURE);
        }

        quarticRQuadratic(nx,ny,x,y,V,p[0],p[1],p[2]);
        return;
    }

    if (strcmp(name, "HarmonicMexicanHat") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0 || p[2] <= 0)
        {
            printf("\n\nHarmonic parameters must be positive and");
            printf(" gaussian height must be positive.\n");
            exit(EXIT_FAILURE);
        }

        HarmonicMexicanHat(nx,ny,x,y,V,p[0],p[1],p[2],p[3],p[4]);
        return;
    }

    if (strcmp(name, "QuarticMexicanHat") == 0)
    {
        if (p[0] <= 0 || p[1] <= 0 || p[2] <= 0)
        {
            printf("\n\nHarmonic parameters must be positive and");
            printf(" gaussian height must be positive.\n");
            exit(EXIT_FAILURE);
        }

        QuarticMexicanHat(nx,ny,x,y,V,p[0],p[1],p[2],p[3],p[4]);
        return;
    }

    printf("\n\n\nERROR: Potential '%s' not implemented.", name);
    printf("\nHave a look in src/linearPotential.c to see");
    printf(" the available ones or to define a new one.\n\n");
    exit(EXIT_FAILURE);
}
