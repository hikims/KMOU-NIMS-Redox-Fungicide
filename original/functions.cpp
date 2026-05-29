#include "functions.hpp"      
#include "parameters.h" 
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>


// --------------------------------------------------------
// LU Decomposition with Neumann Boundary Condition
// --------------------------------------------------------
void LUx(double *a, double *b, double *c, double **L, double *U, int size)
{
    int i;
    L[0][0] = a[0];

    for (i = 1; i < size; i++)  
        L[1][i] = b[i];

    for (i = 0; i < size - 1; i++)  
    {
        U[i] = c[i] / L[0][i];
        L[0][i + 1] = a[i + 1] - L[1][i + 1] * U[i];
    }
}

void solx(double **L, double *U, double *Init, double *z, double *x, int size)
{
    int i;
    z[0] = Init[0] / L[0][0];
    for (i = 1; i < size; i++)
        z[i] = (Init[i] - L[1][i] * z[i - 1]) / L[0][i];

    x[size - 1] = z[size - 1];
    for (i = size - 2; i >= 0; i--)
        x[i] = z[i] - U[i] * x[i + 1];
}

void LUy(double *a, double *b, double *c, double **L, double *U, int size)
{
    int i;
    L[0][0] = a[0];

    for (i = 1; i < size; i++)  
        L[1][i] = b[i];

    for (i = 0; i < size - 1; i++)  
    {
        U[i] = c[i] / L[0][i];
        L[0][i + 1] = a[i + 1] - L[1][i + 1] * U[i];
    }
}

void soly(double **L, double *U, double *Init, double *z, double *x, int size)
{
    int i;
    z[0] = Init[0] / L[0][0];
    for (i = 1; i < size; i++)
        z[i] = (Init[i] - L[1][i] * z[i - 1]) / L[0][i];

    x[size - 1] = z[size - 1];
    for (i = size - 2; i >= 0; i--)
        x[i] = z[i] - U[i] * x[i + 1];
}

// Hill function
double hill_f(const double R, double m, double theta) 
{
    const double Rm = pow(R, m);
    const double th = pow(theta, m);
    return Rm / (Rm + th);
}

// hill_g using constexpr beta from parameters.h (original)
double hill_g(const double R, double m, double theta) 
{
    const double Rm = pow(R, m);
    const double th = pow(theta, m);
    return beta * th / (Rm + th);
}

// hill_g overload with explicit runtime beta
double hill_g(const double R, double m, double theta, double beta_val) 
{
    const double Rm = pow(R, m);
    const double th = pow(theta, m);
    return beta_val * th / (Rm + th);
}

// diffusion term for A
double coefficient_dA(double A, double r, double x, double y, double z, double F, double alpha)
{
    return A + r * (x - 2.0 * y + z) + dt * 0.5 * (-alpha * A * F - l * A);
}

// coefficient_dF using constexpr beta (original)
double coefficient_dF(double F, double r, double x, double y, double z, double R, double m, double theta)
{
    return F + r * (x - 2.0 * y + z) + dt * 0.5 * (hill_f(R, m, theta) * F * (1 - F) - hill_g(R, m, theta) * F);
}

// coefficient_dF overload with explicit runtime beta
double coefficient_dF(double F, double r, double x, double y, double z, double R, double m, double theta, double beta_val)
{
    return F + r * (x - 2.0 * y + z) + dt * 0.5 * (hill_f(R, m, theta) * F * (1 - F) - hill_g(R, m, theta, beta_val) * F);
}
