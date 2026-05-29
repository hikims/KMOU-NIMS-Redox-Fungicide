
#include <vector>

//---------------------------
// FUNCTION PROTOTYPE
//---------------------------
struct MergedClusters {
    std::vector<double> cx;
    std::vector<double> cy;
    std::vector<int>    size;
};

void LUx(double *a, double *b, double *c, double **L, double *U, int size);
void solx(double **L, double *U, double *Init, double *z, double *x, int size);
void LUy(double *a, double *b, double *c, double **L, double *U, int size);
void soly(double **L, double *U, double *Init, double *z, double *x, int size);

double hill_f(const double R, double m, double theta);

// Original hill_g uses constexpr beta from parameters.h
double hill_g(const double R, double m, double theta);
// Overload: explicit runtime beta
double hill_g(const double R, double m, double theta, double beta_val);

double coefficient_dA(double A, double r, double x, double y, double z, double F, double alpha);

// Original coefficient_dF uses constexpr beta from parameters.h
double coefficient_dF(double F, double r, double x, double y, double z, double R, double m, double theta);
// Overload: explicit runtime beta
double coefficient_dF(double F, double r, double x, double y, double z, double R, double m, double theta, double beta_val);
