#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//----------------------------------
// Allocation
//----------------------------------
double** mallocDouble(int row, int col)
{
    double **Mat = (double**) malloc(sizeof(double*) * row);
    if (!Mat) {
        fprintf(stderr, "Memory allocation failed (row pointers)\n");
        exit(1);
    }

    for (int i = 0; i < row; i++)
    {
        Mat[i] = (double*) malloc(sizeof(double) * col);
        if (!Mat[i]) {
            fprintf(stderr, "Memory allocation failed (row %d)\n", i);
            exit(1);
        }

        for (int j = 0; j < col; j++)
            Mat[i][j] = 0.0;
    }
    return Mat;
}

double* mallocSingle(int col)
{
    double *arr = (double*) malloc(sizeof(double) * col);
    if (!arr) {
        fprintf(stderr, "Memory allocation failed (single)\n");
        exit(1);
    }

    for (int j = 0; j < col; j++)
        arr[j] = 0.0;

    return arr;
}

void freeMallocDouble(double **data, int row)
{
    for (int i = 0; i < row; i++)
        free(data[i]);
    free(data);
}

void freeMallocSingle(double *data)
{
    free(data);
}

//----------------------------------
// Initialization
//----------------------------------
void initialization(double **U, int Nx, int Ny, double value)
{
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            U[i][j] = value;
        }
    }
}

void Noise(double **U, int Nx, int Ny, double amp, double minVal, double maxVal)
{
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double noise = amp * (((double)rand() / RAND_MAX) - 0.5);
            U[i][j] += noise;

            if (U[i][j] < minVal) U[i][j] = minVal;
            if (U[i][j] > maxVal) U[i][j] = maxVal;
        }
    }
}

void Circular(double **U, int Nx, int Ny, double dx, double dy, double cx, double cy, double radius, double insideValue, double outsideValue)
{
    double r2 = radius * radius;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double x = i * dx;
            double y = j * dy;
            double dist2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);

            if (dist2 <= r2) U[i][j] = insideValue;
            else             U[i][j] = outsideValue;
        }
    }
}
