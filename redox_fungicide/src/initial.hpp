
//----------------------------------
// FUNCTION PROTOTYPES
//----------------------------------
double** mallocDouble(int row, int col);
void freeMallocDouble(double **data, int row);
double* mallocSingle(int col);
void freeMallocSingle(double *data);
void initialization(double **U, int Nx, int Ny, double value);
void Noise(double **U, int Nx, int Ny, double amp, double minVal, double maxVal);
void Circular(double **U, int Nx, int Ny, double dx, double dy, double cx, double cy, double radius, double insideValue, double outsideValue);
