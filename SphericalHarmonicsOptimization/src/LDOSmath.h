#ifndef LDOSMATH
#define LDOSMATH
void LDOS_gradient(double lambda_0, double *coeffs, const unsigned int num_coeffs, double &rho_s, double *dfdx, cdouble &Chi, double &d_min, int min_mesh, int max_mesh, double resolution=1.0);
void WriteCoeff(double* coeffs, const unsigned int num_coeffs);
#endif