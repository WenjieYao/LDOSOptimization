#ifndef LDOSMATH
#define LDOSMATH
#include <Eigen/Dense>                         //Using Eigen3 for linear algebra
void LDOS_gradient(double lambda_0, double *coeffs, const unsigned int num_coeffs, double &rho_s, double *dfdx, cdouble &Chi, double &d_min, int min_mesh, int max_mesh, double resolution=1.0);
void LDOS_gradient(double lambda_0, Eigen::VectorXd Coeffs, const unsigned int num_coeffs, double &rho_s,Eigen::VectorXd &Dfdx, cdouble &Chi, double &d_min, double resolution=1.0);
void WriteCoeff(double* coeffs, const unsigned int num_coeffs);

//double* Mult(double **a,double *b,unsigned int n);
//double* Mult(double *a,double **b,unsigned int n);
#endif