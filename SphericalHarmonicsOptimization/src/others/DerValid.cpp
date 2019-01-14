#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>                         //Using Eigen3 for linear algebra
#include "scuff-scatter.h"
#include "basisfunctions.h"
#include "LDOSmath.h"

#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ITER 50                             // maximum number for number of coeff sets 
#define Nc num_coeffs                           // simplified number of coefficients

typedef Eigen::VectorXd VectorXd;               // double type vector
typedef Eigen::MatrixXd MatrixXd;               // double type matrix
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  double lambda_0 = 550;                        // wavelength in unit of nm
  double k_0 = 2*pi/lambda_0;                   // wavenumber in unit of 1/nm
  cdouble Chi = 0;                              // chi = epsilon-1
  double d_min = 200;                           // minimum distance in unit of nm
  const unsigned int num_coeffs = 9+1;            // number of coefficients
  
  int min_mesh = 2000;
  int max_mesh = 5000;
  double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0]=d_min/1e3;                     // initial value Ylm(0,0,0,0)^2 = 4pi
  double delta_c = 1e-2;               // increasement of coefficients used
  for(int i=1;i<Nc;++i)
    coeffs[i] = 0.5/i;
  //coeffs[1] = 0.1;
  //coeffs[2] = 0.1;
  WriteCoeff(coeffs,num_coeffs);                // write coefficients to file

  double resolution = 1.5;                      // resolution, larger the finer, default 1
  
  /****************Loop for different coeff LDOS ***************************/
  VectorXd Rho_s(MAX_ITER);                     // array to store all elecric LDOS
  MatrixXd Dfdx(Nc-1,MAX_ITER);                   // matrix to store all derivatives
  VectorXd Dfdx_t(MAX_ITER-1);
  VectorXd Dfdx_s(MAX_ITER-1);
  MatrixXd Coeffs(Nc-1,MAX_ITER);                 // matrix to store all coefficients
  double rho_s;                                 // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc-1];                // gradient of LDOS vs coeffs
  for (int i=0;i<Nc-1;++i)                        // gradient initialzed to 0
    dfdx[i] = 0.0;
  unsigned int pc = 0;
  /*********** k = 0 *******************************************************/
  LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  Rho_s(0) = rho_s;
  //Dfdx(0) = dfdx[pc];
  std::cout << "Rho_s: " << Rho_s(0) << std::endl;
  std::cout << "dfdx: " << dfdx[pc] << std::endl;
  for (int i=0;i<Nc-1;++i){
    Coeffs(i,0) = coeffs[i+1];
    Dfdx(i,0) = dfdx[i];
  }
  /**************************************************************************/
  /*********** Main Loop ****************************************************/
  /**************************************************************************/
  for (int k=1;k<MAX_ITER;++k){
    std::cout << std::endl << "iteratoin k = " << k << " out of MAX_ITER: " << MAX_ITER << std::endl;
    rho_s = 0;
    for (int i=0;i<Nc-1;++i)                      // gradient initialzed to 0
      dfdx[i] = 0.0;
    /* intialize coefficients                  *****************************/
    pc = k%(Nc-1);
    coeffs[pc+1] += delta_c;
    /**************compute the gradient of LDOS at current coeffs **********/
    LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
    for (int i=0;i<Nc-1;++i){
      Coeffs(i,k) = coeffs[i+1];
      Dfdx(i,k) = dfdx[i];
    }

    Rho_s(k) = rho_s;
    Dfdx_t(k-1) = (Dfdx(pc,k)+Dfdx(pc,k-1))/2;
    Dfdx_s(k-1) = (Rho_s(k)-Rho_s(k-1))/delta_c;
    std::cout << "Rho_s: " << Rho_s(k) << std::endl;
    std::cout << "dfdx_t: " << Dfdx_t(k-1) << std::endl;
    std::cout << "dfdx_s: " << Dfdx_s(k-1) << std::endl;

    
  }
  
  std::cout << "Rho_s: " << Rho_s.transpose() << std::endl;
  std::cout << "Dfdx_t: " << Dfdx_t.transpose() << std::endl;
  std::cout << "Dfdx_s: " << Dfdx_s.transpose() << std::endl;
  std::cout << "Coeffs: " << std::endl << Coeffs << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;

  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << " #1 Rho_s" << "\n";
  results_file << " #2 Dfdx_t (Theoretical gradients using SIE) " << "\n";
  results_file << " #3 Dfdx_s (Simulation gradients using forward diff)" << "\n";
  for(int k=0;k<MAX_ITER-1;++k){
    results_file << Rho_s(k) << " " << Dfdx_t(k) << " " << Dfdx_s(k) << "\n";
  }
  results_file << Rho_s(MAX_ITER-1);
  results_file.close();
  delete[] dfdx;
  delete[] coeffs;
}