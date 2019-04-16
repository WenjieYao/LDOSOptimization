#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include "scuff-scatter.h"
#include "basisfunctions.h"
#include "LDOSmath.h"
#include <nlopt.hpp>
#include <vector>


#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ITER 40                            // maximum number for number of coeff sets 
#define Nc num_coeffs                           // simplified number of coefficients

int fcount = 0;
/***************************************************************/
/***************************************************************/
/***************************************************************/
double lambda_0 = 400;                          // wavelength in unit of nm
const unsigned int num_coeffs = 5*5+1;           // number of coefficients
double k_0 = 2*pi/lambda_0;                     // wavenumber in unit of 1/nm
cdouble Chi = 0;                                // chi = epsilon-1
double d_min = 20;                             // minimum distance in unit of nm
int min_mesh = 2000;
int max_mesh = 6000;
double resolution = 20.0;                        // resolution, larger the finer, default 1


double myfunc(std::vector<double> &x, std::vector<double> &grad){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt",std::ios::app);
  ++fcount;
  std::cout << "iteratioins: " << fcount << std::endl;
  
  double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0] =d_min/1e3;
  for(int i=0;i<Nc-1;++i)
    coeffs[i+1] = x[i];
  double rho_s=0;                               // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc-1];              // gradient of LDOS vs coeffs
  for (int i=0;i<Nc-1;++i)                        // gradient initialzed to 0
      dfdx[i] = 0.0;
  LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  if (!grad.empty()) {
    for (int i=0;i<Nc-1;++i)                        // gradient 
        grad[i] = dfdx[i];
  }
  /* print out temporary information *************************************/
  results_file << rho_s << " ";
  std::cout << std::endl << "x: ";
  results_file << coeffs[0] << " ";
  for (int i=0;i<Nc-1;++i){
    std::cout << x[i] << " ";
    results_file << x[i] << " ";
  }
  results_file << "\n";
  std::cout << endl;
  std::cout << "dfdx: ";
  for (int i=0;i<Nc-1;++i){
    std::cout << dfdx[i] << " ";
  }
  std::cout << endl;
  std::cout << "rho: " << rho_s << std::endl << std::endl;

  delete[] coeffs;
  delete[] dfdx;
  return rho_s;
  results_file.close();
}


int main(){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << "#1 rho_s \n";
  results_file << "#2 coeffs \n";
  results_file.close();
  std::vector<double> x0;
  std::vector<double> x;
  std::vector<double> grad0;
  std::vector<double> grad;
  double alpha = 1;
  for (int i=0;i<Nc-1;++i){
    x0.push_back(0);
    x.push_back(0);
    grad0.push_back(0);
    grad.push_back(0);
  }
  x0[0] = 1;
  //x0[5] = 1.0;
  //x0[8] = 1.0;
  double rhos = myfunc(x0,grad0);
  double ng = 0;
  for (auto i = grad0.begin(); i != grad0.cend(); ++i) 
    ng += (*i)*(*i);
  std::cout << "rhos: " << rhos << std::endl;
  std::cout << "x0: " ;
  for (auto i = x0.begin(); i != x0.cend(); ++i) 
    cout << *i << " ";
  std::cout << std::endl;
  results_file.open ("results.txt");
  results_file << "ng:" << ng ;
  results_file << "\n";
  results_file.close();
  for (int k=0;k<MAX_ITER;k++){
    for (int i=0;i<Nc-1;++i)
      x[i] = x0[i]+alpha*grad0[i]/ng;
    myfunc(x,grad);
    alpha += 1;
  }
  
  //std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  std::cout << "ng:" << ng <<std::endl;
  
}