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
#define MAX_ITER 20                            // maximum number for number of coeff sets 
#define Nc num_coeffs                           // simplified number of coefficients

int fcount = 0;
/***************************************************************/
/***************************************************************/
/***************************************************************/
double lambda_0 = 500;                          // wavelength in unit of nm
const unsigned int num_coeffs = 5*5+1;           // number of coefficients
double k_0 = 2*pi/lambda_0;                     // wavenumber in unit of 1/nm
cdouble Chi = 0;                                // chi = epsilon-1
double d_min = 20;                             // minimum distance in unit of nm
const double d_min_c = d_min;
int min_mesh = 4000;
int max_mesh = 6000;
double resolution = 15.0;                        // resolution, larger the finer, default 1
int deriv = 7;

double myfunc(std::vector<double> &x, std::vector<double> &grad){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt",std::ios::app);
  ++fcount;
  std::cout << "iteratioins: " << fcount << std::endl;
  
  double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0] =d_min_c/1e3;
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
  results_file << dfdx[deriv] << " ";
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
  for (int i=0;i<Nc-1;++i){
    x0.push_back(0);
    x.push_back(0);
    grad0.push_back(0);
    grad.push_back(0);
  }
  // intial guess
  if(true){
    std::ifstream file("coeff_initial.txt");
    if (file.is_open()){
    for(int i=0;i<num_coeffs-1;++i){
      file >> x0[i];
    }
    file.close();
    }
  }
  else{
    for(int i=1;i<num_coeffs-1;++i)
      x0[i] = 0.2/pow(i-3.5,2);
  }
  double rhos = 0;//myfunc(x0,grad0);
  double step = 0.01;//rhos/100;
  double alpha = step;
  std::cout << "step:" << step << std::endl;
  std::cout << "rhos: " << rhos << std::endl;
  std::cout << "x0: " ;
  for (int i=0; i < Nc -1 ; ++i) 
    cout << x0[i] << " ";
  std::cout << std::endl;
  for (int k=0;k<MAX_ITER;k++){
    for (int i=0;i<Nc-1;++i)
      x[i] = x0[i]+alpha*(i==deriv);
    myfunc(x,grad);
    alpha += step;
  }
  
  //std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  
  results_file.open ("results.txt",std::ios::app);
  results_file << "step:" << step ;
  results_file << "\n";
  results_file << "rho_limit:" << rho_limit ;
  results_file << "\n";
  results_file << "d_min:" << d_min;
  results_file << "\n";
  results_file.close();
}