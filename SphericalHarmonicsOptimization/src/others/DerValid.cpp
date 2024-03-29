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
#define MAX_ITER 10                            // maximum number for number of coeff sets 
#define Nc num_coeffs                           // simplified number of coefficients

int fcount = 0;
/***************************************************************/
/***************************************************************/
/***************************************************************/
double lambda_0 = 500;                          // wavelength in unit of nm
const unsigned int num_coeffs = 4*4+1;           // number of coefficients
const unsigned int Nx =10;
double k_0 = 2*pi/lambda_0;                     // wavenumber in unit of 1/nm
cdouble Chi = 0;                                // chi = epsilon-1
double d_min = 50;                             // minimum distance in unit of nm
const double d_min_c = d_min;
int min_mesh = 4000;
int max_mesh = 5000;
double resolution = 10.0;                        // resolution, larger the finer, default 1


double myfunc(const std::vector<double> &x, std::vector<double> &grad){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt",std::ios::app);
  ++fcount;
  std::cout << "iteratioins: " << fcount << std::endl;
  
  int cx[Nx];
  int xcount=0;
  for(int i=0;i<Nc-1;++i){                       
    int li=floor(sqrt(i));
    int mi=i-li*li-li;
    if((0==0)&&((li%2)==1)){
      cx[xcount]=i;
      xcount++;
    }
  }
  double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0] =d_min_c/1e3;
  for(int i=0;i<Nc-1;++i)
    coeffs[i+1] = 0;
  for(int i=0;i<Nx;++i)
    coeffs[cx[i]+1] = x[i];
  double rho_s=0;                               // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc-1];              // gradient of LDOS vs coeffs
  for (int i=0;i<Nc-1;++i)                        // gradient initialzed to 0
      dfdx[i] = 0.0;
  LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  if (!grad.empty()) {
    for (int i=0;i<Nx;++i)                        // gradient 
        grad[i] = dfdx[cx[i]];
  }
  /* print out temporary information *************************************/
  results_file << rho_s << " ";
  results_file << Nc << " ";
  std::cout << std::endl << "coeffs: ";
  //results_file << coeffs[0] << " ";
  for (int i=0;i<Nc;++i){
    std::cout << coeffs[i] << " ";
    results_file << coeffs[i] << " ";
  }
  std::cout << endl;
  std::cout << "dfdx: ";
  for (int i=0;i<Nc-1;++i){
    std::cout << dfdx[i] << " ";
    results_file << dfdx[i] << " ";
  }
  results_file << "\n";
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
  for (int i=0;i<Nx;++i){
    x0.push_back(0);
    x.push_back(0);
    grad0.push_back(0);
    grad.push_back(0);
  }
  x0[3] = 1.0;
  //x0[0] = 0.0;
  //x0[8] = 1.0;
  double rhos = myfunc(x0,grad0);
  double ng = 0;
  double step = 0.01;//rhos/100;
  double alpha = step;
  for (int i=0; i < Nx ; ++i) 
    ng += grad0[i]*grad0[i];
  std::cout << "ng:" << ng <<std::endl;
  std::cout << "step:" << step << std::endl;
  std::cout << "rhos: " << rhos << std::endl;
  std::cout << "x0: " ;
  for (int i=0; i < Nx ; ++i) 
    cout << x0[i] << " ";
  std::cout << std::endl;
  for (int k=0;k<MAX_ITER;k++){
    for (int i=0;i<Nc-1;++i)
      x[i] = x0[i]+alpha*grad0[i]/ng;
    myfunc(x,grad);
    alpha += step;
  }
  
  //std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  
  results_file.open ("results.txt",std::ios::app);
  results_file << "ng:" << ng ;
  results_file << "\n";
  results_file << "step:" << step ;
  results_file << "\n";
  results_file << "rho_limit:" << rho_limit ;
  results_file << "\n";
  results_file << "d_min:" << d_min;
  results_file << "\n";
  results_file.close();
}