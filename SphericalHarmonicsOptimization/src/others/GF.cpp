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
#define MAX_ITER 200                            // maximum number for number of coeff sets 
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


double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
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
  std::cout << std::endl << "coeffs: ";
  //results_file << coeffs[0] << " ";
  for (int i=0;i<Nc;++i){
    std::cout << coeffs[i] << " ";
    results_file << coeffs[i] << " ";
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
  nlopt::opt opt(nlopt::LN_BOBYQA, Nx);
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
  std::vector<double> lb(Nx);
  std::vector<double> ub(Nx);
  for (int i=0;i<Nx;++i){
    if(i<3){
      lb[i] = -2;
      ub[i] = 2;
    }
    else{
      lb[i] = -1;
      ub[i] = 1;
    }
  }
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  opt.set_max_objective(myfunc, NULL);

  opt.set_maxeval(MAX_ITER);
  //opt.set_xtol_abs(1e-5);
  //opt.set_xtol_rel(0);
  //opt.get_maxeval();
  //std::cout << "Start converting spheroid to spherical harmonics" << std::endl;
  std::vector<double> x(Nx);
  x[0]=1;
  double maxf;
  std::cout << "Starting optimization... " << std::endl;
  nlopt::result result = opt.optimize(x, maxf);

  std::cout << std::endl << "MAX_ITER: " << opt.get_maxeval() << std::endl;
  std::cout << std::endl << "x: ";
  for (int i=0;i<Nx;++i){
    std::cout << x[i] << " ";
  }
  std::cout << endl;
  std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
}