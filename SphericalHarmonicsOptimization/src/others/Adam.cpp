#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
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
        grad[i] = -dfdx[i];
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
  results_file.close();
  return rho_s;
}

double testfunc(std::vector<double> &x, std::vector<double> &grad){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt",std::ios::app);
  ++fcount;
  std::cout << "iteratioins: " << fcount << std::endl;
  
  double f=0;
  double ferr = (rand()%100);
  double sigma = 10.0;
  ferr /= 1000.0;
  for (int i=0;i<Nc-1;++i)
    f += x[i]*x[i];
  f = exp(-sigma*f)*(1+ferr-0.05);
  if (!grad.empty()) {
    for (int i=0;i<Nc-1;++i)                        // gradient 
        grad[i] = 2*x[i]*f*sigma;
  }
  /* print out temporary information *************************************/
  results_file << f << " ";
  std::cout << std::endl << "x: ";
  for (int i=0;i<Nc-1;++i){
    std::cout << x[i] << " ";
    results_file << x[i] << " ";
  }
  results_file << "\n";
  std::cout << std::endl;
  std::cout << "f: " << f << std::endl << std::endl;
  results_file.close();
  return f;
}

int main(){
  srand((unsigned)time(0));
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << "#1 rho_s \n";
  results_file << "#2 coeffs \n";
  results_file.close();
  //adam parameters
  double beta1 = 0.9;
  double beta2 = 0.999;
  double alpha = 0.01;
  double epsilon = 1e-8;
  //initial valuese
  std::vector<double> x;
  std::vector<double> m;
  double v=0;
  std::vector<double> grad;
  std::vector<double> mc;
  double vc;
  for (int i=0;i<Nc-1;++i){
    x.push_back(0);
    m.push_back(0);
    grad.push_back(0);
    mc.push_back(0);
  }
  x[2] = 1.00;
  //x[6] = 0.1;
  //x[8] = 0.1;
  for (int k=0;k<MAX_ITER;k++){
    myfunc(x,grad);
    double g2=0;
    for (int i=0;i<Nc-1;++i){
      m[i] = beta1*m[i]+(1-beta1)*grad[i];
      g2 += grad[i]*grad[i];
      mc[i] = m[i]/(1-pow(beta1,k+1));
    }
    v = beta2*v+(1-beta2)*g2;
    vc = v/(1-pow(beta2,k+1));
    for (int i=0;i<Nc-1;++i){
      x[i] = x[i] - alpha*mc[i]/(sqrt(vc)+epsilon);
    }
  }
  
  //std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
}