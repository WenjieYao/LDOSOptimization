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
#define MAX_ITER 400                            // maximum number for number of coeff sets 
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
int max_mesh = 5000;
double resolution = 5.0;                        // resolution, larger the finer, default 1


double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
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

double SpheroidRadius(double theta, double r1, double r2){
  return r1*r2/sqrt(r1*r1*cos(theta)*cos(theta)+r2*r2*sin(theta)*sin(theta));
}

void SpheroidToHarmonics(double *coeff, unsigned int num_coeff, double r1, double r2){
  geometry::SphericalHarmonics *sh;  //create spherical harmonic class
  double *temp_coeffs = new double[num_coeff];
  for(int i=0;i<num_coeff;++i)
    temp_coeffs[i] = 1;
  sh=new geometry::SphericalHarmonics(num_coeff, temp_coeffs, -1e6, 1e6);
  double *l=new double[num_coeff];
  double *m=new double[num_coeff];
  sh->get_lm(l,m);
  //create integral theta and phi
  unsigned int Nt = 500;
  unsigned int Np = 2*Nt;
  double *Theta = new double[Nt];
  double *Phi = new double[Np];
  double dt = 1.0*pi/(Nt-1);
  double dp = 2.0*pi/(Np-1);
  for(int i=0;i<Nt;++i)
    Theta[i] = dt*i;
  for(int i=0;i<Np;++i)
    Phi[i] = dp*i;
  //loop for integral
  for(int n=0;n<num_coeff;++n){
    //std::cout << n << " ";
    double cn = 0;
    for(int p=0;p<Np;++p){
      for(int t=0;t<Nt;++t){
        double theta = Theta[t];
        double phi  = Phi[p];
        cn += (1+(m[n]==0)*1)*0.5*sqrt(SpheroidRadius(theta,r1,r2)-d_min/1e3)*sh->Ylm(l[n],m[n],theta,phi)*sin(theta)*dt*dp;
        //cn += sh->Ylm(l[n+1],m[n+1],theta,phi)*sh->Ylm(l[n],m[n],theta,phi)*sin(theta)*dt*dp;
      }
    }
    //std::cout << cn << std::endl;
    coeff[n] = cn;
  }
  delete[] Theta;
  delete[] Phi;
  delete[] l;
  delete[] m;
  delete sh;
}

int main(){
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << "#1 rho_s \n";
  results_file << "#2 coeffs \n";
  results_file.close();
  nlopt::opt opt(nlopt::LD_MMA, Nc-1);

  std::vector<double> lb(Nc-1);
   std::vector<double> ub(Nc-1);
  for (int i=0;i<Nc-1;++i){
    lb[i] = -5;
    ub[i] = 5;
  }
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  opt.set_max_objective(myfunc, NULL);

  opt.set_maxeval(MAX_ITER);
  //opt.set_xtol_abs(1e-5);
  //opt.set_xtol_rel(0);
  //opt.get_maxeval();
 std::cout << "Start converting spheroid to spherical harmonics" << std::endl;
  std::vector<double> x(Nc-1);
  double *coeff = new double[Nc-1];
  double r1 = 128/1e3;
  double r2 = 128/1e3;
  for(int i=0;i<Nc-1;++i)
    coeff[i] = 0;
  coeff[4] = 1.0;
  coeff[8] = 1.0;
  //SpheroidToHarmonics(coeff,Nc-1,r1,r2);
  //x[0]=d_min/1e3*4*pi;                          // initial value 
  for(int i=0;i<Nc-1;++i){
    x[i] = coeff[i];
  }
  double maxf;
  std::cout << "Starting optimization... " << std::endl;
  nlopt::result result = opt.optimize(x, maxf);

  std::cout << std::endl << "MAX_ITER: " << opt.get_maxeval() << std::endl;
  std::cout << std::endl << "x: ";
  for (int i=0;i<Nc-1;++i){
    std::cout << x[i] << " ";
  }
  std::cout << endl;
  std::cout << "rho: " << maxf << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  delete coeff;
}