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

#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ITER 100                             // maximum number for number of coeff sets
#define Nc 2
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  double lambda_0 = 550;                        // wavelength in unit of nm
  double k_0 = 2*pi/lambda_0;                   // wavenumber in unit of 1/nm
  cdouble Chi = 0;                              // chi = epsilon-1
  double d_min = 203.4;                           // minimum distance in unit of nm
  
  int min_mesh = 1000;
  int max_mesh = 5000;
  double r1 = d_min/1000;
  double r2 = d_min/1000;
  double delta_c = 2.0/1e3;                     // increasement of coefficients used

  double resolution = 1.5;                      // resolution, larger the finer, default 1
  
  /****************Loop for different coeff LDOS ***************************/
  double Rho_s[MAX_ITER];                       // array to store all elecric LDOS
  double Dfdr1[MAX_ITER];                       // matrix to store all derivatives
  double Dfdr2[MAX_ITER];
  double Dfdr_t[MAX_ITER-1];
  double Dfdr_s[MAX_ITER-1];
  double R1[MAX_ITER];                          // matrix to store all coefficients
  double R2[MAX_ITER];

  double rho_s;                                 // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc];                // gradient of LDOS vs coeffs
  for (int i=0;i<Nc;++i)                        // gradient initialzed to 0
    dfdx[i] = 0.0;
  unsigned int pc = 0;
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << " #1 Rho_s" << "\n";
  results_file << " #2 Dfdr1_t (Theoretical gradients using SIE) " << "\n";
  results_file << " #2 Dfdr2_t (Theoretical gradients using SIE) " << "\n";
  results_file << " #3 Dfdr_s (Simulation gradients using forward diff)" << "\n";
  results_file.close();
  /*********** k = 0 *******************************************************/
  LDOS_gradient(lambda_0, r1, r2, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  Rho_s[0] = rho_s;
  //Dfdx(0) = dfdx[pc];
  std::cout << "Rho_s: " << Rho_s[0] << std::endl;
  std::cout << "dfdx: " << dfdx[pc] << std::endl;
  
  R1[0] = r1;
  R2[0] = r2;
  Dfdr1[0] = dfdx[0];
  Dfdr2[0] = dfdx[1];  
  /**************************************************************************/
  /*********** Main Loop ****************************************************/
  /**************************************************************************/
  for (int k=1;k<MAX_ITER;++k){
    std::cout << std::endl << "iteratoin k = " << k << " out of MAX_ITER: " << MAX_ITER << std::endl;
    rho_s = 0;
    for (int i=0;i<Nc;++i)                      // gradient initialzed to 0
      dfdx[i] = 0.0;
    /* intialize coefficients                  *****************************/
    //pc = k%Nc;
    if (pc==0)
        r1 += delta_c;
    else
        r2 += delta_c;

    /**************compute the gradient of LDOS at current coeffs **********/
    LDOS_gradient(lambda_0, r1, r2, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
    
    R1[k] = r1;
    R2[k] = r2;
    Dfdr1[k] = dfdx[0];
    Dfdr2[k] = dfdx[1];

    Rho_s[k] = rho_s;
    if (pc==0)
        Dfdr_t[k-1] = (Dfdr1[k]+Dfdr1[k-1])/2;
    else
        Dfdr_t[k-1] = (Dfdr2[k]+Dfdr2[k-1])/2;

    Dfdr_s[k-1] = (Rho_s[k]-Rho_s[k-1])/delta_c;
    std::cout << "Rho_s: " << Rho_s[k] << std::endl;
    std::cout << "dfdr_t: " << Dfdr_t[k-1] << std::endl;
    std::cout << "dfdr_s: " << Dfdr_s[k-1] << std::endl;
    ofstream results_file; 
    results_file.open ("results.txt",std::ios::app);
    results_file << Rho_s[k-1] << " " << Dfdr1[k-1] << " " << Dfdr2[k-1] << " " << Dfdr_s[k-1] << "\n";
    results_file.close();
  }
  
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  results_file.open ("results.txt",std::ios::app);
  results_file << Rho_s[MAX_ITER-1] << " " << Dfdr1[MAX_ITER-1] << " " << Dfdr2[MAX_ITER-1];
  results_file.close();
  delete[] dfdx;
}