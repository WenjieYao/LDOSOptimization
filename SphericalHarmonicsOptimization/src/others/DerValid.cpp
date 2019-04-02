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
#define MAX_ITER1 1                             // maximum number for number of coeff sets
#define MAX_ITER2 50                             // maximum number for number of coeff sets
#define MAX_ITER MAX_ITER1*MAX_ITER2                             // maximum number for number of coeff sets
#define Nc num_coeffs
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  double lambda_0 = 550;                        // wavelength in unit of nm
  double k_0 = 2*pi/lambda_0;                   // wavenumber in unit of 1/nm
  cdouble Chi = 0;                              // chi = epsilon-1
  double d_min = 10;                           // minimum distance in unit of nm
  
  int min_mesh = 1000;
  int max_mesh = 5000;
  double r1 = d_min/1e3;
  double r2 = 0.0;
  double delta_c1 = 0.0;
  double delta_c2 = 0.05;                     // increasement of coefficients used

  double resolution = 1.0;                      // resolution, larger the finer, default 1
  
  /****************Loop for different coeff LDOS ***************************/
  double Rho_s[MAX_ITER1][MAX_ITER2];                       // array to store all elecric LDOS
  double Dfdr1[MAX_ITER1][MAX_ITER2];                       // matrix to store all derivatives
  double Dfdr2[MAX_ITER1][MAX_ITER2];
  double R1[MAX_ITER1][MAX_ITER2];                          // matrix to store all coefficients
  double R2[MAX_ITER1][MAX_ITER2];

  for(int k1=0;k1<MAX_ITER1;++k1){
    for(int k2=0;k2<MAX_ITER2;++k2){
      R1[k1][k2] = r1+k1*delta_c1;
      R2[k1][k2] = r2+k2*delta_c2;
    }
  }
  double rho_s;                                 // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc];                // gradient of LDOS vs coeffs
  for (int i=0;i<Nc;++i)                        // gradient initialzed to 0
    dfdx[i] = 0.0;
  
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << " #1 r1" << "\n";
  results_file << " #2 r2" << "\n";
  results_file << " #3 Rho_s" << "\n";
  results_file << " #4 Dfdr1_t (Theoretical gradients using SIE) " << "\n";
  results_file << " #5 Dfdr2_t (Theoretical gradients using SIE) " << "\n";
  results_file.close();
  /*********** k = 0 *******************************************************
  LDOS_gradient(lambda_0, r1, r2, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  Rho_s[0] = rho_s;
  //Dfdx(0) = dfdx[pc];
  std::cout << "Rho_s: " << Rho_s[0] << std::endl;
  std::cout << "dfdx: " << dfdx[pc] << std::endl;

  Dfdr1[0] = dfdx[0];
  Dfdr2[0] = dfdx[1];  
  /**************************************************************************/
  /*********** Main Loop ****************************************************/
  /**************************************************************************/
  for (int k1=0;k1<MAX_ITER1;++k1){
    for(int k2=0;k2<MAX_ITER2;++k2){
    std::cout << std::endl << "iteratoin k = " << (k1)*MAX_ITER1+k2 << " out of MAX_ITER: " << MAX_ITER1*MAX_ITER2 << std::endl;
    rho_s = 0;
    for (int i=0;i<Nc;++i)                      // gradient initialzed to 0
      dfdx[i] = 0.0;
    /* intialize coefficients                  *****************************/
    double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0] =d_min/1e3;
  for(int i=0;i<Nc-1;++i)
    coeffs[i+1] = x[i];
    //pc = k%Nc;
    r1 = R1[k1][k2];
    r2 = R2[k1][k2];

    /**************compute the gradient of LDOS at current coeffs **********/
    LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,min_mesh,max_mesh,resolution);
  
    Dfdr1[k1][k2] = dfdx[0];
    Dfdr2[k1][k2] = dfdx[1];

    Rho_s[k1][k2] = rho_s;
    /*
    if (pc==0)
        Dfdr_t[k-1] = (Dfdr1[k]+Dfdr1[k-1])/2;
    else
        Dfdr_t[k-1] = (Dfdr2[k]+Dfdr2[k-1])/2;

    Dfdr_s[k-1] = (Rho_s[k]-Rho_s[k-1])/delta_c;
    std::cout << "Rho_s: " << Rho_s[k] << std::endl;
    std::cout << "dfdr_t: " << Dfdr_t[k-1] << std::endl;
    std::cout << "dfdr_s: " << Dfdr_s[k-1] << std::endl;
    */
    std::cout << "Rho_s: " << Rho_s[k1][k2] << std::endl;
    ofstream results_file; 
    results_file.open ("results.txt",std::ios::app);
    results_file << r1 << " " << r2 << " " << Rho_s[k1][k2] << " " << Dfdr1[k1][k2] << " " << Dfdr2[k1][k2] << "\n";
    results_file.close();
  }
  }
  
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  /*
  results_file.open ("results.txt",std::ios::app);
  results_file << Rho_s[MAX_ITER-1] << " " << Dfdr1[MAX_ITER-1] << " " << Dfdr2[MAX_ITER-1];
  results_file.close();
  */
  delete[] dfdx;
}