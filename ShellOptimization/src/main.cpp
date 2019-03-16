#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include "scuff-scatter.h"
#include "LDOSmath.h"

#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ITER1 30                             // maximum number for number of coeff sets
#define MAX_ITER2 1                             // maximum number for number of coeff sets
#define MAX_ITER MAX_ITER1*MAX_ITER2                             // maximum number for number of coeff sets
#define Nc 2
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  double lambda_0 = 550;                        // wavelength in unit of nm
  double k_0 = 2*pi/lambda_0;                   // wavenumber in unit of 1/nm
  cdouble Chi = 0;                              // chi = epsilon-1
  
  int min_mesh = 1000;
  int max_mesh = 5000;
  double r1 = 210.0/1e3;
  double r2 = 203.5/1e3;
  double delta_c = 5.0/1e3;                     // increasement of coefficients used

  double resolution = 1.5;                      // resolution, larger the finer, default 1
  
  /****************Loop for different coeff LDOS ***************************/
  double Rho_s[MAX_ITER1][MAX_ITER2];                       // array to store all elecric LDOS
  double R1[MAX_ITER1][MAX_ITER2];                          // matrix to store all coefficients
  double R2[MAX_ITER1][MAX_ITER2];

  for(int k1=0;k1<MAX_ITER1;++k1){
    for(int k2=0;k2<MAX_ITER2;++k2){
      R1[k1][k2] = r1+k1*delta_c;
      R2[k1][k2] = r2+k2*delta_c;
    }
  }
  double rho_s;                                 // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc];                // gradient of LDOS vs coeffs
  
  ofstream results_file;                          //write to file
  results_file.open ("results.txt");
  results_file << " #1 r1" << "\n";
  results_file << " #2 r2" << "\n";
  results_file << " #3 Rho_s" << "\n";
  results_file.close();
  /**************************************************************************/
  /*********** Main Loop ****************************************************/
  /**************************************************************************/
  for (int k1=0;k1<MAX_ITER1;++k1){
    for(int k2=0;k2<MAX_ITER2;++k2){
    std::cout << std::endl << "iteratoin k = " << (k1)*MAX_ITER2+k2 << " out of MAX_ITER: " << MAX_ITER1*MAX_ITER2 << std::endl;
    rho_s = 0;
    /* intialize coefficients                  *****************************/
    //pc = k%Nc;
    r1 = R1[k1][k2];
    r2 = R2[k1][k2];

    /**************compute the gradient of LDOS at current coeffs **********/
    LDOS_eval(lambda_0, r1, r2, rho_s ,min_mesh,max_mesh,resolution);

    Rho_s[k1][k2] = rho_s;

    std::cout << "Rho_s: " << Rho_s[k1][k2] << std::endl;
    ofstream results_file; 
    results_file.open ("results.txt",std::ios::app);
    results_file << r1 << " " << r2 << " " << Rho_s[k1][k2] << "\n";
    results_file.close();
  }
  }
}