#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <string.h>
#include <iomanip>
#include "scuff-scatter.h"
#include "basisfunctions.h"

#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ITER 200                            // maximum number for number of coeff sets 
#define Nc num_coeffs                           // simplified number of coefficients

int main(){
    double lambda_0 = 500;                          // wavelength in unit of nm
    const unsigned int num_coeffs = 4*4+1;           // number of coefficients
    const unsigned int Nx =10;
    double k_0 = 2*pi/lambda_0;                     // wavenumber in unit of 1/nm
    cdouble Chi = 0;                                // chi = epsilon-1
    double d_min = 50;                             // minimum distance in unit of nm
    
    std::vector<double> x;
    for (int i=0;i<Nx;++i)
    x.push_back(0);
    std::ifstream file("x_initial.txt");
    for(int i=0;i<Nx;++i)
      file >> x[i];
    file.close();
    
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
    coeffs[0] =d_min/1e3;
    for(int i=0;i<Nc-1;++i)
        coeffs[i+1] = 0;
    for(int i=0;i<Nx;++i)
        coeffs[cx[i]+1] = x[i];

    /***************************************************************/
    /* Initialze Spherical Harmonic funciton and Create Mesh *******/
    /***************************************************************/
    geometry::SphericalHarmonics *sh;  //create spherical harmonic class
    sh=new geometry::SphericalHarmonics(num_coeffs, coeffs, -1e6, 1e6);
    int Nt = 180;
    int Np = 360;
    double dt = pi/Nt;
    double dp = 2*pi/Np;
    double fdiff = 0;
    double forigin = 0;
    for(int it=0;it<Nt;++it){
        for(int ip=0;ip<Np;++ip){
            double theta = it*dt;
            double phi = ip*dp;
            double d_temp = sh->radius(theta,phi)*1e3;
            fdiff += 1.0/pow(k_0*d_temp,3) + 1.0/(k_0*d_temp);
            forigin += 1.0/pow(k_0*d_min,3) + 1.0/(k_0*d_min);
        }
    }
    std::cout << "difference by percertange: " << fdiff/forigin*100 << std::endl;
    
    return 0;
}