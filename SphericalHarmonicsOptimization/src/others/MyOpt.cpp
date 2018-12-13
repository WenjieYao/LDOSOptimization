#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>                         //Using Eigen3 for linear algebra
#include "scuff-scatter.h"
#include "basisfunctions.h"
#include "LDOSmath.h"

#define pi 3.1415926535897  
#define Z_0 376.730313                          // Free space impedance, unit : Ohms
#define MAX_ALPHA 10                            // maximum number for alpha search
#define MAX_ITER 1                            // maximum number for shape optimization 
#define Nc num_coeffs                           // simplified number of coefficients

typedef Eigen::VectorXd VectorXd;               // double type vector
typedef Eigen::MatrixXd MatrixXd;               // double type matrix
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{
  double lambda_0 = 550;                        // wavelength in unit of nm
  double k_0 = 2*pi/lambda_0;                   // wavenumber in unit of 1/nm
  cdouble Chi = 0;                              // chi = epsilon-1
  double d_min = 200;                           // minimum distance in unit of nm
  const unsigned int num_coeffs = 9;            // number of coefficients
  
  double *coeffs = new double[Nc];              // array to store coeffs
  coeffs[0]=d_min/1e3*2*sqrt(pi);               // initial value 
  for(int i=1;i<Nc;++i)
    coeffs[i] = 0.00;
  WriteCoeff(coeffs,num_coeffs);                // write coefficients to file

  double resolution = 0.5;                      // resolution, larger the finer
  
  /****************Loop to find the maximum LDOS ***************************/
  VectorXd Rho_s(MAX_ITER+1);                     // array to store all elecric LDOS
  MatrixXd Coeffs(Nc,MAX_ITER+1);                 // matrix to store all coefficients
  /*********** k = 0 *******************************************************/
  double rho_s=0;                                 // electric LDOS at center computed with 3 scatter simulation
  double *dfdx = new double[Nc];                  // gradient of LDOS vs coeffs
  for (int i=0;i<Nc;++i)                          // gradient initialzed to 0
    dfdx[i] = 0.0;
  /* temporary parameters used in optimization *****************************/
  VectorXd x_k(Nc);                               // x_k
  VectorXd dx_k(Nc);                              // delta x_k (p_k)
  VectorXd dfdx_k(Nc);                            // df/dx(x_k) 
  MatrixXd H_k = MatrixXd::Constant(Nc,Nc,0.0);   // Hessian matrix
  MatrixXd I = MatrixXd::Constant(Nc,Nc,0.0);     // Identity matrix
  for (int i=0;i<Nc;++i){
    I(i,i) = 1.0;
    H_k(i,i) = 1.0;                               // initialzied with 1  
  }
  H_k(0,0) = 1./coeffs[0];                        // initialzied with 1/coeff    
  VectorXd y_k(Nc);
  double f_temp = 0;                              // temporary f_(k+1)
  VectorXd dfdx_temp = VectorXd::Constant(Nc,0.0);// temporary df/dx(k+1)
  VectorXd x_temp = VectorXd::Constant(Nc,0.0);   // temporary x_k used in line search

  /**************compute the gradient of LDOS at current coeffs **********/
  LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,resolution);
  Rho_s(0) = rho_s;
  std::cout << "rho_s: " << rho_s << std::endl;

  for (int i=0;i<Nc;++i){
    Coeffs(i,0) = coeffs[i];
    x_k(i) = coeffs[i];
    dfdx_k(i) = dfdx[i];
  }
    
  /* parameters that used in line search for alpha *************************/
  double c_1 = 1e-4;
  double c_2 = 1.5 ;                             
  double tau = 0.5;
  double alpha = 1.0;                           // initial value for alpha
  /* loop for line search ***************************************************/
  std::cout << "Start line searching..." << std::endl;
  for (int p=0;p<MAX_ALPHA;++p){
    std::cout << "p: " << p << std::endl << std::endl;

    dx_k = alpha*H_k*dfdx_k;
    while(dx_k.norm()>sqrt(Nc)){
      alpha *= tau;
      dx_k = alpha*H_k*dfdx_k;
    }
    x_temp = x_k+dx_k;
    LDOS_gradient(lambda_0,x_temp,Nc,f_temp,dfdx_temp, Chi,d_min,resolution);
    std::cout <<f_temp<<" " << rho_s << std::endl;//-c_1*alpha*(dx_k.transpose()*dfdx_k).sum()<< std::endl;
    //std::cout <<abs((dx_k.transpose()*dfdx_temp).sum())<<" " << c_2*abs((dx_k.transpose()*dfdx_k).sum()) <<endl;
    /*
    if (f_temp>=(rho_s-c_1*alpha*(dx_k.transpose()*dfdx_k).sum())
      &&(abs((dx_k.transpose()*dfdx_temp).sum())<=c_2*abs((dx_k.transpose()*dfdx_k).sum()))
      &&(rho_s-c_1*alpha*(dx_k.transpose()*dfdx_k).sum()>0)){
      std::cout << "Sucessful line search!" << std::endl;
      break;
      }
    //*/
    if (f_temp>=rho_s){
      std::cout << "Sucessful line search!" << std::endl;
      break;
      }
    else
      alpha *= tau;
    if (p==(MAX_ALPHA-1)){
      std::cout << " Line searching failed to converge..." << std::endl;
      std::cout << " alpha: " << alpha << std::endl;
      std::cout << "Rho_s: " << Rho_s << std::endl;
      std::cout << "Coeffs: " << std::endl << Coeffs << std::endl;
      exit (EXIT_FAILURE);
    }
  }
  for (int i=0;i<Nc;++i)
    coeffs[i] = x_temp(i);
  y_k = dfdx_temp-dfdx_k;

  /**************************************************************************/
  /*********** Main Loop ****************************************************/
  /**************************************************************************/
  for (int k=1;k<MAX_ITER;++k){
    std::cout << "k: " << k << std::endl;
    /* Initialize temporary parameters to 0 *********************************/
    rho_s=f_temp;                                        //f_k

    /**************compute the gradient of LDOS at current coeffs **********/
    //LDOS_gradient(lambda_0, coeffs, num_coeffs, rho_s, dfdx, Chi,d_min ,resolution);
    Rho_s(k) = f_temp;
    std::cout << "rho_s: " << f_temp << std::endl;
    /* construct Hessian matrix H_k ******************************************/
    for (int i=0;i<Nc;++i){
      Coeffs(i,k) = x_temp(i);
      x_k(i) = x_temp(i);
      dfdx_k(i) = dfdx_temp(i);
    }
    
    H_k = (I-dx_k*y_k.transpose()/(dx_k.transpose()*y_k).sum())*H_k
          * (I-y_k*dx_k.transpose()/(y_k.transpose()*dx_k).sum())
          + (dx_k*dx_k.transpose())/(y_k.transpose()*dx_k).sum();
    /* loop for line search ***************************************************/
    std::cout << "Start line searching..." << std::endl;
    alpha = 1.0;
    for (int p=0;p<MAX_ALPHA;++p){
      std::cout << std::endl << " k: " << k << "   p: " << p << std::endl;

      dx_k = alpha*H_k*dfdx_k;
      while(dx_k.norm()>sqrt(Nc)){
        alpha *= tau;
        dx_k = alpha*H_k*dfdx_k;
      }
      x_temp = x_k+dx_k;
      LDOS_gradient(lambda_0,x_temp,Nc,f_temp,dfdx_temp, Chi,d_min,resolution);   
      std::cout <<f_temp<<" " << rho_s << std::endl;//-c_1*alpha*(dx_k.transpose()*dfdx_k).sum()<< std::endl;
      //std::cout <<abs((dx_k.transpose()*dfdx_temp).sum())<<" " << c_2*abs((dx_k.transpose()*dfdx_k).sum()) <<endl;
      /*
      if (f_temp>=(rho_s-c_1*alpha*(dx_k.transpose()*dfdx_k).sum())
        &&(abs((dx_k.transpose()*dfdx_temp).sum())<=c_2*abs((dx_k.transpose()*dfdx_k).sum()))
        &&(rho_s-c_1*alpha*(dx_k.transpose()*dfdx_k).sum()>0)){
        std::cout << "Sucessful line search!" << std::endl;
        break;
        }
      //*/
      if (f_temp>=rho_s){
        std::cout << "Sucessful line search!" << std::endl;
        break;
        }
      else
        alpha *= tau;
      if (p==(MAX_ALPHA-1)){
        std::cout << " Line searching failed to converge..." << std::endl;
        std::cout << " alpha: " << alpha << std::endl;
        std::cout << "Rho_s: " << Rho_s << std::endl;
        std::cout << "Coeffs: " << std::endl << Coeffs<< std::endl;
        exit (EXIT_FAILURE);
      }
    }
    /*************************************************************************/
    for (int i=0;i<Nc;++i)
      coeffs[i] = x_temp(i);
    y_k = dfdx_temp-dfdx_k;
    
    if ((Rho_s(k)-Rho_s(k-1)/Rho_s(k))<1e-2)
      break;
  }

  Rho_s(MAX_ITER) = f_temp;
  for (int i=0;i<Nc;++i)
    Coeffs(i,MAX_ITER) = x_temp(i);
  std::cout << "Rho_s: " << Rho_s << std::endl;
  std::cout << "Coeffs: " << std::endl << Coeffs << std::endl;
  double rho_limit = pow(k_0*d_min,-3)*abs(Chi*Chi)/imag(Chi)*(1+pow(k_0*d_min,2));
  std::cout << "d_min: " << d_min << std::endl;
  std::cout << "rho_limit: " << rho_limit << std::endl;
  delete[] dfdx;
  delete[] coeffs;
}