#include <stdio.h>
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
#define Z_0 376.730313    // Free space impedance, unit : Ohms
/***************************************************************/
/***************************************************************/
/***************************************************************/

void WriteCoeff(double *coeffs, const unsigned int num_coeffs){
  ofstream coeff_file;                          //write to file
  coeff_file.open ("coeff.txt");
  for(int i=0;i<num_coeffs;++i){
    //coeffs[i] *= 2*sqrt(pi);
    coeff_file << coeffs[i] << " ";
  }
  coeff_file.close();
}


void LDOS_gradient(double lambda_0, double *coeffs, const unsigned int num_coeffs, 
                  double &rho_s, double *dfdx, cdouble &Chi, double &d_min, 
                  int min_mesh, int max_mesh, double resolution)
{
  /***************************************************************/
  /* Initialze Frequency (Wavelength) and write to file **********/
  /***************************************************************/
  //double lambda_0 = 300; //wavelength in unit of nm
  double omega_0 = 2*pi/lambda_0*1e3; //omega in unit of 3e14 rad/sec
  double rho_0 = omega_0*omega_0/2/pi/pi; //free space electric LDOS
  ofstream omega_file;  //write to file
  omega_file.open ("OmegaValues.dat");
  omega_file << omega_0;
  omega_file.close();

  /***************************************************************/
  /* Initialze Spherical Harmonic funciton and Create Mesh *******/
  /***************************************************************/
  WriteCoeff(coeffs,num_coeffs);
  double el = 0;
  for (int i=0;i<num_coeffs;++i)
    el += abs(coeffs[i]);
  el = el/100.0/resolution; //maximum edge length
  std::string MeshCommand = "./SHMesher coeff.txt "
      +std::to_string(num_coeffs)+" "+std::to_string(el)+" Ellipsoid.msh "
      +std::to_string(min_mesh)+" "+std::to_string(max_mesh);
  system(MeshCommand.c_str()); //create mesh
  std::cout << "Mesh complete." << std::endl;
  geometry::SphericalHarmonics *sh;  //create spherical harmonic class
  sh=new geometry::SphericalHarmonics(num_coeffs, coeffs, -1e6, 1e6);
  double *l=new double[num_coeffs];
  double *m=new double[num_coeffs];
  sh->get_lm(l,m);
  //std::cout <<"r(0,0): " << sh->radius(0.0,0.0) << std::endl;  
  //std::cout <<"Ylm(0,0,0,0): " << sh->Ylm(0,0,0.0,0.0) << std::endl;  
  /***************************************************************/
  /* Get Eps and Mu of the specified material  *******************/
  /***************************************************************/
  MatProp *MP=new MatProp("FILE_Silver.txt"); //outside silver
  if (MP->ErrMsg)
   ErrExit(MP->ErrMsg);
  cdouble Eps_1; 
  cdouble Mu_1;
  MP->GetEpsMu( cdouble(omega_0,0.0), &Eps_1, &Mu_1);
  Chi = Eps_1-cdouble(1.0,0.0);
  cdouble Eps_2(1.0,0.0); //inside air
  cdouble Mu_2(1.0,0.0); 
  //std::cout << "Eps: " << real(Eps_2) << "+" << imag(Eps_2) << "i" << std::endl;
  //std::cout << "Mu: " << real(Mu) << "+" << imag(Mu) << "i" << std::endl;
  /***************************************************************/
  /* Initialze Parameters ****************************************/
  /***************************************************************/  
  char *GeoFile=(char *)"SphericalHarmonics.scuffgeo";
  char *OmegaFile=(char *)"OmegaValues.dat";
  double psLoc[3]={0.0,0.0,0.0};             int npsLoc=1;
  int npsStrength=1;            
  char *EPFile=(char *)"EPFile";
  /***************************************************************/
  /* Start Main Computing ****************************************/
  /***************************************************************/  
  rho_s=0.0;  //electric LDOS at center computed with 3 scatter simulation
  //double *dfdx = new double[num_coeffs]; //gradient initialzed to 0
  for (int i=0;i<num_coeffs-1;++i)
    dfdx[i] = 0;
  /* loop through 3 electric dipole source ***********************/
  HMatrix **PSDMatrix = new HMatrix*[3];
  for (int ps=0;ps<3;++ps){
    //HMatrix *PSDMatrix[ps]=0;    //scatter matrix to store frequency,fields
    cdouble EH[6];           //scattered field at the center
    //initialize point source
    cdouble psStrength[3]={0.0,0.0,0.0};
    for (int i=0;i<3;++i)
      if (i==ps)
        psStrength[i] = 1.0;
    std::cout << "Start scattering processing for point source : " << ps << std::endl;
    PSDMatrix[ps]=scuff_scatter(GeoFile, OmegaFile, psLoc, npsLoc, psStrength, npsStrength,EPFile,EH);
    rho_s += imag(EH[ps]);
    // loop every panel and get the fields
    
    for (int np=0;np<PSDMatrix[ps]->NR;++np){
      cdouble dfdx_temp=0;
      dfdx_temp=((Eps_2 - Eps_1) * (PSDMatrix[ps]->GetEntry(np,9) * PSDMatrix[ps]->GetEntry(np,9)
                    + PSDMatrix[ps]->GetEntry(np,10) * PSDMatrix[ps]->GetEntry(np,10) 
                    + PSDMatrix[ps]->GetEntry(np,11) * PSDMatrix[ps]->GetEntry(np,11))
                    + (1. / Eps_1 - 1. / Eps_2) * Z_0 * Z_0
                    * PSDMatrix[ps]->GetEntry(np,4) * PSDMatrix[ps]->GetEntry(np,4)); 
      dfdx_temp *= PSDMatrix[ps]->GetEntry(np,3)/pi/omega_0;
      double theta=0;
      double phi=0;
      sh->theta_phi(real(PSDMatrix[ps]->GetEntry(np,0)),real(PSDMatrix[ps]->GetEntry(np,1)),real(PSDMatrix[ps]->GetEntry(np,2)),theta,phi);
      // store the derivative in every coefficient respectively
      double d_temp = sh->radius(theta,phi);
      for (int ic=0;ic<num_coeffs-1;++ic){
        double Ynx;
        Ynx = sh->Ylm(l[ic],m[ic],theta,phi);
        dfdx[ic] += 2*Ynx*imag(dfdx_temp)/rho_0*sqrt(d_temp-coeffs[0])/sqrt(1+sh->drdt(theta,phi)*sh->drdt(theta,phi)+sh->drdp(theta,phi)*sh->drdp(theta,phi));
      }
      d_temp *= 1e3;
      if (d_min>d_temp)
        d_min=d_temp;
    }
    //*/
    delete PSDMatrix[ps];
  }
  rho_s = rho_s/pi/omega_0/rho_0+1;
  
  delete[] PSDMatrix;
  delete MP;
  delete[] l;
  delete[] m;
  delete sh;
}

void LDOS_gradient(double lambda_0, double r1, double r2, 
                  double &rho_s, double *dfdx, cdouble &Chi, double &d_min, 
                  int min_mesh, int max_mesh, double resolution)
{
  /***************************************************************/
  /* Initialze Frequency (Wavelength) and write to file **********/
  /***************************************************************/
  //double lambda_0 = 300; //wavelength in unit of nm
  double omega_0 = 2*pi/lambda_0*1e3; //omega in unit of 3e14 rad/sec
  double rho_0 = omega_0*omega_0/2/pi/pi; //free space electric LDOS

  /***************************************************************/
  /* Create Mesh *******/
  /***************************************************************/
  double el = sqrt(r1*r1+r2*r2)/50.0/resolution; //maximum edge length
  std::string MeshCommand = "./GeoMesher ellipsoid "
      +std::to_string(r1)+" "+std::to_string(r1)+" "+std::to_string(r2)+" "
      +std::to_string(el)+" 0 "+std::to_string(min_mesh)+" "+std::to_string(max_mesh);
  std::cout << MeshCommand.c_str() << std::endl;
  system(MeshCommand.c_str()); //create mesh
  std::cout << "Mesh complete." << std::endl;
  /***************************************************************/
  /* Get Eps and Mu of the specified material  *******************/
  /***************************************************************/
  MatProp *MP=new MatProp("FILE_Silver.txt"); //outside silver
  if (MP->ErrMsg)
   ErrExit(MP->ErrMsg);
  cdouble Eps_1; 
  cdouble Mu_1;
  MP->GetEpsMu( cdouble(omega_0,0.0), &Eps_1, &Mu_1);
  Chi = Eps_1-cdouble(1.0,0.0);
  cdouble Eps_2(1.0,0.0); //inside air
  cdouble Mu_2(1.0,0.0); 
  //std::cout << "Eps: " << real(Eps_2) << "+" << imag(Eps_2) << "i" << std::endl;
  //std::cout << "Mu: " << real(Mu) << "+" << imag(Mu) << "i" << std::endl;
  /***************************************************************/
  /* Initialze Parameters ****************************************/
  /***************************************************************/  
  char *GeoFile=(char *)"SphericalHarmonics.scuffgeo";
  char *OmegaFile=(char *)"OmegaValues.dat";
  double psLoc[3]={0.0,0.0,0.0};             int npsLoc=1;
  int npsStrength=1;            
  char *EPFile=(char *)"EPFile";
  /***************************************************************/
  /* Start Main Computing ****************************************/
  /***************************************************************/  
  rho_s=0.0;  //electric LDOS at center computed with 3 scatter simulation
  //double *dfdx = new double[num_coeffs]; //gradient initialzed to 0
  for (int i=0;i<2;++i)
    dfdx[i] = 0;
  /* loop through 3 electric dipole source ***********************/
  HMatrix **PSDMatrix = new HMatrix*[3];
  for (int ps=0;ps<3;++ps){
    //HMatrix *PSDMatrix[ps]=0;    //scatter matrix to store frequency,fields
    cdouble EH[6];           //scattered field at the center
    //initialize point source
    cdouble psStrength[3]={0.0,0.0,0.0};
    for (int i=0;i<3;++i)
      if (i==ps)
        psStrength[i] = 1.0;
    std::cout << "Start scattering processing for point source : " << ps << std::endl;
    PSDMatrix[ps]=scuff_scatter(GeoFile, OmegaFile, psLoc, npsLoc, psStrength, npsStrength,EPFile,EH);
    rho_s += imag(EH[ps]);
    // loop every panel and get the fields
    
    for (int np=0;np<PSDMatrix[ps]->NR;++np){
      cdouble dfdx_temp=0;
      dfdx_temp=((Eps_2 - Eps_1) * (PSDMatrix[ps]->GetEntry(np,9) * PSDMatrix[ps]->GetEntry(np,9)
                    + PSDMatrix[ps]->GetEntry(np,10) * PSDMatrix[ps]->GetEntry(np,10) 
                    + PSDMatrix[ps]->GetEntry(np,11) * PSDMatrix[ps]->GetEntry(np,11))
                    + (1. / Eps_1 - 1. / Eps_2) * Z_0 * Z_0
                    * PSDMatrix[ps]->GetEntry(np,4) * PSDMatrix[ps]->GetEntry(np,4)); 
      dfdx_temp *= PSDMatrix[ps]->GetEntry(np,3)/pi/omega_0;
      double theta=0;
      double phi=0;
      double xx=real(PSDMatrix[ps]->GetEntry(np,0));
      double yy=real(PSDMatrix[ps]->GetEntry(np,1));
      double zz=real(PSDMatrix[ps]->GetEntry(np,2));
      double rr=sqrt(xx*xx+yy*yy+zz*zz);
      // store the derivative in every coefficient respectively
      double dr1 = rr/r1/r1/r1*(xx*xx+yy*yy);
      double dr2 = rr/r2/r2/r2*zz*zz;
      double normal = 1/rr/sqrt((xx*xx+yy*yy)/r1/r1/r1/r1+zz*zz/r2/r2/r2/r2);
      dfdx[0] += dr1*imag(dfdx_temp)/rho_0*normal;
      dfdx[1] += dr2*imag(dfdx_temp)/rho_0*normal;
      rr *= 1e3;
      if (d_min>rr)
        d_min=rr;
    }
    //*/
    delete PSDMatrix[ps];
  }
  rho_s = rho_s/pi/omega_0/rho_0+1;
  
  delete[] PSDMatrix;
  delete MP;
}