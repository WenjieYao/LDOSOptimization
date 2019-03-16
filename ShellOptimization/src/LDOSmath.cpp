#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <iostream>
#include <fstream>
#include <string>
#include "scuff-scatter.h"
#include "LDOSmath.h"

#define pi 3.1415926535897  
#define Z_0 376.730313    // Free space impedance, unit : Ohms
/***************************************************************/
/***************************************************************/
/***************************************************************/

void LDOS_eval(double lambda_0, double r1, double r2, 
                  double &rho_s, int min_mesh, int max_mesh, double resolution)
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
  /* Create Mesh *******/
  /***************************************************************/
  double el = r1/50.0/resolution; //maximum edge length
  std::string MeshCommand1 = "./GeoMesher ellipsoid "
      +std::to_string(r1)+" "+std::to_string(r1)+" "+std::to_string(r1)+" "
      +std::to_string(el)+" 0 "+std::to_string(min_mesh)+" "+std::to_string(max_mesh);
  //std::cout << MeshCommand.c_str() << std::endl;
  system(MeshCommand1.c_str()); //create mesh
  system("mv Ellipsoid.msh Ellipsoid1.msh");

  el = r2/50.0/resolution; //maximum edge length
  std::string MeshCommand2 = "./GeoMesher ellipsoid "
      +std::to_string(r2)+" "+std::to_string(r2)+" "+std::to_string(r2)+" "
      +std::to_string(el)+" 0 "+std::to_string(min_mesh)+" "+std::to_string(max_mesh);
  //std::cout << MeshCommand.c_str() << std::endl;
  system(MeshCommand2.c_str()); //create mesh
  system("mv Ellipsoid.msh Ellipsoid2.msh");
  std::cout << "Mesh complete." << std::endl;
  /***************************************************************/
  /* Initialze Parameters ****************************************/
  /***************************************************************/  
  char *GeoFile=(char *)"Shell.scuffgeo";
  char *OmegaFile=(char *)"OmegaValues.dat";
  double psLoc[3]={0.0,0.0,0.0};             int npsLoc=1;
  int npsStrength=1;            
  char *EPFile=(char *)"EPFile";
  /***************************************************************/
  /* Start Main Computing ****************************************/
  /***************************************************************/  
  rho_s=0.0;  //electric LDOS at center computed with 3 scatter simulation
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
    delete PSDMatrix[ps];
  }
  rho_s = rho_s/pi/omega_0/rho_0+1;
  
  delete[] PSDMatrix;
}