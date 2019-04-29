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
#include "LDOSmath.h"

#define pi 3.1415926535897  
#define Z_0 376.730313    // Free space impedance, unit : Ohms
/***************************************************************/
/***************************************************************/
/***************************************************************/

void WriteCoeff(double *coeffs, const unsigned int num_coeffs){
  std::ofstream coeff_file;                          //write to file
  coeff_file.open ("coeff.txt");
  for(int i=0;i<num_coeffs;++i){
    //coeffs[i] *= 2*sqrt(pi);
    coeff_file << coeffs[i] << " ";
  }
  coeff_file.close();
}

bool CompareCoeff(double *coeffs, const unsigned int num_coeffs){
  std::ifstream file("coeff.txt");
  double *Excoeffs = new double[num_coeffs];
  if (file.is_open()) {
    for(int i=0;i<num_coeffs;++i){
      file >> Excoeffs[i];
    }
    file.close();
  }
  double diff = 0;
  for(int i=0;i<num_coeffs;++i)
    diff += abs(Excoeffs[i]-coeffs[i]);
  bool rv = false;
  if(diff>0.1)
    rv = true;
  delete[] Excoeffs;
  return rv;
}

void ReadMSH(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<int> &v1, std::vector<int> &v2, std::vector<int> &v3){
    std::ifstream file("Ellipsoid.msh");
    int ni=0;
    int nNodes=0;
    int nElements=0;
    if (file.is_open()) {
        std::string line;
        while(std::getline(file,line)){
            ni++;
            std::stringstream ss(line);
            if(ni<5)
                continue;
            else if(ni==5)
                ss >> nNodes;
            else if((ni<6+nNodes)&&(nNodes!=0)){
                int n_temp;
                double x_temp,y_temp,z_temp;
                ss >> n_temp >> x_temp >> y_temp >> z_temp;
                //std::cout << x_temp << " " << y_temp << " " << z_temp << std::endl;
                x.push_back(x_temp);
                y.push_back(y_temp);
                z.push_back(z_temp);
            }
            else if((ni==(8+nNodes))&&(nNodes!=0))
                ss >> nElements;
            else if((ni>(8+nNodes))&&(ni<=(8+nNodes+nElements))){
                int n_temp;
                int vt1,vt2,vt3;
                ss >> n_temp >> n_temp >> n_temp >> n_temp >> n_temp >> vt1 >> vt2 >> vt3;
                //std::cout << vt1 << " " << vt2 << " " << vt3 << std::endl;
                v1.push_back(vt1);
                v2.push_back(vt2);
                v3.push_back(vt3);
            }    
        }
    file.close();
  }
}

void WriteMSH(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<int> &v1, std::vector<int> &v2, std::vector<int> &v3){
    std::ofstream file("temp.msh");
    int nNodes=x.size();
    int nElements=v1.size();
    if (file.is_open()) {
        file << std::setprecision(16) << "$MeshFormat\n"; 
        file << "2.2 0 8\n";
        file << "$EndMeshFormat\n"; 
        file << "$Nodes\n";
        file << nNodes << "\n";
        for(int i=0;i<nNodes;++i)
            file << (i+1) << " " << x[i] << " " << y[i] << " " << z[i] << "\n";
        file << "$EndNodes\n";
        file << "$Elements\n";
        file << nElements << "\n";
        for(int i=0;i<nElements;++i)
            file << (i+1) << " 2 2 2 2 " << v1[i] << " " << v2[i] << " " << v3[i] << "\n";
        file << "$EndElements\n";
    file.close();
  }
  system("rm Ellipsoid.msh");
  system("mv temp.msh Ellipsoid.msh");
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
  std::ofstream omega_file;  //write to file
  omega_file.open ("OmegaValues.dat");
  omega_file << omega_0;
  omega_file.close();

  /***************************************************************/
  /* Initialze Spherical Harmonic funciton and Create Mesh *******/
  /***************************************************************/
  geometry::SphericalHarmonics *sh;  //create spherical harmonic class
  sh=new geometry::SphericalHarmonics(num_coeffs, coeffs, -1e6, 1e6);
  double *l=new double[num_coeffs];
  double *m=new double[num_coeffs];
  sh->get_lm(l,m);

  if(CompareCoeff(coeffs,num_coeffs)){
  WriteCoeff(coeffs,num_coeffs);
  double el = 0;
  for (int i=0;i<num_coeffs;++i)
    el += abs(coeffs[i]);
  el = el/100.0/resolution; //maximum edge length
  std::string MeshCommand = "./SHMesher coeff.txt "
      +std::to_string(num_coeffs)+" "+std::to_string(el)+" Ellipsoid.msh "
      +std::to_string(min_mesh)+" "+std::to_string(max_mesh)+" "+std::to_string(0.02+d_min/1e3);
  system(MeshCommand.c_str()); //create mesh
  }
  else{
    std::vector<double> x,y,z;
    std::vector<int> v1,v2,v3;
    ReadMSH(x,y,z,v1,v2,v3);
    for(int i=0;i<x.size();++i){
        double theta,phi;
        sh->theta_phi(x[i],y[i],z[i],theta,phi);
        double strech = sh->radius(theta,phi)/sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        x[i] *= strech;
        y[i] *= strech;
        z[i] *= strech;
    }
    WriteMSH(x,y,z,v1,v2,v3);
  }

  std::cout << "Mesh complete." << std::endl;
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
  d_min=100;
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
      double d_sum =0;
      double Ynx=0;
      for (int ic=0;ic<num_coeffs-1;++ic){
        Ynx = sh->Ylm(l[ic],m[ic],theta,phi);
        d_sum += coeffs[ic+1]*Ynx;
        }
      for (int ic=0;ic<num_coeffs-1;++ic){
        Ynx = sh->Ylm(l[ic],m[ic],theta,phi);
        dfdx[ic] += 2*Ynx*imag(dfdx_temp)/rho_0*d_sum/sqrt(1+pow(sh->drdt(theta,phi)/d_temp,2)+pow(sh->drdp(theta,phi)/d_temp/sin(theta),2));
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
  //std::cout << MeshCommand.c_str() << std::endl;
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