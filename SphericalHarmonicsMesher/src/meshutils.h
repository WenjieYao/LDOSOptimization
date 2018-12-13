
#ifndef MESHUTILS_H
#define MESHUTILS_H

#include <string>
//#include <vector>

int GmshNumTri(const char *gmsh_file);
// Convert OFF file format to GMSH file format
void ConvertOFFtoGMSH(const char *OFFFileIn, const char *GMSHFileOut, int verbose=0);

class Mesh
{
    public:
        /*** Variables ***/
        const char* file;
        std::string type; // file = filename containing mesh, type = comsol/gmsh
        double volume, surfArea;
        int numNodes, numTriangles;
        double **nodeData; // (x,y,z,nx,ny,nz) for each node
        double **triangles;  // (node1,node2,node3,nx,ny,nz,cx,cy,cz,area) where ci = i^th centroid, area = panel area for each triangle

        /*** Methods ***/
        Mesh(const char* filename, std::string meshType);
        Mesh(const Mesh& copyMesh);
        ~Mesh();
        Mesh CreateCoating(const char* filename, double t); 
        void CalcVolSA();
        void SetType(std::string meshType); 
        void ProcessFile(double **triangles, double **nodeData);
        void ProcessGmshFile(double **triangles, double **nodeData);
        void ProcessComsolFile(double **triangles, double **nodeData);
        void GetSizesFromFile();
        void GetSizesFromGmshFile();
        void GetSizesFromComsolFile();
        void WriteFile(); // Write to file in either Comsol or Gmsh style
        void WriteGmshFile();
        void WriteComsolFile();
        void WriteNodeNormals(const char *);
        //std::vector<int> GetVAR();
};

#endif
