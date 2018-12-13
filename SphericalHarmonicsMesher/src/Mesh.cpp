
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <vector>

#include "meshutils.h"

using std::vector;

/*** Constructor(filename,meshType) ***/
Mesh::Mesh(const char* filename, std::string meshType)
{
    file = filename;
    SetType(meshType);
    Mesh::GetSizesFromFile(); // numNodes, numTriangles
    std::cout << "got sizes from file: numnodes = " << numNodes << "  numtri: " << numTriangles << std::endl;
    // allocate memory for nodeData and triangles
    nodeData = new double*[numNodes];
    triangles = new double*[numTriangles];
    for(int tc=0; tc<numNodes; tc++) 
        nodeData[tc] = new double[6];		
    for(int tc=0; tc<numTriangles; tc++) 
        triangles[tc] = new double[10];		

    // read node and triangle data from file, into nodeData and triangles
    volume = 0.0;
    surfArea = 0.0;
    std::cout << "before proc" << std::endl;
    Mesh::ProcessFile(nodeData, triangles);
    std::cout << "after proc" << std::endl;
}

/*** Copy constructor ***/
Mesh::Mesh(const Mesh& copyMesh) 
{
    file = copyMesh.file;
    type = copyMesh.type;
    volume = copyMesh.volume;
    surfArea = copyMesh.surfArea;
    numNodes = copyMesh.numNodes;
    numTriangles = copyMesh.numTriangles;
    //copy values of pointers nodeData and triangles
    int tc, j;
    nodeData = new double*[numNodes];
    triangles = new double*[numTriangles];
    for(tc=0; tc<numNodes; tc++) {
        nodeData[tc] = new double[6];
        for(j=0; j<6; j++)
            nodeData[tc][j] = copyMesh.nodeData[tc][j];
    }			
    for(tc=0; tc<numTriangles; tc++) {
        triangles[tc] = new double[10];
        for(j=0; j<10; j++) 
            triangles[tc][j] = copyMesh.triangles[tc][j];
    }
}

/*** Destructor ***/
Mesh::~Mesh() {
    for(int i=0; i<numNodes; i++) {
        delete[] triangles[i];
        delete[] nodeData[i];
    }
    delete[] triangles;
    delete[] nodeData;
}

/*
// Get Vert Avg Radius
vector<int> Mesh::GetVAR()
{
vector<int> VAR(numNodes);
vector<int> numNbrs(0, numNodes);
int tri1, tri2, tri3;
double dist1, dist2;
for(int i=0; i<numTriangles; ++i) {
for(int j=0; j<3; ++j) {
tri1 = triangles[i][j];
tri2 = triangles[i][(j+1)%3];
tri3 = triangles[i][(j+2)%3];
dist1 = sqrt( pow( nodeData[tri1][0] - nodeData[tri2][0], 2)
+ pow( nodeData[tri1][1] - nodeData[tri2][1], 2)
+ pow( nodeData[tri1][2] - nodeData[tri2][2], 2) );
dist2 = sqrt( pow( nodeData[tri1][0] - nodeData[tri3][0], 2)
+ pow( nodeData[tri1][1] - nodeData[tri3][1], 2)
+ pow( nodeData[tri1][2] - nodeData[tri3][2], 2) );
VAR[tri1] += dist1 + dist2;
numNbrs[tri1] += 2;
}
}
for(int i=0; i<numNodes; ++i) {
VAR[i] /= numNbrs[i];
}
return VAR;
}
*/

// Create a coating (in a new mesh object, returned) of thickness t
Mesh Mesh::CreateCoating(const char* filename, double t) 
{
    // First copy this object (changing to the new filename)
    Mesh a = *this; 
    a.file = filename;

    // Then extrude the xyz
    int i;
    for(i=0; i<a.numNodes; i++) {
        a.nodeData[i][0] += t * a.nodeData[i][3]; // x = x + nx*t
        a.nodeData[i][1] += t * a.nodeData[i][4]; // y = y + ny*t
        a.nodeData[i][2] += t * a.nodeData[i][5]; // z = z + nz*t
    }

    // Recalculate triangle data (nx,ny,nz,Cx,Cy,Cz,area all change)
    int node0, node1, node2;
    double x0, x1, x2, y0, y1, y2, z0, z1, z2, Cx, Cy, Cz;
    double ux, uy, uz, vx, vy, vz, nx, ny, nz, normN;
    for(i=0; i<a.numTriangles; i++) {
        node0 = a.triangles[i][0];
        node1 = a.triangles[i][1];
        node2 = a.triangles[i][2];

        // Get data points for this triangle
        x0 = a.nodeData[node0][0]; x1 = a.nodeData[node1][0]; x2 = a.nodeData[node2][0];
        y0 = a.nodeData[node0][1]; y1 = a.nodeData[node1][1]; y2 = a.nodeData[node2][1];
        z0 = a.nodeData[node0][2]; z1 = a.nodeData[node1][2]; z2 = a.nodeData[node2][2];

        // Centroids
        Cx = (x0+x1+x2)/3.0;
        Cy = (y0+y1+y2)/3.0;
        Cz = (z0+z1+z2)/3.0;

        // Normals
        ux = x1-x0; uy = y1-y0; uz = z1-z0;
        vx = x2-x0; vy = y2-y0; vz = z2-z0;
        nx = uy*vz - uz*vy;
        ny = uz*vx - ux*vz;
        nz = ux*vy - uy*vx;
        normN = sqrt(nx*nx+ny*ny+nz*nz) + 1e-8;
        nx /= normN; ny /= normN; nz /= normN;

        // normal vector for this triangle
        a.triangles[i][3] = nx;
        a.triangles[i][4] = ny;
        a.triangles[i][5] = nz;
        // Centroids
        a.triangles[i][6] = Cx;
        a.triangles[i][7] = Cy;
        a.triangles[i][8] = Cz;
        // Triangle (panel) area
        a.triangles[i][9] = 0.5*normN;
    }

    // Calculate volume and surface area
    a.CalcVolSA();

    // Write to file
    a.WriteFile();

    // Return the new object
    return a;
}

void Mesh::CalcVolSA() 
{
    volume = 0.0;
    surfArea = 0.0;
    for(int i=0; i<numTriangles; i++) {
        surfArea += triangles[i][9];
        volume += triangles[i][9] * triangles[i][6] * triangles[i][3]; // area * Cx * nx
    }
}

/*** Set mesh type as "Comsol" or "gmsh" ***/
void Mesh::SetType(std::string meshType) 
{
    if( (meshType.compare("Comsol")==0) || (meshType.compare("Gmsh")==0) ) {
        type = meshType;
    } else {
        // Throw error
        std::cerr << "Mesh type must be either Comsol or Gmsh \n";
    }
}

/*** Process file, pass to correct method by type ***/
void Mesh::ProcessFile(double **nodeData, double **triangles)
{
    if( type.compare("Gmsh")==0 )
        Mesh::ProcessGmshFile(nodeData, triangles);
    if( type.compare("Comsol")==0 )
        Mesh::ProcessComsolFile(nodeData, triangles);
}

/*** Process Gmsh file ***/
void Mesh::ProcessGmshFile(double **nodeData, double **triangles)
{
    std::string line;
    std::ifstream infile(file);
    if( infile.is_open() ) {
        // Read through the header, get the number of nodes
        while(infile >> line) {
            if( line.compare("$Nodes")==0 ) {
                infile >> line;
                break;
            }
        }

        // populate nodeData 1st 3 columns (x,y,z)
        int nodeCnt = 0;
        while( infile >> nodeCnt >> nodeData[nodeCnt][0] 
                >> nodeData[nodeCnt][1] >> nodeData[nodeCnt][2] ) {
            nodeData[nodeCnt-1][3] = nodeData[nodeCnt-1][4] = nodeData[nodeCnt-1][5] = 0;
            if( nodeCnt>=numNodes )
                break;
        }

        // Get the number of elements
        while( infile >> line ) {
            if( !line.compare("$Elements") ) {
                infile >> line;
                break;
            }
        }

        // read in triangle data
        int triCnt = 0;
        int triType, numTag, tagNo, node0, node1, node2;
        double x0, x1, x2, y0, y1, y2, z0, z1, z2, ux, uy, uz, vx, vy, vz;
        double nx, ny, nz, Cx, Cy, Cz, normN;
        while(triCnt < numTriangles-1) {
            infile >> triCnt >> triType >> numTag;
            for(int tagcc=0; tagcc<numTag; tagcc++)
                infile >> tagNo;
            infile >> node0 >> node1 >> node2;

            triCnt--; node0--; node1--; node2--; // 0-based indexing
            triangles[triCnt][0] = node0;
            triangles[triCnt][1] = node1;
            triangles[triCnt][2] = node2;

            // Get data points for this triangle
            x0 = nodeData[node0][0]; x1 = nodeData[node1][0]; x2 = nodeData[node2][0];
            y0 = nodeData[node0][1]; y1 = nodeData[node1][1]; y2 = nodeData[node2][1];
            z0 = nodeData[node0][2]; z1 = nodeData[node1][2]; z2 = nodeData[node2][2];

            // Centroids
            Cx = (x0+x1+x2)/3.0;
            Cy = (y0+y1+y2)/3.0;
            Cz = (z0+z1+z2)/3.0;

            // Normals
            ux = x1-x0; uy = y1-y0; uz = z1-z0;
            vx = x2-x0; vy = y2-y0; vz = z2-z0;
            nx = uy*vz - uz*vy;
            ny = uz*vx - ux*vz;
            nz = ux*vy - uy*vx;
            normN = sqrt(nx*nx+ny*ny+nz*nz) + 1e-8;
            nx /= normN; ny /= normN; nz /= normN;

            // add normal vector to node0
            nodeData[node0][3] += nx;
            nodeData[node0][4] += ny;
            nodeData[node0][5] += nz;
            // add normal vector to node1
            nodeData[node1][3] += nx;
            nodeData[node1][4] += ny;
            nodeData[node1][5] += nz;
            // add normal vector to node2
            nodeData[node2][3] += nx;
            nodeData[node2][4] += ny;
            nodeData[node2][5] += nz;
            // normal vector for this triangle
            triangles[triCnt][3] = nx;
            triangles[triCnt][4] = ny;
            triangles[triCnt][5] = nz;
            // Centroids
            triangles[triCnt][6] = Cx;
            triangles[triCnt][7] = Cy;
            triangles[triCnt][8] = Cz;
            // Triangle (panel) area
            triangles[triCnt][9] = 0.5*normN;

            surfArea += 0.5*normN;
            volume += 0.5*normN * Cx * nx;
        }

        infile.close();

        // ensure normal vectors in node data are unit vectors
        for(int i=0; i<numNodes; i++) {
            normN = sqrt( nodeData[i][3]*nodeData[i][3] + nodeData[i][4]*nodeData[i][4]
                    + nodeData[i][5]*nodeData[i][5] );
            nodeData[i][3] /= normN;
            nodeData[i][4] /= normN;
            nodeData[i][5] /= normN;
        }

    } else { 
        std::cerr << "Error trying to open file \n";
    }
}

/*** Process Comsol file ***/
void Mesh::ProcessComsolFile(double **nodeData, double **triangles) 
{
    char line[256];
    std::ifstream infile(file);
    if( infile.is_open() ) {
        // Read through the header
        while( infile.getline(line,256) ) 
            if( !strcmp(line,"# Mesh point coordinates") )
                break;

        // populate nodeData 1st 3 columns (x,y,z)
        int nodeCnt = 0;
        while( infile >> nodeData[nodeCnt][0] >> nodeData[nodeCnt][1] 
                >> nodeData[nodeCnt][2] ) {
            nodeData[nodeCnt][3] = nodeData[nodeCnt][4] = nodeData[nodeCnt][5] = 0;
            nodeCnt++;
            if( nodeCnt>=numNodes )
                break;
        }

        // Get to the number of elements
        while( infile.getline(line,256) ) {
            if( !strcmp(line,"3 # number of nodes per element") ) {
                infile.getline(line,256);
                infile.getline(line,256);
                break;
            }
        }

        // read in triangle data
        int triCnt = 0;
        int node0, node1, node2;
        double x0, x1, x2, y0, y1, y2, z0, z1, z2, ux, uy, uz, vx, vy, vz;
        double nx, ny, nz, Cx, Cy, Cz, normN;
        while(triCnt < numTriangles-1) {
            infile >> node0 >> node1 >> node2;

            triangles[triCnt][0] = node0;
            triangles[triCnt][1] = node1;
            triangles[triCnt][2] = node2;

            // Get data points for this triangle
            x0 = nodeData[node0][0]; x1 = nodeData[node1][0]; x2 = nodeData[node2][0];
            y0 = nodeData[node0][1]; y1 = nodeData[node1][1]; y2 = nodeData[node2][1];
            z0 = nodeData[node0][2]; z1 = nodeData[node1][2]; z2 = nodeData[node2][2];

            // Centroids
            Cx = (x0+x1+x2)/3.0;
            Cy = (y0+y1+y2)/3.0;
            Cz = (z0+z1+z2)/3.0;

            // Normals
            ux = x1-x0; uy = y1-y0; uz = z1-z0;
            vx = x2-x0; vy = y2-y0; vz = z2-z0;
            nx = uy*vz - uz*vy;
            ny = uz*vx - ux*vz;
            nz = ux*vy - uy*vx;
            normN = sqrt(nx*nx+ny*ny+nz*nz) + 1e-8;
            nx /= normN; ny /= normN; nz /= normN;

            // add normal vector to node0
            nodeData[node0][3] += nx;
            nodeData[node0][4] += ny;
            nodeData[node0][5] += nz;
            // add normal vector to node1
            nodeData[node1][3] += nx;
            nodeData[node1][4] += ny;
            nodeData[node1][5] += nz;
            // add normal vector to node2
            nodeData[node2][3] += nx;
            nodeData[node2][4] += ny;
            nodeData[node2][5] += nz;
            // normal vector for this triangle
            triangles[triCnt][3] = nx;
            triangles[triCnt][4] = ny;
            triangles[triCnt][5] = nz;
            // Centroids
            triangles[triCnt][6] = Cx;
            triangles[triCnt][7] = Cy;
            triangles[triCnt][8] = Cz;
            // Triangle (panel) area
            triangles[triCnt][9] = 0.5*normN;

            surfArea += 0.5*normN;
            volume += 0.5*normN * Cx * nx;

            triCnt++; // 0-based indexing
            if(triCnt >= numTriangles)
                break;
        }

        infile.close();

        // ensure normal vectors in node data are unit vectors
        for(int i=0; i<numNodes; i++) {
            normN = sqrt( nodeData[i][3]*nodeData[i][3] + nodeData[i][4]*nodeData[i][4]
                    + nodeData[i][5]*nodeData[i][5] );
            nodeData[i][3] /= normN;
            nodeData[i][4] /= normN;
            nodeData[i][5] /= normN;
        }
    }
}

void Mesh::GetSizesFromFile()
{
    if( !type.compare("Gmsh") ) 
        Mesh::GetSizesFromGmshFile();
    if( !type.compare("Comsol") )
        Mesh::GetSizesFromComsolFile();
}

void Mesh::GetSizesFromGmshFile()
{
    std::string line;
    std::ifstream infile(file);
    if( infile.is_open() ) {
        // Read through the header, get the number of nodes
        while(infile >> line) {
            if( !line.compare("$Nodes") ) 
                infile >> numNodes;
            if( !line.compare("$Elements") ) {
                infile >> numTriangles;
                break;
            }
        }
        infile.close();
    }

}

void Mesh::GetSizesFromComsolFile()
{
    char line[256];
    std::ifstream infile(file);
    if( infile.is_open() ) {
        // Read through the header, get the number of nodes
        while( infile.getline(line,256) ) {
            if( strcmp(line,"3 # sdim")==0 )
                infile >> numNodes;
            if( strcmp(line,"3 # number of nodes per element")==0 ) {
                infile >> numTriangles;
                break;
            }
        }
        infile.close();
    }
}

void Mesh::WriteFile()
{
    if( !type.compare("Gmsh") ) 
        Mesh::WriteGmshFile();
    if( !type.compare("Comsol") )
        Mesh::WriteComsolFile();
}

void Mesh::WriteGmshFile() // write Gmsh mesh to filename
{
    std::ofstream outdata(file);
    if( outdata.is_open() ) {
        outdata << std::setprecision(16) << "$MeshFormat \n2.2	0	8 \n"
            << "$EndMeshFormat \n$Nodes \n" << numNodes << " \n";
        int i;
        for(i=0; i<numNodes; i++) {
            outdata << i+1 << " " << nodeData[i][0] << " " << nodeData[i][1]
                << " " << nodeData[i][2] << " \n";
        }
        outdata << "$EndNodes \n$Elements \n" << numTriangles << "\n";

        // Note the conversion back to 1-based indexing
        for(i=0; i<numTriangles; i++) {
            outdata << i+1 << " 2 2 2 2 " << triangles[i][0]+1 << " " 
                << triangles[i][1]+1 << " " << triangles[i][2]+1 << "\n";
        }
        outdata << "$EndElements \n";
        outdata.close();
    }
}

// Note that the Comsol file will NOT be readable by Comsol, but will be readable as a Comsol file in my codes
void Mesh::WriteComsolFile() // write Comsol mesh to filename
{
    std::ofstream outdata(file);
    if( outdata.is_open() ) {
        outdata << "# Created by Mesh.cpp \n";
        outdata << "# Note that this file is NOT readable by Comsol, and should "
            << "only be used within the Mesh() class for optimization purposes\n\n";
        outdata << "# -------- Object 0 -------- \n";
        outdata << "3 # sdim \n";
        outdata << numNodes << " # number of mesh points \n";
        outdata << "0 # lowest mesh point index \n\n# Mesh point coordinates\n";
        int i;
        for(i=0; i<numNodes; i++) {
            outdata << nodeData[i][0] << " " << nodeData[i][1]
                << " " << nodeData[i][2] << "\n";
        }
        outdata << "\n3 # number of nodes per element\n";
        outdata << numTriangles << " # number of elements\n# Elements\n";
        for(i=0; i<numNodes; i++) {
            outdata << triangles[i][0] << " " << triangles[i][1] << " " 
                << triangles[i][2] << "\n";
        }
        outdata.close();
    }
}

void Mesh::WriteNodeNormals(const char *filename)
{
    std::cout << "File: " << filename << std::endl;
    std::ofstream outdata(filename);
    if( outdata.is_open() ) {
        for(int i=0; i<numNodes; ++i) 
            outdata << nodeData[i][3] << " " << nodeData[i][4] << " " << nodeData[i][5] << std::endl;
    }
    outdata.close();
}
