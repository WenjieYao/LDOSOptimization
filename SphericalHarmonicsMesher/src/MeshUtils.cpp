
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

int GmshNumTri(const char *gmsh_file) {
    std::ifstream stream;
    stream.open(gmsh_file);
    std::string line;
    std::string keyword = "$Elements";
    int num_tri;
    while (getline(stream, line))
        if (line.find(keyword) != std::string::npos)
            stream >> num_tri;
    return num_tri;

}

void ConvertOFFtoGMSH(const char *OFFFileIn, const char *GMSHFileOut, int verbose) {
    // read in OFF File
    std::ifstream infile;
    infile.open(OFFFileIn);
    if(infile.is_open()) {
        int numVert, numPanel, tmp, i, cnt, *panel;
        double *vert;
        std::string line;
        getline(infile, line); // throw away first line
        infile >> numVert >> numPanel >> tmp;
        if(verbose)
            std::cout << "numVert: " << numVert << "\nnumPanel: " << numPanel << std::endl;
        getline(infile, line);
        getline(infile, line);
        vert = new double[3*numVert];
        panel = new int[3*numPanel];
        for(i=0; i<3*numVert; ++i)
            infile >> vert[i];
        cnt = 0;
        for(i=0; i<4*numPanel; ++i) {
            if(i%4==0)
                infile >> tmp;
            else
                infile >> panel[cnt++];
        }
        infile.close();

        // check orientation
        int reorient;
        double *x0, *x1, *x2, xc[3], n[3];
        //if( (reorient==1) && verbose )
        //	std::cout << "Flipping normals to face outward" << std::endl;

        // write out GMSH file
        std::ofstream outfile;
        outfile.open(GMSHFileOut);
        if( outfile.is_open() ) {
            outfile << std::setprecision(16) << "$MeshFormat \n2.2	0	8\n"
                << "$EndMeshFormat \n$Nodes \n" << numVert << " \n";
            for(i=0; i<numVert; ++i)
                outfile << i+1 << " " << vert[3*i] << " " << vert[3*i+1] 
                    << " " << vert[3*i+2] << "\n";
            outfile << "$EndNodes\n$Elements\n" << numPanel << "\n";
            for(i=0; i<numPanel; ++i) {
                reorient = 0;
                x0 = vert + 3*panel[3*i];
                x1 = vert + 3*panel[3*i+1];
                x2 = vert + 3*panel[3*i+2];
                xc[0] = (x0[0] + x1[0] + x2[0])/3.;
                xc[1] = (x0[1] + x1[1] + x2[1])/3.;
                xc[2] = (x0[2] + x1[2] + x2[2])/3.;
                n[0] = (x1[1]-x0[1])*(x2[2]-x0[2]) 
                    - (x1[2]-x0[2])*(x2[1]-x0[1]);
                n[1] = (x1[2]-x0[2])*(x2[0]-x0[0]) 
                    - (x1[0]-x0[0])*(x2[2]-x0[2]);
                n[2] = (x1[0]-x0[0])*(x2[1]-x0[1]) 
                    - (x1[1]-x0[1])*(x2[0]-x0[0]);
                if( (n[0]*xc[0] + n[1]*xc[1] + n[2]*xc[2])<0 )
                    reorient = 1;

                outfile << i+1 << " 2 2 2 2 " << panel[3*i]+1 << " " 
                    << panel[3*i+1+reorient]+1 << " " << panel[3*i+2-reorient]+1 << "\n";
            }
            outfile << "$EndElements\n";
            outfile.close();
        } else {
            std::cout << "Unable to write to GMSH file\n";
        }

    } else {
        std::cout << "Unable to open OFF file\n";
    }

}

