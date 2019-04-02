
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "cgalDistMesh.h"
#include "basisfunctions.h"
#include "meshutils.h"

double p1, p2, p3, p4; // global because of required signature for dist_function
std::string filename;
double *coeffs;
unsigned num_coeffs;
geometry::BaseGeo *bg;

void read_coeffs_from_file(std::string filename, double *c) {
    std::ifstream infile(filename.c_str());
    if (infile.is_open()) {
        unsigned i = 0;
        while (infile >> c[i++]) {};
        infile.close();
    }
}

double min(double a, double b) {
    return (a<b) ? a : b;
}


double spherical_harmonics_dist_function(double x, double y, double z) {
    return bg->level_set_function(x, y, z);
}

// args = "type", right now "box" or "ellipsoid"
// 		then 3 parameters (i.e. rx, ry, rz for ellipsoid, p1, p2, p3 for box)
// 		then maxEdgeLength
int main(int argc, char *argv[]) {
    const char *offFile = "out.off";

    // TODO: add rounding for the cone
    if(argc < 5) {
        std::cerr << "SphericalHarmonics:  runMesher coeff_file num_coeffs el output_file [min_tri max_tri] \n";
        return -1;
    }

    double bound;
    filename = std::string(argv[1]);
    num_coeffs = atoi(argv[2]);
    double  d_bound = strtod(argv[7], 0);
    std::cout << "d_bound: " << d_bound << std::endl;
    coeffs = new double[num_coeffs];
    read_coeffs_from_file(filename, coeffs);
    /*
    std::cout << "coefficients: ";
    for( int pi=0; pi<num_coeffs ; ++pi){
        std::cout << coeffs[pi] << " ";
    }
    std::cout << std::endl;
    //*/
    const double maxEdgeLength = strtod(argv[3], 0);
    //const double maxSquaredRadius = 0.01;
    std::string output_file = std::string(argv[4]);
    const char *gmshFile = output_file.c_str();

    int min_tri = (argc==8) ? atoi(argv[5]) : 50;
    int max_tri = (argc==8) ? atoi(argv[6]) : 20000;
    mesh_inputs mi = {min_tri, max_tri, 0, maxEdgeLength, true, offFile, gmshFile};
    mesh_outputs mo;

    bg = new geometry::SphericalHarmonics(num_coeffs, coeffs, -1e6, 1e6);
    bound = bg->bound_squared_radius();
    std::cout << "bound: " << bound << std::endl;
    cgalDistMesh(spherical_harmonics_dist_function, bound, &mi, &mo, d_bound);
    ConvertOFFtoGMSH(offFile, gmshFile, 1);
    system("rm out.off");
    return 0;
}
