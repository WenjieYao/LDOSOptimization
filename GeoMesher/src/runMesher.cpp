
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "cgalDistMesh.h"
#include "meshutils.h"
#include "basisfunctions.h"

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

// for meshing an ellipsoid
double ellipsoid_dist_function(double x, double y, double z) {
    return x*x/p1/p1 + y*y/p2/p2 + z*z/p3/p3 - 1.;
}

double shell_dist_function(double x, double y, double z) {
    double r1 = p1;
    double r2 = p2;
    double theta_max = p3;
    double r = sqrt(x*x + y*y + z*z);
    double theta = acos(z/(r+1e-8));
    double phi = 0;
    if (theta < theta_max)
        phi = -1. * min(r-r1, r2-r);
    else
        phi = theta - theta_max;
    return phi;
}

double cone_dist_function(double x, double y, double z) {
    double r = p1;
    double h = p2;
    x = sqrt(x*x + y*y); // rotate into x-z plane
    double dist;
    if (z<=0 && z>=-h) {
        double x_edge = 1. / (r/h + h/r) * (r/h * x - z);
        double z_edge = 1. / (r/h + h/r) * (h/r * z - x);
        double dist_edge = sqrt((x-x_edge)*(x-x_edge) + (z-z_edge)*(z-z_edge));
        if (x < -z*r/h) {
            dist = -1. * min(z+h, dist_edge);
        } else {
            dist = dist_edge;
        }
    } else if (z>0)
        dist = z;
    else
        dist = -h - z;
    return dist;
}

typedef struct {
    double pt[3];// pt on plane
    double n[3]; // normal vector
} plane;

//double tetr_dist_function(double xi, double yi, double zi) {
    //double el = p1; // edge length
    //double v1 = el / sqrt(3.);
    //double v2 = -sqrt(3.) * el / 6.;
    //double v3 = el / 2.;
    //double v4 = sqrt(6.) * a / 3.;
    //plane p1 = { {v1,0,0}, {
//}

// p1 = larger radius, p2 = smaller radius
double ring_dist_function (double xi, double yi, double zi) {
    return geometry::ls_function::circle(p2, sqrt(xi*xi + yi*yi) - p1, zi);
}

double circ_cyl_dist_function (double xi, double yi, double zi) {
    zi = fabs(zi);
    double rho = sqrt(xi*xi + yi*yi);
    double r = p1;
    double Lz = p2;
    double r_rnd = p4;
    //std::cout << "r Lz r_rnd: " << r << " " << Lz << " " << r_rnd << std::endl;
    double phi = geometry::ls_function::circle(r, xi, yi);
    double phiz = geometry::ls_function::rounded_rect(2.*r, Lz, r_rnd, rho, zi);
    phi = geometry::ls_function::intersection(phi, phiz);
    return phi;
}

double box_dist_function(double x, double y, double z) {
    x = fabs(x); y = fabs(y); z = fabs(z); // first quadrant
    double Lx = p1;
    double Ly = p2;
    double Lz = p3;
    double r_rnd = p4;
    double xyz[3] = {x, y, z};

    double phi_xy = geometry::ls_function::rounded_rect(Lx, Ly, r_rnd, x, y);
    double phi_xz = geometry::ls_function::rounded_rect(Lx, Lz, r_rnd, x, z);
    double phi_yz = geometry::ls_function::rounded_rect(Ly, Lz, r_rnd, y, z);
    double phi = geometry::ls_function::intersection(phi_xy, phi_xz);
    phi = geometry::ls_function::intersection(phi, phi_yz);
    return phi;
}


double tri_cyl_dist_function(double x, double y, double z) {
    double b = p1, h = p2, t = p3, r_rnd = p4;

    x = x + h/2.; // shift by h/2
    y = fabs(y);
    z = fabs(z);
    
    if (z > t/2.)
        return z - t/2.;
    if (x < 0)
        return -x;
    if (y > b/2. - b*x/(2.*h))
        return y - b/2. + b*x/(2.*h);
    
    double xl = (4.*h*h*x - 2.*h*b*y + h*b*b) / (4.*h*h + b*b);
    double dist = (xl - x) * sqrt(1. + 4.*h*h/(b*b));
    dist = -1. * min(dist, min(x, t/2.-z));

    // round the two symmetric vertices
    xl = r_rnd * (1. + 1. / sqrt(1. + 4.*h*h/(b*b)));
    double yl = b/2. - b*xl/(2.*h);
    double y0 = yl + 2.*h*(r_rnd - xl)/b;
    if ((y>y0) && (y>yl+2.*h*(x-xl)/b))
        dist = sqrt((x-r_rnd)*(x-r_rnd) + (y-y0)*(y-y0)) - r_rnd;

    // round the third vertex
    yl = r_rnd / sqrt(1. + b*b/(4.*h*h));
    xl = h * (1. - 2.*yl/b);
    double x0 = xl - r_rnd / sqrt(1. + 4.*h*h/(b*b));
    if ((x>x0+b*(y-yl)/(2.*h)) && (x>x0-b*(y+yl)/(2.*h)))
        dist = sqrt((x-x0)*(x-x0) + y*y) - r_rnd;

    // round along the x=0, z=+-h/2 edges
    double z0 = t/2. - r_rnd;
    if ((x<r_rnd) && (z>t/2.-r_rnd))
        dist = sqrt((x-r_rnd)*(x-r_rnd) + (z-z0)*(z-z0)) - r_rnd;

    // round along other edge (actually two, symmetric)
    xl = (4.*h*h*x - 2.*b*h*y + b*b*h) / (b*b + 4.*h*h);
    yl = b/2. * (1. - xl/h);
    x0 = xl - r_rnd / sqrt(1. + 4*h*h/(b*b));
    y0 = yl + 2. * h * (x0 - xl) / b;
    if ((y>y0-b*(x-x0)/(2.*h)) && (z>t/2.-r_rnd))
        dist = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) - r_rnd;

    // round the intersection of first vertex and z=+-h/2 edges
    xl = r_rnd * (1. + 1. / sqrt(1. + 4.*h*h/(b*b)));
    yl = b/2. - b*xl/(2.*h);
    y0 = yl + 2.*h*(r_rnd - xl)/b;
    if ((y>y0) && (y>yl+2.*h*(x-xl)/b) && (z>z0))
        dist = sqrt((x-r_rnd)*(x-r_rnd) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) - r_rnd;

    // round the intersection of the third vertex and z=+-h/2 edges
    yl = r_rnd / sqrt(1. + b*b/(4.*h*h));
    xl = h * (1. - 2.*yl/b);
    x0 = xl - r_rnd / sqrt(1. + 4.*h*h/(b*b));
    if ((x>x0+b*(y-yl)/(2.*h)) && (x>x0-b*(y+yl)/(2.*h)) && (z>z0))
        dist = sqrt((x-x0)*(x-x0) + y*y + (z-z0)*(z-z0)) - r_rnd;
    
    return dist;

}

// args = "type", right now "box" or "ellipsoid"
// 		then 3 parameters (i.e. rx, ry, rz for ellipsoid, p1, p2, p3 for box)
// 		then maxEdgeLength
int main(int argc, char *argv[]) {
    const char *offFile = "out.off";
    const char *gmshFile = "Ellipsoid.msh";

    // TODO: add rounding for the cone
    if(argc < 6) {
        std::cerr << "Usage: runMesher type p1 p2 p3 edge_length [p4] [min_tri max_tri] \n"
            << "    ellipsoid: runMesher ellipsoid rx ry rz el [p4 ignored] [min_tri max_tri] \n"
            << "    box:       runMesher box       lx ly lz el r_rnd [min_tri max_tri] \n"
            << "    tri. cyl:  runMesher tricyl    b  h  t  el r_rnd [min_tri max_tri] \n" 
            << "    circ. cyl: runMesher circcyl   r  h  [p3 ignored] el r_rnd [min_tri max_tri] \n" 
            << "    spherocyl: runMesher sphcyl    w  L  [p3 ignored] el [p5 ignored] [min_tri max_tri] \n" 
            << "    cone:  runMesher cone  r  h  [p3 ignored] el [p4 ignored] [min_tri max_tri] \n" 
            << "    ring:  runMesher ring  r0 r1 [p3 ignored] el [p4 ignored] [min_tri max_tri] \n"
            << "    shell: runMesher shell r0 r1 theta_max    el [p4 ignored] [min_tri max_tri] \n";
            return -1;
    }

    double bound;
    
    p1 = strtod(argv[2], 0);
    p2 = strtod(argv[3], 0);
    p3 = strtod(argv[4], 0);
    bound = 2*(p1*p1+p2*p2+p3*p3);
    std::cout << "bound: " << bound << std::endl;
    const double maxEdgeLength = strtod(argv[5], 0);
    //const double maxSquaredRadius = 0.01;

    int min_tri = (argc==9) ? atoi(argv[7]) : 50;
    int max_tri = (argc==9) ? atoi(argv[8]) : 20000;
    mesh_inputs mi = {min_tri, max_tri, 1, maxEdgeLength, true, offFile, gmshFile};
    mesh_outputs mo;

    if (strcmp(argv[1],"ellipsoid")==0) {
        cgalDistMesh(ellipsoid_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1],"box")==0) {
        p4 = (argc>6) ? strtod(argv[6], 0) : 0;
        cgalDistMesh(box_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1],"tricyl")==0) {
        p4 = (argc>6) ? strtod(argv[6], 0) : 0;
        cgalDistMesh(tri_cyl_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1],"circcyl")==0) {
        p4 = (argc>6) ? strtod(argv[6], 0) : 0;
        cgalDistMesh(circ_cyl_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1],"sphcyl")==0) {
        p1 = p1/2.; // radius = width/2
        // p2 fixed
        p4 = p1;    // r_rnd = radius
        cgalDistMesh(circ_cyl_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1], "cone")==0) {
        cgalDistMesh(cone_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1], "ring")==0) {
        cgalDistMesh(ring_dist_function, bound, &mi, &mo);
    } else if (strcmp(argv[1], "shell")==0) {
        bound = 2.*p2*p2;
        cgalDistMesh(shell_dist_function, bound, &mi, &mo);
    } else {
        std::cerr << "Unknown geometry: " << argv[1] << "\n";
        return -1;
    }
    ConvertOFFtoGMSH(offFile, gmshFile, 1);
    system("rm out.off");
}
