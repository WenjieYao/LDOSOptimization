
#include <cmath>
#include <iostream>
#include "BFInclude/LevelSetFunctions.h"

namespace geometry {

namespace ls_function {

double min(double a, double b) {
    return (a<b) ? a : b;
}
double max(double a, double b) {
    return (a>b) ? a : b;
}

// set operations for level-set functions (i.e. sign functions)
double intersection(double phi1, double phi2) {
    return max(phi1, phi2);
}
double union_two(double phi1, double phi2) {
    return min(phi1, phi2);
}

// vector operations
double vec_dot_product(const double a[2], const double b[2]) {
    return a[0]*b[0] + a[1]*b[1];
}
void vec_subtract(const double a[2], const double b[2], double c[2]) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
}
double vec_length(const double x[2]) {
    return sqrt(vec_dot_product(x,x));
}
void vec_normalize(double x[2]) {
    double xlength = vec_length(x);
    x[0] /= xlength;
    x[1] /= xlength;
}
double tri_det(const double v0[2], const double v1[2]) {
    return v0[0] * v1[1] - v0[1] * v1[0];
}

// basic shapes
double circle(double r, double x, double y) {
    return sqrt(x*x + y*y) - r;
}

// source mathworld.wolfram.com/TriangleInterior.html
double triangle(const double x0[2], const double x1[2], const double x2[2], double x, double y) {
    double v1[2] = {x1[0]-x0[0], x1[1]-x0[1]};
    double v2[2] = {x2[0]-x0[0], x2[1]-x0[1]};
    double v[2] = {x,y};
    double a = (tri_det(v,v2) - tri_det(x0,v2)) / tri_det(v1,v2);
    double b = -(tri_det(v,v1) - tri_det(x0,v1)) / tri_det(v1,v2);
    if ((a>0) && (b>0) && (a+b<1))
        return -a * b * vec_length(v1) * vec_length(v2);
    else
        return fabs(a) * fabs(b);
}

double rectangle(const double x0[2], const double v0[2], const double v1[2], double x, double y) {
    double xy[2] = {x - x0[0] - 0.5 * v0[0] - 0.5 * v1[0],
                    y - x0[1] - 0.5 * v0[1] - 0.5 * v1[1]};
    // abs --> first quadrant in v0, v1 basis
    double a = fabs(vec_dot_product(xy, v0) / vec_dot_product(v0, v0));
    double b = fabs(vec_dot_product(xy, v1) / vec_dot_product(v1, v1));
    double dist;
    if ((a<0.5) && (b<0.5)) // inside
        dist = -1. * min((0.5 - a) * vec_length(v0), (0.5 - b) * vec_length(v1));
    else {
        double a0 = min(a,0.5);
        double b0 = min(b,0.5);
        dist = sqrt((a-a0)*(a-a0)*vec_dot_product(v0,v0) + (b-b0)*(b-b0)*vec_dot_product(v1,v1));
    }
    return dist;
}
double rectangle(double Lx, double Ly, double x, double y) {
    double x0[2] = {-0.5*Lx, -0.5*Ly};
    double v0[2] = {Lx,0};
    double v1[2] = {0,Ly};
    return rectangle(x0, v0, v1, x, y);
}

double rounded_rect(const double x0[2], const double v0[2], const double v1[2], double r_rnd, double x, double y) {
    double phi = rectangle(x0, v0, v1, x, y);
    double xy[2] = {x,y};
    double v0l = vec_length(v0);
    double v1l = vec_length(v1);
    double pos_cyl[2] = 
        {x0[0] + v0[0] * (1. - r_rnd / v0l) + v1[0] * (1 - r_rnd / v1l),
         x0[1] + v0[1] * (1. - r_rnd / v0l) + v1[1] * (1 - r_rnd / v1l)};
    phi = cylinder_rnd(phi, pos_cyl, v0, v1, r_rnd, xy);
    return phi;
}

double rounded_rect(double Lx, double Ly, double r_rnd, double x, double y) {
    double x0[2] = {-0.5*Lx, -0.5*Ly};
    double v0[2] = {Lx,0};
    double v1[2] = {0,Ly};
    return rounded_rect(x0, v0, v1, r_rnd, x, y);
}

void tri_circ_center(const double x0[2], const double x1[2], const double x2[2], double radius, double xc[2], double v1p[2], double v2p[2]) {
    double v1[2], v2[2];
    vec_subtract(x1, x0, v1);
    vec_subtract(x2, x0, v2);
    vec_normalize(v1);
    vec_normalize(v2);
    double v1v2 = vec_dot_product(v1, v2);
    double a = radius / sqrt(1. - v1v2 * v1v2);
    xc[0] = x0[0] + a * (v1[0] + v2[0]);
    xc[1] = x0[1] + a * (v1[1] + v2[1]);
    v1p[0] = a * (v1v2 * v2[0] - v1[0]);
    v1p[1] = a * (v1v2 * v2[1] - v1[1]);
    v2p[0] = a * (v1v2 * v1[0] - v2[0]);
    v2p[1] = a * (v1v2 * v1[1] - v2[1]);
}

// would be nice to not redo computations of circle points every time
//   once the triangles, rectangles, etc. are actual BaseGeos, this won't be a problem.
double rounded_tri(const double x0[2], const double x1[2], const double x2[2], double r_rnd, double x, double y) {
    double xy[2] = {x, y};
    double phi = triangle(x0, x1, x2, x, y);
    double xc[2], v1p[2], v2p[2];
    tri_circ_center(x0, x1, x2, r_rnd, xc, v1p, v2p);
    phi = cylinder_rnd(phi, xc, v1p, v2p, r_rnd, xy);
    tri_circ_center(x1, x0, x2, r_rnd, xc, v1p, v2p);
    phi = cylinder_rnd(phi, xc, v1p, v2p, r_rnd, xy);
    tri_circ_center(x2, x0, x1, r_rnd, xc, v1p, v2p);
    phi = cylinder_rnd(phi, xc, v1p, v2p, r_rnd, xy);
    return phi;
}

// note that xy only need to be the position coordinates in the basis v0_cyl, v1_cyl
double cylinder_rnd(double phi, const double x0[2], const double v0_cyl[2], const double 
        v1_cyl[2], double r_rnd, const double xy[2]) {
    double pos_vec[2] = {xy[0] - x0[0], xy[1] - x0[1]};
    if ((tri_det(pos_vec, v0_cyl)/tri_det(v1_cyl, v0_cyl) > 0) 
     && (tri_det(pos_vec, v1_cyl)/tri_det(v0_cyl, v1_cyl) > 0)) {
        double phi_cyl = circle(r_rnd, pos_vec[0], pos_vec[1]);
        phi = intersection(phi, phi_cyl);
    }
    return phi;
}

} // namespace ls_function
} // namespace geometry

// cylinder centered at xyz0 with radius r and translation symmetry along cyl_axis (infinitely extruded)
//double circ_cylinder_rot_axis(double r, axis_type cyl_axis, double x, double y, double z) {
    //if (cyl_axis==ZAXIS)
        //return circle(r, x, y);
    //else if (cyl_axis==XAXIS)
        //return circle(r, y, z);
    //else if (cyl_axis==YAXIS)
        //return circle(r, x, z);
    //else
        //return -1;
//}
