
// all of the functions listed here return the level set function as defined by the function name and inputs

#ifndef LEVELSET_GEO_H 
#define LEVELSET_GEO_H 

namespace geometry {

namespace ls_function {

typedef enum axis_type {XAXIS, YAXIS, ZAXIS} axis_type;

double intersection(double phi1, double phi2);
double union_two(double phi1, double phi2);

double circle(double r, double x, double y);
double rectangle(double Lx, double Ly, double x, double y);
double rectangle(const double x0[2], const double v0[2], const double v1[2], double x, double y);
double triangle(const double x0[2], const double x1[2], const double x2[2], double x, double y);

double rounded_rect(double Lx, double Ly, double r_rnd, double x, double y);
double rounded_rect(const double x0[2], const double v0[2], const double v1[2], double r_rnd, double x, double y);
double rounded_tri(const double x0[2], const double x1[2], const double x2[2], double r_rnd, double x, double y);
//double circ_cylinder_rot_axis(double r, axis_type cyl_axis, double x, double y, double z);

// may need to do sphere rounding also

double cylinder_rnd(double phi, const double pos_cyl[2], const double v0_cyl[2], const double 
        v1_cyl[2], double r_rnd, const double xy[2]);
//double cylinder_rnd(double phi, const double pos_cyl[3], const double v0_cyl[3], 
        //const double v1_cyl[3], double r_rnd, axis_type cyl_axis, const double xyz[3]);
} // namespace ls_function
} // namespace geometry

#endif
