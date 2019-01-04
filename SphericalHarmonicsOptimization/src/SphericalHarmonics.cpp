
//#include <fstream>
#include <iostream>
#include "BFInclude/SphericalHarmonics.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/constants/constants.hpp>

namespace geometry {

void SphericalHarmonics::init_lm() {
    unsigned i = 0;
    l[i] = 0;
    m[i] = 0;
    for (i = 1; i < N; ++i) {
        if (m[i-1] != l[i-1]) {
            m[i] = m[i-1] + 1;
            l[i] = l[i-1];
        } else {
            l[i] = l[i-1] + 1;
            m[i] = -l[i];
        }
    }
}

SphericalHarmonics::~SphericalHarmonics() {
    delete[] t_rand;
    delete[] p_rand;
    delete[] l;
    delete[] m;
}

// note: may want to convert this to a vector matrix multiply
double SphericalHarmonics::radius(double t, double p) const {
    double r = 0;
    double temp_r = 0;
    for (unsigned i=1; i<N; ++i)
        temp_r += c[i] * Ylm(l[i-1], m[i-1], t, p);
    r = c[0] + temp_r*temp_r;
    return r;
}

double SphericalHarmonics::level_set_function(double x, double y, double z) const {
    return sqrt(x*x + y*y + z*z) - radius(x,y,z);
}

double SphericalHarmonics::bound_squared_radius() const {
    double bound = 0;
    for (unsigned i=0; i<num_constr_points; ++i) {
        double r = radius(t_rand[i], p_rand[i]);
        if (2.*r*r > bound)
            bound = 2. * r * r;
    }
    return bound;
}

double SphericalHarmonics::radius(double x, double y, double z) const {
    double t, p;
    theta_phi(x,y,z,t,p);
    return radius(t,p);
}

double SphericalHarmonics::drdc(unsigned i, double t, double p) const {
    return Ylm(l[i], m[i], t, p);
}

double SphericalHarmonics::drdt(double t, double p) const {
    double f = 0;
    for (unsigned i=1; i<N; ++i)
        f += c[i] * dYdTheta(l[i-1], m[i-1], t, p);
    f *= sqrt(radius(t,p)-c[0])*2;
    return f;
}

double SphericalHarmonics::drdp(double t, double p) const {
    double f = 0;
    for (unsigned i=1; i<N; ++i)
        f += c[i] * dYdPhi(l[i-1], m[i-1], t, p);
    f *= sqrt(radius(t,p)-c[0])*2;
    return f;
}

double SphericalHarmonics::grad_phi_mag(double t, double p) const {
    double r = radius(t, p);
    double dphidr = 1;
    double dphidt = -1. * drdt(t,p) / r;
    double dphidp = -1. * drdp(t,p) / (r*sin(t) + 1e-16);
    return sqrt(dphidr * dphidr + dphidt * dphidt + dphidp * dphidp);
}

//double *SphericalHarmonics::get_params(double *r) const {
    //for (unsigned i=0; i<N; ++i)
        //r[i] = c[i];
    //return r + N;
//}

void SphericalHarmonics::get_lm(double *li, double *mi) const {
    for (unsigned i=0; i<N; ++i) {
        li[i] = l[i];
        mi[i] = m[i];
    }
}

//double *SphericalHarmonics::get_bounds(double *minmax) const {
    //*(minmax++) = c_min;
    //*(minmax++) = c_max;
    //return minmax;
//}

//const double *SphericalHarmonics::update_params(const double *r) {
    //for (unsigned i=0; i<N; ++i)
        //c[i] = r[i];
    //return r + N;
//}

// dx_i / dc_j
double *SphericalHarmonics::surf_derivs_params(const double *xyz, unsigned n_xyz, double *dxdp) const {
    double theta=0, phi=0;
    for (unsigned i=0; i<n_xyz; ++i) {
        theta_phi(xyz[3*i], xyz[3*i+1], xyz[3*i+2], theta, phi);
        double gpm = grad_phi_mag(theta, phi);
        for (unsigned j=0; j<N; ++j) {
            *(dxdp++) = Ylm(l[j], m[j], theta, phi) / gpm;
        }
    }
    return dxdp;
}

void SphericalHarmonics::add_constraint(const char *constraint_type, double val) {
    if (!strcmp(constraint_type, "MIN_RADIUS")) {
        set_min_radius(val);
    } else if (!strcmp(constraint_type, "MAX_RADIUS")) {
        set_max_radius(val);
    } else {
        std::cerr << "Unknown constraint type added: " << constraint_type << std::endl;
        exit(EXIT_FAILURE);
    }
    //if (t_rand == NULL) {
        //t_rand = new double[num_constr_points];
        //p_rand = new double[num_constr_points];
        //geometry::genRandomSpherePoints(num_constr_points, t_rand, p_rand);
    //}
}

double *SphericalHarmonics::constraint_vals(unsigned c_ind, double *c_iter) const {
    unsigned r_index = c_ind % num_constr_points;
    if (c_ind < num_constr_points && r_min_active)
        c_iter = min_rad_constraint_vals(r_index, c_iter);
    else
        c_iter = max_rad_constraint_vals(r_index, c_iter);
    return c_iter;
}

double *SphericalHarmonics::min_rad_constraint_vals(unsigned r_ind, double *c_iter) const {
    *c_iter = r_min - radius(t_rand[r_ind], p_rand[r_ind]);
    return ++c_iter;
}

double *SphericalHarmonics::max_rad_constraint_vals(unsigned r_ind, double *c_iter) const {
    *c_iter = radius(t_rand[r_ind], p_rand[r_ind]) - r_max;
    return ++c_iter;
}

double *SphericalHarmonics::min_rad_constraint_derivs(unsigned r_ind, double *g_iter) const {
    for (unsigned i=0; i<N; ++i)
        *(g_iter++) = -1. * drdc(i, t_rand[r_ind], p_rand[r_ind]);
    return g_iter;
}

// TODO: should not recompute if already computed min_rad derivs
double *SphericalHarmonics::max_rad_constraint_derivs(unsigned r_ind, double *g_iter) const {
    for (unsigned i=0; i<N; ++i)
        *(g_iter++) = drdc(i, t_rand[r_ind], p_rand[r_ind]);
    return g_iter;
}

double *SphericalHarmonics::constraint_derivs(unsigned c_ind, double *g_iter) const {
    unsigned r_index = c_ind % num_constr_points;
    if (c_ind < num_constr_points && r_min_active)
        g_iter = min_rad_constraint_derivs(r_index, g_iter);
    else
        g_iter = max_rad_constraint_derivs(r_index, g_iter);
    return g_iter;
}

void SphericalHarmonics::set_min_radius(double val) {
    num_cfx += num_constr_points;
    r_min = val;
    r_min_active = true;
}

void SphericalHarmonics::set_max_radius(double val) {
    num_cfx += num_constr_points;
    r_max = val;
    r_max_active = true;
}

// loop over relevant portions of cval and grad functions, filling in values
void SphericalHarmonics::constraint_function(double *cval_iter, double *grad_iter) const {
    for (std::size_t i=0; i<num_cfx; ++i)
        cval_iter = constraint_vals(i, cval_iter);
    if (grad_iter != NULL)
        for (std::size_t i=0; i<num_cfx; ++i)
            grad_iter = constraint_derivs(i, grad_iter);
}

// Below: static methods to compute spherical harmonic basis functions
double SphericalHarmonics::Ylm(unsigned l, int m, double t, double p) {
    if (m<=0)
        return boost::math::spherical_harmonic_r(l, m, t, p);
    else
        return boost::math::spherical_harmonic_i(l, m, t, p);
}

double SphericalHarmonics::dYdPhi(unsigned l, int m, double t, double p) {
    if (m<=0)
        return -1. * m * boost::math::spherical_harmonic_i(l, m, t, p);
    else 
        return m * boost::math::spherical_harmonic_r(l, m, t, p);
}

double SphericalHarmonics::dYdTheta(unsigned l, int m, double t, double p) {
    double x = cos(t);
    double val = 0.;
    int li = (int) l;
    if (m < li)
        val += 0.5 * boost::math::legendre_p(l, m+1, x);
    if (m > -1 * li)
        val += -0.5 * (li + m) * (li - m + 1) * boost::math::legendre_p(l, m-1, x);

    val *= sqrt( (2.*li + 1.) / (4. * boost::math::double_constants::pi) 
            * boost::math::tgamma_ratio(li - m + 1, li + m + 1));
    if (m<=0)
        val *= cos(m * p);
    else
        val *= sin(m * p);
    return val;
}

void SphericalHarmonics::theta_phi(double x, double y, double z, double &t, double &p) {
    if (x==0 && y==0 && z==0) {
        t = 0;
        p = 0;
    } else {
        t = acos(z / sqrt(x*x + y*y + z*z));
        p = atan2(y, x);
    }
}

} // namespace Geometry
