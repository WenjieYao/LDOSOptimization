
#ifndef SPHERICALHARMONICS_GEO_H
#define SPHERICALHARMONICS_GEO_H

#include "BaseGeo.h"

namespace geometry {

void genRandomSpherePoints(unsigned n, double *t, double *p);

class SphericalHarmonics : public BaseGeo {
    public:
        SphericalHarmonics(unsigned n, double *coeff, double coeff_min, double coeff_max) :
            BaseGeo(n, coeff, coeff_min, coeff_max), num_constr_points(500), 
            r_min(0), r_max(0), r_min_active(false), r_max_active(false), l(new double[n]), 
            m(new double[n]), t_rand(new double[500]), p_rand(new double[500]) { init_lm(); 
                genRandomSpherePoints(num_constr_points, t_rand, p_rand); };
        ~SphericalHarmonics();
        double *surf_derivs_params(const double *xyz, unsigned n_xyz, double *dxdp) const;

        // MIN_RADIUS or MAX_RADIUS
        void add_constraint(const char *constraint_type, double val);
        double *constraint_vals(unsigned c_index, double *c_iter) const;
        double *constraint_derivs(unsigned c_index, double *g_iter) const;
        void constraint_function(double *constr_vals, double *grad) const;

        double bound_squared_radius() const;
        double level_set_function(double x, double y, double z) const;
        
        double radius(double theta, double phi) const;
        double radius(double x, double y, double z) const;
        void get_lm(double *l, double *m) const;
        static double Ylm(unsigned l, int m, double t, double p);
        static double dYdPhi(unsigned l, int m, double t, double p);
        static double dYdTheta(unsigned l, int m, double t, double p);
        static void theta_phi(double x, double y, double z, double &t, double &p);

    private:
        void set_min_radius(double val);
        void set_max_radius(double val);
        void init_lm();
        double drdc(unsigned i, double theta, double phi) const;
        double drdt(double theta, double phi) const;
        double drdp(double theta, double phi) const;
        double grad_phi_mag(double theta, double phi) const;
        double *min_rad_constraint_vals(unsigned r_ind, double *c_iter) const;
        double *max_rad_constraint_vals(unsigned r_ind, double *c_iter) const;
        double *min_rad_constraint_derivs(unsigned r_ind, double *c_iter) const;
        double *max_rad_constraint_derivs(unsigned r_ind, double *c_iter) const;

        unsigned num_constr_points;
        double r_min, r_max;
        bool r_min_active, r_max_active;
        double *l, *m;
        double *t_rand;
        double *p_rand;
};

} // namespace geometry

#endif
