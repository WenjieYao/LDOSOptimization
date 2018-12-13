
// Should I call this base class an AdaptiveStructure?
//
// Many of the functions below return pointers to doubles, which are generally
// one-past the end of the input vector, which has been iterated through.  This
// is similar to C++ iterators, and helps when building geometry classes out of
// multiple other geometry classes (i.e. CoatedParticle, potentially composite
// particles...)
//
// pure virtual functions must be implemented (above surf_deriv_params)
// bound_squared_radius() and level_set_function() may only be required for CGAL 
//   (this text goes above those two functions)
//
//  Add an optional BaseGeo constructor that takes a file instead of a double array!
#ifndef BASE_GEO_H
#define BASE_GEO_H

#include <iostream>
#include <cstdlib>

namespace geometry {

class BaseGeo {
    public:
        BaseGeo(int Nin, const double *coeff, double coeff_min, double coeff_max) : N(Nin), 
                c(new double[Nin]), c_min(coeff_min), c_max(coeff_max), num_cfx(0) {
            for (unsigned i=0; i<N; ++i)
                c[i] = coeff[i];
        }
        virtual ~BaseGeo() {delete[] c;};

        // default implementations, can be overriden (except get_number_of_params)
        unsigned get_number_of_params() const {return N;};
        virtual double *get_params(double *r) const;
        virtual double *get_bounds(double *minmax) const;
        virtual const double *update_params(const double *r);

        // return surface gradient as a function of geometrical param
        virtual double *surf_derivs_params(const double *xyz, unsigned n_xyz, double *dxdp) 
            const =0;

        // provide the boundary through a level set function, with max. radius
        virtual double bound_squared_radius() const =0;
        virtual double level_set_function(double x, double y, double z) const =0;

        // default implementation, can be overriden (note that by default, no constraints)
        virtual double *constraint_vals(unsigned c_index, double *c_iter) const {return c_iter;};
        virtual double *constraint_derivs(unsigned c_index, double *g_iter) const {return g_iter;};
        unsigned get_number_of_constraint_fxs() const {return num_cfx;};
        virtual void constraint_function(double *constr_vals, double *grad) const {};
        virtual void add_constraint(const char *type, double val) {
            std::cerr << "No constraints can be added to this geometry."
                      << "\n\tAttempted constraint: " << type << std::endl;
            exit(EXIT_FAILURE);
        }
        bool within_constraints() const;
        //void write_scuff_file() {  write_preamble(); this->write_to_scuff_file(); close_file(); };

    protected:
        unsigned N;
        double *c;
        double c_min;
        double c_max;
        
        unsigned num_cfx;
        //Permittivity epsilon;

    private:
        //void write_preamble() {}; // to write...
            //void close_file() {};
    };

} // namespace geometry

#endif
