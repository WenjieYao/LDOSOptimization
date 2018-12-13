
#include "basisfunctions.h"
//#include <iostream>

namespace geometry {

//BaseGeo::~BaseGeo() {
    //delete[] c;
//}

double *BaseGeo::get_params(double *r) const {
    for (unsigned i=0; i<N; ++i)
        r[i] = c[i];
    return r + N;
};

// This applies to all surfaces
double *BaseGeo::get_bounds(double *minmax) const {
    *(minmax++) = c_min;
    *(minmax++) = c_max;
    return minmax;
}

const double *BaseGeo::update_params(const double *r) {
    for (unsigned i=0; i<N; ++i)
        c[i] = r[i];
    return r + N;
}

bool BaseGeo::within_constraints() const {
    if (num_cfx>0) {
        double *cvals = new double[num_cfx];
        constraint_function(cvals, NULL);
        for (unsigned i=0; i<num_cfx; ++i) {
            if (cvals[i]>0) {
                delete[] cvals;
                return false;
            }
        }
        delete[] cvals;
    }
    return true;
}

} // namespace geometry
