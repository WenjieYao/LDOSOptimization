
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/constants/constants.hpp>

namespace geometry {

void genRandomSpherePoints(unsigned numPoints, double *theta, double *phi) {
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> ZeroOne(0.0,1.0);

    double t_min = 0, p_min = 0;
    double t_max = boost::math::double_constants::pi;
    double p_max = boost::math::double_constants::two_pi;
    for(unsigned i=0; i<numPoints; ++i) {
        theta[i] = acos((cos(t_min) - cos(t_max)) * ZeroOne(rng) + cos(t_max) );
        phi[i] = (p_max - p_min) * ZeroOne(rng) + p_min;
    }
}

} // namespace geometry
