#include <mathf/specfunc.hpp>
#include <cmath>

namespace qpp{

  double logfact(int n) {
    double f = 0e0;
    for (int k = 1; k<=n; k++)
      f += std::log(1e0*k);
    return f;
  }
  
  double atanxy(double x, double y) {
    
    const double eps=1e-8;
    double phi;
    
    if ( std::abs(x) <= eps ) {
      phi = pi / 2;
      if ( y < 0 )
	phi += pi;
    }
    else {
      phi = std::atan(y/x);
      if ( x  < 0 )
	phi += pi;
    }
    return phi;
  }
};

