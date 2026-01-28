#ifndef _QPP_POTENTIALS_3B_H
#define  _QPP_POTENTIALS_3B_H

#include <Eigen/Dense>
#include <mathf/lace3d.hpp>
#include <data/types.hpp>
#include <geom/ngbr.hpp>
#include <consts.hpp>
#include <classic/potentials.hpp>

namespace qpp{

  template <class REAL>
  class three_harmonic : public potential_3body<REAL>{

    REAL theta0, k, k3, k4;
    
  public:
    
    using potential_3body<REAL>::atom1;
    using potential_3body<REAL>::atom2;
    using potential_3body<REAL>::atom3;
    using potential_3body<REAL>::r12max;
    using potential_3body<REAL>::r13max;
      
    virtual void e123(REAL r12, REAL r13, REAL theta,
		      REAL &e, std::vector<REAL> & d1e, std::vector<REAL> & d2e,
		      bool do_d1e = false, bool do_d2e = false) const
    {
      REAL dt = theta - theta0;
      e =  dt*dt*(k/2 + k3*dt/6 +k4*dt*dt/24);
      if (do_d1e){
	d1e[0] = d1e[1] = REAL(0);
	d1e[2] = (k + k3*dt/2 + k4*dt*dt/6)*dt;
      }
      if (do_d2e)
	// fixme - implement this
	{}
    }

    three_harmonic(){}

    three_harmonic(const STRING_EX & _atom1, const STRING_EX & _atom2, const STRING_EX & _atom3,
		   REAL _k, REAL _theta0, REAL _k3, REAL _k4, REAL _r12max, REAL _r13max){
      atom1 = _atom1;
      atom2 = _atom2;
      atom3 = _atom3;
      r12max = _r12max;
      r13max = _r13max;
      // degrees to rad
      theta0 = _theta0*pi/180;
      k = _k;
      k3 = _k3;
      k4 = _k4;
      //std::cout << "three body harmonic potential constructed\n";
    }

    three_harmonic(ISTREAM & is){ is >> *this; }
    
    virtual ISTREAM & load(ISTREAM & is)
    {
      ISTREAM *ret = &( is >> atom1 >> atom2 >> atom3 >> k >> k3 >> k4 >> r12max >> r13max );
      return *ret;
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef three_harmonic<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<three_harmonic<REAL>, potential_3body<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, const STRING_EX &,
	     REAL, REAL, REAL, REAL, REAL, REAL>())
	.def_readonly("theta0", &SELF::theta0)
	.def_readonly("k", &SELF::k)
	.def_readonly("k3", &SELF::k3)
	.def_readonly("k4", &SELF::k4)
	;
    }

#endif
    
  };

}

#endif
