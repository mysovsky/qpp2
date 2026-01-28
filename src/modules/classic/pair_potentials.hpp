#ifndef _QPP_PAIR_POTENTIALS_H
#define  _QPP_PAIR_POTENTIALS_H

#include <Eigen/Dense>
#include <mathf/lace3d.hpp>
#include <data/types.hpp>
#include <geom/ngbr.hpp>
#include <consts.hpp>
#include <classic/potentials.hpp>

namespace qpp{
  

// ---------------------------------------------
      
  template <class REAL>
  class buckingham_potential : public pair_potential<REAL>{
  public:
        
    REAL A, rho, C;
    
    using pair_potential<REAL>::atom1;
    using pair_potential<REAL>::atom2;
    using pair_potential<REAL>::rmax;

    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const
    {
      REAL r3 = r*r*r;
      REAL Aexp =  A*std::exp(-r/rho);
      REAL Cr6 = C/(r3*r3);
      E =  Aexp - Cr6;
      if (do_d1e)
	D1E = -Aexp/rho + 6e0*Cr6/r;
      if (do_d2e)
	D2E = Aexp/(rho*rho) - 42e0*Cr6/(r*r);
    }

    buckingham_potential(const STRING_EX & _at1, const STRING_EX & _at2, 
			 REAL _A, REAL _rho, REAL _C, REAL _rmax)
    {
      atom1 = _at1;
      atom2 = _at2;
      rmax = _rmax;
      
      A = _A;
      C = _C;
      rho = _rho;

      /*
      A = /hartree_to_ev;
      rho /= bohr_to_angs;
      
      REAL bohr6 = bohr_to_angs;
      bohr6 *= bohr6;
      bohr6 *= bohr6*bohr6;
      
      C /= hartree_to_ev * bohr6;
      
      rmax /= bohr_to_angs;
      */
    }
    
    buckingham_potential(ISTREAM &is){ is >> *this; }
    
    virtual ISTREAM& load(ISTREAM &is)
    {
      ISTREAM *ret = &(is >> atom1 >> atom2 >> A >> rho >> C  >> rmax);
      
      // [A] = eV, [rho] = A, [C] = eV*A^6;
      //
      // get atomic units:

      /*
      A /= hartree_to_ev;
      rho /= bohr_to_angs;
      
      REAL bohr6 = bohr_to_angs;
      bohr6 *= bohr6;
      bohr6 *= bohr6*bohr6;
      
      C /= hartree_to_ev * bohr6;
      
      rmax /= bohr_to_angs;
      */
      
      return *ret;
    }
      
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef buckingham_potential<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<buckingham_potential<REAL>, pair_potential<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, REAL, REAL, REAL, REAL>())
	.def_readonly("A", &SELF::A)
	.def_readonly("rho", &SELF::rho)
	.def_readonly("C", &SELF::C)
	;
    }

#endif
    
  };

  // ---------------------------------------------------------------
  
  template <class REAL>
  class buckingham4_potential : public pair_potential<REAL>{
  public:
        
    REAL A, rho, C;
    REAL ra, rm, rb;

  private:
    
    REAL a2,b2,c2,d2,f2;
    REAL a3,b3,c3,d3;

    void init_poly(){
      
      REAL rb4 = rb*rb;
      rb4 *= rb4;

      REAL den3 = rb4*(rb-rm);
      den3 *= den3;

      a3 = -2*C*(8*rb - 7*rm)/den3;
      b3 = 3*C*(9*rb*rb-7*rm*rm)/den3;
      c3 = -6*C*rm*rb*(9*rb-8*rm)/den3;
      d3 = -2*C*rb*rb*(6*rb*rb - 21*rb*rm + 14*rm*rm)/den3;

      REAL u1a = A*exp(-ra/rho);
      REAL u1ap = -u1a/rho;
      REAL u1app = -u1ap/rho;

      REAL u3m = ((a3*rm + b3)*rm + c3)*rm+d3;
      REAL u3mpp = 6*a3*rm + 2*b3;

      REAL rma = rm-ra;
      REAL rma2 = rma*rma;
      REAL rma3 = rma2*rma;
      REAL rma4 = rma3*rma;
      REAL rma5 = rma4*rma;

      a2 = -0.5*((u1app-u3mpp)*rma2 + 6*u1ap*rma + 12*(u1a - u3m))/rma5;
      b2 = -((u1app-1.5*u3mpp)*rma2 + 7*u1ap*rma + 15*(u1a - u3m))/rma4;
      c2 = 0.5*u3mpp/rma - (u1a-u3m)/rma3 - a2*rma2 + b2*rma;
      d2 = 0.5*u3mpp;
      f2 = u3m;
    }
    
  public:
    
    using pair_potential<REAL>::atom1;
    using pair_potential<REAL>::atom2;
    using pair_potential<REAL>::rmax;

    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const
    {
      if (r < ra){
	REAL Aexp =  A*std::exp(-r/rho);
	E = Aexp;
	if (do_d1e)
	  D1E = -Aexp/rho;
	if (do_d2e)
	  D2E = Aexp/(rho*rho);
      }
      else if (r < rm){
	REAL rrm = r - rm;
	E = (((a2*rrm + b2)*rrm + c2)*rrm + d2)*rrm*rrm + f2;
	if (do_d1e)
	  D1E = (((5*a2*rrm + 4*b2)*rrm + 3*c2)*rrm + 2*d2)*rrm;
	if (do_d2e)
	  D2E = ((20*a2*rrm + 12*b2)*rrm + 6*c2)*rrm + 2*d2;
      }
      else if (r < rb){
	E = ((a3*r + b3)*r + c3)*r + d3;
	if (do_d1e)
	  D1E = (3*a3*r + 2*b3)*r + c3;
	if (do_d2e)
	  D2E = 6*a3*r + 2*b3;
      }
      else {
	REAL r3 = r*r*r;
	REAL Cr6 = C/(r3*r3);
	E =  -Cr6;
	if (do_d1e)
	  D1E = 6e0*Cr6/r;
	if (do_d2e)
	  D2E = - 42e0*Cr6/(r*r);
      }
    }

    buckingham4_potential(const STRING_EX & _at1, const STRING_EX & _at2, 
			  REAL _A, REAL _rho, REAL _C, REAL _ra, REAL _rm, REAL _rb, REAL _rmax)
    {
      atom1 = _at1;
      atom2 = _at2;
      rmax = _rmax;
      
      A = _A;
      C = _C;
      rho = _rho;
      ra = _ra;
      rm = _rm;
      rb = _rb;

      init_poly();
      
      /*
      A = /hartree_to_ev;
      rho /= bohr_to_angs;
      
      REAL bohr6 = bohr_to_angs;
      bohr6 *= bohr6;
      bohr6 *= bohr6*bohr6;
      
      C /= hartree_to_ev * bohr6;
      
      rmax /= bohr_to_angs;
      */
    }
    
    buckingham4_potential(ISTREAM &is){ is >> *this; }
    
    virtual ISTREAM& load(ISTREAM &is)
    {
      ISTREAM *ret = &(is >> atom1 >> atom2 >> A >> rho >> C  >> ra >> rm >> rb >> rmax);
      init_poly();
      return *ret;
    }
      
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef buckingham4_potential<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<buckingham4_potential<REAL>, pair_potential<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, REAL, REAL, REAL, REAL, REAL, REAL, REAL>())
	.def_readonly("A", &SELF::A)
	.def_readonly("rho", &SELF::rho)
	.def_readonly("C", &SELF::C)
	.def_readonly("ra", &SELF::ra)
	.def_readonly("rm", &SELF::rm)
	.def_readonly("rb", &SELF::rb)
	;
    }

#endif
    
  };

  // ---------------------------------------------
      
  template <class REAL>
  class morse_potential : public pair_potential<REAL>{
  public:


    REAL D,a,r0;
    
    using pair_potential<REAL>::atom1;
    using pair_potential<REAL>::atom2;
    using pair_potential<REAL>::rmax;

    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const
    {
      REAL ex =  std::exp(-a*(r-r0));
      REAL t = 1e0 -ex;
      E = D*t*t - D;
      if (do_d1e)
	D1E = 2e0*a*D*t*ex;
      if (do_d2e)
	D2E = 2e0*a*a*D*ex*(2e0 - ex);
    }

    morse_potential(const STRING_EX & _at1, const STRING_EX & _at2, 
		    REAL _D, REAL _a, REAL _r0, REAL _rmax)
    {
      atom1 = _at1;
      atom2 = _at2;
      rmax = _rmax;

      D = _D;
      a = _a;
      r0 = _r0;
      /*
	D /= hatree_to_ev;
	a *= bohr_to_angs;
	r0 /= bohr_to_angs;
	
	rmax /= bohr_to_angs;
      */
    }
    
    morse_potential(ISTREAM &is){ is >> *this; }

    virtual ISTREAM& load(ISTREAM &is)
    {
      ISTREAM *ret = &(is >> atom1 >> atom2 >> D >> a >> r0 >> rmax);
      
      // [D] = eV, [a] = A^-1, [r0] = A;
      //
      // get atomic units:
      
      /*
	D /= hatree_to_ev;
	a *= bohr_to_angs;
	r0 /= bohr_to_angs;
	
	rmax /= bohr_to_angs;
      */
      
      return *ret;
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef morse_potential<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<morse_potential<REAL>, pair_potential<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, REAL, REAL, REAL, REAL>())
	.def_readonly("D", &SELF::D)
	.def_readonly("a", &SELF::a)
	.def_readonly("r0", &SELF::r0)
	;
    }

#endif
    
  };

  // ---------------------------------------------
      
  template <class REAL>
  class cutcoulomb_potential : public pair_potential<REAL>{
  public:


    REAL U,a;
    
    using pair_potential<REAL>::atom1;
    using pair_potential<REAL>::atom2;
    using pair_potential<REAL>::rmax;

    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const
    {
      REAL den2 = 1e0 + r*r/(a*a), den =  std::sqrt(den2);
      E = U/den;
      if (do_d1e)
	D1E = -U*r/(a*a*den*den2);
      if (do_d2e)
	D2E = U*(2*r*r/(a*a)-1e0)/(den2*den2*den);
    }

    cutcoulomb_potential(const STRING_EX & _at1, const STRING_EX & _at2, 
		    REAL _U, REAL _a, REAL _rmax)
    {
      atom1 = _at1;
      atom2 = _at2;
      rmax = _rmax;

      U = _U;
      a = _a;
    }
    
    cutcoulomb_potential(ISTREAM &is){ is >> *this; }

    virtual ISTREAM& load(ISTREAM &is)
    {
      ISTREAM *ret = &(is >> atom1 >> atom2 >> U >> a  >> rmax);
      
      // [D] = eV, [a] = A^-1, [r0] = A;
      //
      // get atomic units:
      
      /*
	D /= hatree_to_ev;
	a *= bohr_to_angs;
	r0 /= bohr_to_angs;
	
	rmax /= bohr_to_angs;
      */
      
      return *ret;
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef cutcoulomb_potential<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<cutcoulomb_potential<REAL>, pair_potential<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, REAL, REAL, REAL>())
	.def_readonly("U", &SELF::U)
	.def_readonly("a", &SELF::a)
	;
    }

#endif
    
  };

  // ---------------------------------------------
      
  template <class REAL>
  class spring_potential : public pair_potential<REAL>{
  public:

    REAL k;
    
    using pair_potential<REAL>::atom1;
    using pair_potential<REAL>::atom2;
    using pair_potential<REAL>::rmax;
    using pair_potential<REAL>::too_close;

    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const
    {
      E = 0.5*k*r*r;
      if (do_d1e)
	D1E = k*r;
      if (do_d2e)
	D2E = k;
    }

    spring_potential(const STRING_EX & _at1, const STRING_EX & _at2, 
		     REAL _k, REAL _rmax)
    {
      atom1 = _at1;
      atom2 = _at2;
      rmax = _rmax;

      k = _k;
      /*
	D /= hatree_to_ev;
	a *= bohr_to_angs;
	r0 /= bohr_to_angs;
	
	rmax /= bohr_to_angs;
      */
    }
    
    spring_potential(ISTREAM &is){ is >> *this; }

    virtual ISTREAM& load(ISTREAM &is)
    {
      ISTREAM *ret = &(is >> atom1 >> atom2 >> k  >> rmax);
      
      // [D] = eV, [a] = A^-1, [r0] = A;
      //
      // get atomic units:
      
      /*
	D /= hatree_to_ev;
	a *= bohr_to_angs;
	r0 /= bohr_to_angs;
	
	rmax /= bohr_to_angs;
      */
      
      return *ret;
    }


    virtual void energy_and_derivs(REAL & E,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E, bool do_d1e,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2D, bool do_d2e,
				   geometry<REAL> & geom,
				   const neighbours_table<REAL> & ngbr,
				   const std::vector<int> & active_atoms,
				   const std::vector<int> & frozen_atoms,
				   const std::vector<int> & core_shells) const
    {
      std::vector<unsigned char> atomregs(geom.nat(),0);
      // 1 - in active region
      // 2 - in frozen region
      // 0 - neither
      if (active_atoms.size() == 0)
	// empty active_atoms means all atoms belong to it
	for (auto i = atomregs.begin();  i!=atomregs.end(); i++) *i = 1;
      else
	{
	  for (int i:active_atoms)
	    atomregs[i] = 1;
	  for (int i:frozen_atoms)
	    atomregs[i] = 2;
	}
      if (core_shells.size()!=geom.nat())
	throw std::runtime_error("spring potential requested by core-shell pairs are not defined\n");
      for (int i=0; i<geom.nat(); i++)
	if (atomregs[i] > 0 && geom.atom(i) == atom1){
	  int j = core_shells[i];
	  //if (j>i)
	  //  continue;
	  if (atomregs[j] == 0)
	    continue;
	  if (atomregs[i] == 2 && atomregs[j] == 2)
	    continue;
	  index J;
	  bool found = false;
	  for (int k = 0; k<ngbr.n(i); k++)
	    if (ngbr(i,k)(0) == j){
	      J = ngbr(i,k);
	      found = true;
	      break;
	    }
	  if (!found){
	    std::cout << "Warning: core(shell) #"<< i << " lost its counterpart\n";
	    continue;
	  }	    
	  vector3<REAL> rvec = geom.pos(J) - geom.pos(i);
	  REAL r = rvec.norm();
	  REAL e, de, d2e;
	  e12(r,e,de,d2e,do_d1e,do_d2e);
	  E += e;
	  if (do_d1e && r > too_close)
	    {
	      vector3<REAL> devec = de*rvec/r;
	      D1E.template block<1,3>(i,0) -= devec;
	      D1E.template block<1,3>(j,0) += devec;
	      if (do_d2e)
		{
		  // fixme - implement this
		}
	    }
	}
      
      /*
      for (int i=0; i<geom.nat(); i++)
	if (atomregs[i] > 0 && geom.atom(i) == atom1)
	  for (int j=0; j<ngbr.n(i); j++)
	    {
	      index k = ngbr(i,j);
	      if (atomregs[k] == 0 || geom.atom(k)!=atom2)
		continue;
	      // fixme - what about periodic cells?
	      if (atomregs[i] == 2 && atomregs[k] ==2 )
		continue;
	      if (atom1==atom2 && i>=k)
		continue;
	    }
      */
      std::cout << "E= "  << E << "\n";
    }


#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    typedef spring_potential<REAL> SELF;
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<spring_potential<REAL>, pair_potential<REAL> >(m, pyname)
	.def(py::init<const STRING_EX &, const STRING_EX &, REAL, REAL>())
	.def_readonly("k", &SELF::k)
	;
    }

#endif
    
  };  

};

#endif
