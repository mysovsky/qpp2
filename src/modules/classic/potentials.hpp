#ifndef _QPP_POTENTIALS_H
#define  _QPP_POTENTIALS_H

#include <typeinfo> 
//#include <Eigen/Dense>
#include <mathf/lace3d.hpp>
#include <data/types.hpp>
#include <geom/xgeom.hpp>
#include <geom/ngbr.hpp>
#include <consts.hpp>
//#include <pyqpp/py_indexed_property.hpp>

namespace qpp{
  
  template <class REAL>
  class classical_potential{
  public:

    static REAL too_close;

    typedef classical_potential<REAL> SELF;
    
    virtual void energy_and_derivs(REAL & E,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E, bool do_d1e,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2D, bool do_d2e,
				   geometry<REAL> & geom,
				   const neighbours_table<REAL> & ngbr,
				   const std::vector<int> & active_atoms,
				   const std::vector<int> & frozen_atoms,
				   const std::vector<int> & core_shells) const {}
    
    REAL energy(geometry<REAL> & geom,
		const neighbours_table<REAL> & ngbr,
		const std::vector<int> & active_atoms = {},
		const std::vector<int> & frozen_atoms = {},
		const std::vector<int> & core_shells = {} ) const{
      REAL E = 0e0;
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d1e;
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d2e;
      energy_and_derivs(E, empty_d1e, false, empty_d2e, false, geom, ngbr,
			active_atoms, frozen_atoms, core_shells);
      return E;
    }

    void d1e(REAL & E, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E,
	     geometry<REAL> & geom,
	     const neighbours_table<REAL> & ngbr,
	     const std::vector<int> & active_atoms = {},
	     const std::vector<int> & frozen_atoms = {},
	     const std::vector<int> & core_shells = {}  ) const{
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d2e;
      energy_and_derivs(E, D1E, true, empty_d2e, false, geom, ngbr,
			active_atoms, frozen_atoms, core_shells);
    }

    void d2e(REAL & E, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E,
	     Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E,
	     geometry<REAL> & geom,
	     const neighbours_table<REAL> & ngbr,
	     const std::vector<int> & active_atoms = {},
	     const std::vector<int> & frozen_atoms = {},
	     const std::vector<int> & core_shells = {}  ) const{
      energy_and_derivs(E, D1E, true, D2E, true, geom, ngbr,
			active_atoms, frozen_atoms, core_shells);
    }

    virtual bonding_table<REAL> distances() const { return bonding_table<REAL>(); }
    
    // fixme - temporary
    virtual ISTREAM & load(ISTREAM &) =0;

    // debug
    virtual void debug(){}

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    py::object py_energy(geometry<REAL> & geom,
			 const neighbours_table<REAL> & ngbr,
			 bool do_d1e, bool do_d2e,
			 const std::vector<int> & active_atoms,
			 const std::vector<int> & frozen_atoms,
			 const std::vector<int> & core_shells )
    {
      if (!do_d1e)
	return py::cast(energy(geom,ngbr, active_atoms, frozen_atoms, core_shells));
      else if (!do_d2e)
	{
	  REAL E = 0e0;
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E =
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(geom.nat(),3);
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d2e;
	  energy_and_derivs(E, D1E, true, empty_d2e, false, geom, ngbr,
			    active_atoms, frozen_atoms, core_shells);

	  //std::cout << "D1E before return\n";
	  //std::cout << D1E<< "\n";
	    
	  
	  return py::make_tuple(E,D1E);
	}
      else
	{
	  REAL E = 0e0;
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E=
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(geom.nat(),3);
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D2E =
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(3*geom.nat(),3*geom.nat());
	  energy_and_derivs(E, D1E, true, D2E, true, geom, ngbr, active_atoms, frozen_atoms, core_shells);
	  return py::make_tuple(E,D1E,D2E);	  
	}
    }

    static void py_export(py::module m, const char * pyname) {
      py::class_<classical_potential<REAL> >(m, pyname, py::dynamic_attr())
	.def("energy",    &SELF::py_energy, py::arg("geom"), py::arg("ngbr"),
	     py::arg("do_d1e") = false, py::arg("do_d2e") = false,
	     py::arg("active_atoms") = std::vector<int>(), py::arg("frozen_atoms") = std::vector<int>(),
	     py::arg("core_shells") = std::vector<int>(),
	     "Energy contribution from this potential")
	.def("distances", &SELF::distances, "Returns maximum distances of this potential in the form of bonding_table")
	.def_readwrite_static("too_close", & SELF::too_close)
	;
    }

#endif
    
  };

  template <class REAL>
  REAL classical_potential<REAL>::too_close = 1e-8;
  
  // ---------------------- Potential calculator ------------------

  template <class REAL>
  class mm_calculator{

    //std::vector< std::shared_ptr<classical_potential<REAL> > > pot;
    std::vector<classical_potential<REAL>*> pot;
    
    bool regions_defined;

    std::vector<int> active_atoms, frozen_atoms;
    std::vector<int> active_regions, frozen_regions;
    std::vector<int> core_shells;        
    
  public:

    mm_calculator(){ regions_defined = false;}

    void add_potential(classical_potential<REAL> & p)
    {
      pot.push_back(&p);
    }

    void set_active_atoms(const std::vector<int> & _atoms){active_atoms = _atoms; regions_defined = false;}

    void set_frozen_atoms(const std::vector<int> & _atoms){frozen_atoms = _atoms; regions_defined = false;}

    void set_active_regions(const std::vector<int> & _reg){active_regions = _reg; regions_defined = true;}

    void set_frozen_regions(const std::vector<int> & _reg){frozen_regions = _reg; regions_defined = true;}

    void set_core_shells(const std::vector<int> & _cs){core_shells = _cs;}

    void energy_and_derivs(REAL & E,
			   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E, bool do_d1e,
			   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E, bool do_d2e,
			   geometry<REAL> & geom)
    {
      std::cout << "mm_calculator::energy_and_derivs\n";
      if (regions_defined){
	// mark active and frozen atoms by corresponding regions
	//must be xgeometry
	if (!geom.is_xgeometry())
	  TypeError("Partition into regions is requested for the geometry with no region fiedls (not xgeometry)");
	xgeometry<REAL>  *xgeom = (xgeometry<REAL>*)(&geom);

	std::cout << xgeom->name << " " << xgeom->nat();
	for (int i=0; i<xgeom->nfields(); i++) std::cout << xgeom->field_name(i) << " " << (int)(xgeom->field_type(i));
	std::cout << "\n";
	
	std::vector<int> iregs;
	for (int i=0; i<xgeom->nfields(); i++)
	  if ( xgeom->field_name(i).substr(0,3) == "reg" && xgeom->field_type(i)& basic_types::type_int ){
	    iregs.push_back(i);
	    std::cout << xgeom->field_name(i) <<" " <<  i<< "\n";
	  }

	std::cout << "IREG:\n";
	for(int i:iregs)
	  std::cout << i << " ";
	std::cout << "\n";
	  
	active_atoms.clear();
	frozen_atoms.clear();
	for (int i=0; i<xgeom->nat(); i++){
	  bool active = false, frozen = false;
	  for (int k:iregs)
	    if (xgeom->field_types[k] ==  basic_types::type_int ){
	      int r = xgeom->template xfield<int>(k,i);
	      if (std::find(active_regions.begin(), active_regions.end(), r) !=
		  active_regions.end()) {
		active = true;
		break;
	      }
	    }
	    else if  (xgeom->field_types[k] ==  basic_types::type_int + basic_types::type_array){
	      std::vector<int> reg = xgeom->template xfield<std::vector<int> >(k,i);
	      if (have_common(reg,active_regions)){
		  active = true;
		  break;
		}
		}
	      for (int k:iregs)
		if (std::find(frozen_regions.begin(), frozen_regions.end(), xgeom->template xfield<int>(k,i)) !=
		    frozen_regions.end()) {
		  frozen = true;
		  break;
		}
	      if (active)
		active_atoms.push_back(i);
	      else if (frozen)
		frozen_atoms.push_back(i);
	    
	    }
	  
	}

	std::cout << "active atoms: ";
	for (int i:active_atoms)
	  std::cout << i << " ";
	std::cout << "\nfrozen atoms: ";
	for (int i:frozen_atoms)
	  std::cout << i << " ";
	std::cout << "\n";
      

	auto bt = bonding_table<REAL>();
	for (auto p:pot)
	  bt.merge(p->distances());

	std::cout << py::str(bt.to_dict()) << "\n";

	geom.build_typetable();
	neighbours_table<REAL>  ngbr(geom,bt);
	ngbr.build();

	for (auto p:pot){
	  std::cout << typeid(pot).name();
	  for (auto i:  active_atoms) std::cout << i  << " ";
	  std::cout << "\n";
	
	  p->energy_and_derivs(E, D1E, do_d1e, D2E, do_d2e, geom, ngbr, active_atoms, frozen_atoms, core_shells);
	}
      }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      py::object py_energy(geometry<REAL> & geom, bool do_d1e, bool do_d2e)
      {     
	if (!do_d1e){
	  REAL E = 0e0;
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d1e;
	  Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d2e;
	  energy_and_derivs(E, empty_d1e, false, empty_d2e, false, geom);
	  return py::cast(E);
	}
	else if (!do_d2e)
	  {
	    REAL E = 0e0;
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E =
	      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(geom.nat(),3);	    
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> empty_d2e;
	    energy_and_derivs(E, D1E, true, empty_d2e, false, geom);

	    //std::cout << "D1E before return\n";
	    //std::cout << D1E << "\n";
	  
	    return py::make_tuple(E,D1E);
	  }
	else
	  {
	    REAL E = 0e0;
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E =
	      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(geom.nat(),3);
	    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D2E = 
	      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(3*geom.nat(),3*geom.nat());
	    energy_and_derivs(E, D1E, true, D2E, true, geom);
	    return py::make_tuple(E,D1E,D2E);	  
	  }
      }

      static void py_export(py::module m, const char * pyname) {
	py::class_<mm_calculator<REAL> >(m, pyname)
	  .def(py::init<>())
	  .def("set_active_atoms",   &mm_calculator<REAL>::set_active_atoms)
	  .def("set_frozen_atoms",   &mm_calculator<REAL>::set_frozen_atoms)
	  .def("set_active_regions", &mm_calculator<REAL>::set_active_regions)
	  .def("set_frozen_regions", &mm_calculator<REAL>::set_frozen_regions)
	  .def_readwrite("active_atoms",  &mm_calculator<REAL>::active_atoms)
	  .def_readwrite("frozen_atoms",  &mm_calculator<REAL>::frozen_atoms)
	  .def("add_potential",      &mm_calculator<REAL>::add_potential)
	  .def("energy",  &mm_calculator<REAL>::py_energy,
	       py::arg("geom"), py::arg("do_d1e") = false, py::arg("do_d2e") = false )
	  .def("set_core_shells", &mm_calculator<REAL>::set_core_shells)
	  ;
      }
    
#endif
  };    

  // ---------------------------------------------------------------

  //template<class REAL>
  //inline Eigen::Matrix<REAL, 1, 3> eig3(const vector3<REAL> &v){
  //  return Eigen::Matrix<REAL, 1, 3>({v[0],v[1],v[2]});
  //}
      
  // ---------------------- Pair potentials ------------------------

  template <class REAL>
  class pair_potential : public classical_potential<REAL>{
  public:

    using classical_potential<REAL>::too_close;
    
    typedef pair_potential<REAL> SELF;
    
    using classical_potential<REAL>::distances;
    
    virtual void e12(REAL r, REAL & E, REAL & D1E, REAL & D2E, bool do_d1e = false, bool do_d2e = false) const{}
    
    /*
    virtual REAL e12(REAL r) const {}
    
    virtual REAL de12(REAL r) const {}

    virtual REAL d2e12(REAL r) const {}
    */
    
    STRING_EX atom1, atom2;
    REAL rmax;

    virtual void energy_and_derivs(REAL & E,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E, bool do_d1e,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E, bool do_d2e,
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
	      vector3<REAL> rvec = geom.pos(k) - geom.pos(i);
	      REAL r = rvec.norm();
	      if (r>rmax)
		continue;
	      REAL e, de, d2e;
	      e12(r,e,de,d2e,do_d1e,do_d2e);
	      E += e;
	      if (do_d1e && r > too_close)
		{
		  vector3<REAL> devec = de*rvec/r;
		  D1E.template block<1,3>(i,0) -= devec;
		  D1E.template block<1,3>(k,0) += devec;
		  if (do_d2e)
		    {
		      // fixme - implement this
		    }
		}
	    }
      std::cout << "E= "  << E << "\n";
    }
    
    virtual bonding_table<REAL> distances() const
    {
      bonding_table<REAL> b;
      b.set_pair(atom1,atom2,rmax);
      return b;
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    py::object py_e12(REAL r, bool do_d1e, bool do_d2e)
    {
      REAL E12, D1E12, D2E12;
      e12(r,E12, D1E12, D2E12, do_d1e, do_d2e);
      if (!do_d1e)
	return py::cast(E12);
      else if (!do_d2e)
	return py::make_tuple(E12,D1E12);
      else
	return py::make_tuple(E12,D1E12,D2E12);	
    }
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<SELF, classical_potential<REAL> >(m, pyname)
	.def("e12",      &SELF::py_e12, py::arg("r"), py::arg("do_d1e") = false, py::arg("do_d2e") = false,
	     "Energy and its derivatives for two atoms at distance r")
	//.def("de12",     &SELF::de12, "1st derivative of energy of two atoms at distance r")
	//.def("d2e12",    &SELF::d2e12, "2nd derivative of energy of two atoms at distance r")
	.def("atoms",  [](SELF & self) -> py::tuple {return py::make_tuple(self.atom1,self.atom2);})
	.def_readonly("rmax", &SELF::rmax)
	;
    }

#endif
  };
  
  // ---------------------- 3-particle potentials ------------------------
  
  template <class REAL>
  class potential_3body : public classical_potential<REAL>{
  public:
    
    using classical_potential<REAL>::too_close;
    using classical_potential<REAL>::distances;
    typedef potential_3body<REAL> SELF;
    
    STRING_EX atom1, atom2, atom3;
    REAL r12max, r13max;

    virtual void e123(REAL r12, REAL r13, REAL theta,
		      REAL &e, std::vector<REAL> & d1e, std::vector<REAL> & d2e,
		      bool do_d1e = false, bool do_d2e = false) const {}

    virtual void energy_and_derivs(REAL & E,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E, bool do_d1e,
				   Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E, bool do_d2e,
				   geometry<REAL> & geom,
				   const neighbours_table<REAL> & ngbr,
				   const std::vector<int> & active_atoms,
				   const std::vector<int> & frozen_atoms,
				   const std::vector<int> & core_shells) const
    {
      std::cout << "Three body energy called\n";
      
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

      for (int i=0; i<geom.nat(); i++)
	std::cout << i << " " << atomregs[i] << "\n";
      
      
      E = 0e0;
      // fixme - no need to rebuild ngbr table for each potential
      for (int i=0; i<geom.nat(); i++)
	if (atomregs[i] > 0 && geom.atom(i)==atom1)
	  for (int j1=0; j1<ngbr.n(i); j1++)
	    {
	      int j = ngbr(i,j1);
	      if ( atomregs[j] == 0 || geom.atom(j)!=atom2 )
		continue;
	      //if ( atom1==atom2 && i>=j)
		// fixme - what about periodic cells?
	      //	continue;
	      for (int k1=0; k1<ngbr.n(i); k1++)
		{
		  int k = ngbr(i,k1);
		  if ( atomregs[k] == 0 || geom.atom(k)!=atom3 )		      
		    continue;
		  if (atomregs[i] == 2 && atomregs[j] == 2 && atomregs[k] ==2)
		    continue;
		  //if ( atom1==atom3 && i>=k)
		  // continue;
		  if (atom2==atom3 && j>=k)
		    continue;
		  vector3<REAL> r12v = geom.pos(j) - geom.pos(i);
		  vector3<REAL> r13v = geom.pos(k) - geom.pos(i);
		  REAL
		    r12 = r12v.norm(),
		    r13 = r13v.norm(),
		    r12r13 = r12v.dot(r13v),
		    cos_tet = r12r13/(r12*r13);
		  if ( r12 > r12max || r13 > r13max )
		    continue;
		  REAL theta;
		  if (cos_tet>=1e0)
		    theta = 0e0;
		  else if (cos_tet<-1e0)
		    theta = pi;
		  else
		    theta = std::acos(cos_tet);
		  vector3<REAL> n12 = r12v/r12, n13 = r13v/r13;
		  
		  REAL e;
		  std::vector<REAL> d1e123(3), d2e123(6);
		  e123(r12,r13,theta, e, d1e123, d2e123, do_d1e, do_d2e);
		  E += e;
		  if (do_d1e)
		    {
		      vector3<REAL>
			dtet_r2 = n13 - n12.dot(n13)*n12,
			dtet_r3 = n12 - n12.dot(n13)*n13;
		      REAL sin = std::sin(theta);
		      if (sin < too_close){
			dtet_r2 *= 0e0;
			dtet_r3 *= 0e0;
		      }
		      else{
			dtet_r2 /= - r12*sin;
			dtet_r3 /= - r13*sin;
		      }
		      //Eigen::Vector3d
		      vector3<REAL>
			DEr2 = d1e123[0]*n12 + d1e123[2]*dtet_r2,
			DEr3 = d1e123[1]*n13 + d1e123[2]*dtet_r3;
		      D1E.template block<1,3>(j,0) += DEr2;
		      D1E.template block<1,3>(k,0) += DEr3;
		      D1E.template block<1,3>(i,0) -= DEr2 + DEr3;
		      if (do_d2e)
			// fixme - implement this
			{
			}
		    }
		}
	    }
    }
    
    virtual bonding_table<REAL> distances() const
    {
      bonding_table<REAL> b;
      b.set_pair(atom1,atom2,r12max);
      b.set_pair(atom1,atom3,r13max);
      return b;
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
  
  py::object py_e123(REAL r12, REAL r13, REAL theta, bool do_d1e, bool do_d2e)
    {
      REAL E;
      std::vector<REAL> d1e(3), d2e(6);
      e123(r12, r13, theta, E, d1e, d2e, do_d1e, do_d2e);
      if (!do_d1e)
	return py::cast(E);
      else if (!do_d2e)
	return py::make_tuple(E, d1e);
      else
	return py::make_tuple(E, d1e, d2e);	
    }
    
    static void py_export(py::module m, const char * pyname) {
      py::class_<SELF, classical_potential<REAL> >(m, pyname)
	.def("e123",      &SELF::py_e123, py::arg("r12"), py::arg("r13"), py::arg("theta"),
	     py::arg("do_d1e") = false, py::arg("do_d2e") = false,
	     "Energy and its derivatives for three atoms positioned by r12, r13 and theta")
	.def("atoms",  [](SELF & self) -> py::tuple {return py::make_tuple(self.atom1,self.atom2,self.atom3);})
	.def_readonly("r12max", &SELF::r12max)
	.def_readonly("r13max", &SELF::r13max)
	;
    }
  
#endif

  };

  
  template <class REAL>
  std::istream& operator>>(std::istream &is, pair_potential<REAL> &pp)
  {
    return pp.load(is);
  } 

  
};

#endif
