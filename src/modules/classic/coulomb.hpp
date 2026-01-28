#ifndef _QPP_COULOMB_H
#define  _QPP_COULOMB_H

#include <geom/xgeom.hpp>
#include <symm/gen_cell.hpp>
#include <consts.hpp>
#include <Eigen/Dense>

/*
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
//#include <pybind11/pybind11.h>
//#include <pybind11/operators.h>
#include <pybind11/stl.h>
//PYBIND11_MAKE_OPAQUE(const std::vector<std::vector<float> > &);
//PYBIND11_MAKE_OPAQUE(const std::vector<std::vector<double> > &);
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#endif
*/
namespace qpp{

  template<class REAL>
  bool find_core_shells(std::vector<int> & pairs, xgeometry<REAL> & geom, REAL maxdistance){
    bool res = false;
    pairs.resize(geom.nat());
    std::fill(pairs.begin(), pairs.end(), -1);
    geom.build_typetable();
    int n = geom.typetable()->n_atom_types();
    std::vector<int> pair_types(n,-1);
    for (int i=0; i<n; i++)
      for (int j=0; j<i; j++){
	STRING_EX s1 = geom.typetable()->atomic_type(i);
	STRING_EX s2 = geom.typetable()->atomic_type(j);
	int k = common_begin(s1,s2);
	if (k>0 and oneof<STRING_EX>( s1.substr(k),{"cor","core"}))
	  if (oneof<STRING_EX>( s2.substr(k), {"shl","shel","shell"})){
	    pair_types[i] = j;
	    pair_types[j] = i;
	    res = true;
	  }
	if (k>0 and oneof<STRING_EX>( s2.substr(k),{"cor","core"}))
	  if (oneof<STRING_EX>( s1.substr(k), {"shl","shel","shell"})){
	    pair_types[i] = j;
	    pair_types[j] = i;
	    res = true;
	  }
      }

    /*
    std::cout << "find_cores_shells\n";
    for (int t=0; t<geom.n_types(); t++)
      std::cout << t << " " << geom.atom_of_type(t) << "\n";
    for (int p:pair_types) std::cout << p << " ";
    std::cout << "\n";
    */
    
    bonding_table<REAL> bt;
    bt.default_distance = maxdistance;
    neighbours_table<REAL> nt(geom,bt);
    nt.build();
    
    for (int i = 0; i<geom.nat(); i++)
      if (pair_types[geom.typetable()->type(i)]!=-1){
	int t = pair_types[geom.typetable()->type(i)];
	
	//std::cout << i << " ";
	//for (int k=0; k<nt.n(i); k++) std::cout << nt(i,k) << geom.atom(nt(i,k)(0));
	//std::cout << "\n";
	
	for (int k=0; k<nt.n(i); k++){
	  int j = nt(i,k)(0);
	  if (geom.typetable()->type(j)==t){
	    pairs[i] = j;
	    pairs[j] = i;
	    res = true;
	  }
	}
      }	
    return res;
  }
  
  template<class REAL>
  void mmcharges(std::vector<REAL> & charge, xgeometry<REAL> & geom, const std::vector<int> & regions){
    int i1 = -1, i2 = -1, ir1 = -1, ir2 = -1;
    for (int i = 0; i<geom.nfields(); i++){
      if (geom.field_name(i)=="qmm1")
	i1 = i;
      if (geom.field_name(i)=="qmm2")
	i2 = i;
      if (geom.field_name(i)=="reg1")
	ir1 = i;
      if (geom.field_name(i)=="reg2")
	ir2 = i;
    }
    if (i1==-1 || i2 == -1)
      throw std::runtime_error("\"qmm1\" and \"qmm2\" fields not found\n");
    charge.resize(geom.nat());
    for (int j = 0; j < geom.nat(); j++){
      REAL q = 0e0;
      if (oneof<int>(geom.template xfield<int>(ir1,j),regions))
	q += geom.template xfield<REAL>(i1,j);
      if (oneof<int>(geom.template xfield<int>(ir2,j),regions))
	q += geom.template xfield<REAL>(i2,j);
      charge[j] = q;
    }
  }
    
  template <class REAL>
  struct coulomb_point_charges{

    // length units in which the coordinates of charges are provided
    // energy units in which to return energy
    STRING_EX length_units, energy_units;

    REAL length_scale, energy_scale;
    
    std::vector<REAL> charges;
    std::vector<vector3<REAL> > coords;
    std::vector<int> core_shell;
    
    REAL too_close;
    REAL core_shell_distance;
    bool shell_model;

    typedef coulomb_point_charges<REAL> SELF;
    
    coulomb_point_charges()
    {
      length_units = "bohr";
      energy_units = "au";
      too_close = globals::too_close;
      setscales();
      core_shell_distance=globals::too_close;
    }

    void setscales(){
      STRING_EX lu = tolower(length_units),
	eu = tolower(energy_units);
      if (lu == "bohr")
	length_scale = 1e0;
      else if (lu == "angstrom")
	length_scale = 1e0/ang_to_bohr;
      else if (lu == "nm" || lu == "nanometer")
	length_scale = 10e0/ang_to_bohr;
      if (eu == "au" || eu == "hartree" )
	energy_scale = 1e0;
      else if (eu == "ev")
	energy_scale = hartree_to_ev;
      else if (eu == "kjm")
	energy_scale = 2625.5;
      else if (eu == "kcalm")
	energy_scale = 627.5;
    }

    coulomb_point_charges( xgeometry<REAL> & geom, const std::vector<int> & reglist = {},
			   REAL _too_close = (REAL)globals::too_close){
      length_units = "bohr";
      energy_units = "au";
      too_close = _too_close;

      std::vector<int> atoms;
      if (reglist.size()==0)
	for (int i=0; i < geom.nat(); i++) atoms.push_back(i);
      else {
	std::vector<int> iregs;
	for (int i=0; i<geom.nfields(); i++)
	  if ( geom.field_name(i).substr(0,3) == "reg" && geom.field_type(i) == basic_types::type_int )
	    iregs.push_back(i);
	for (int i=0; i<geom.nat(); i++){
	  bool active = false;
	  for (int k:iregs) {
	    int r = geom.template xfield<int>(k,i);
	    if (std::find(reglist.begin(), reglist.end(), r) != reglist.end()) {
	      active = true;
	      break;
	    }
	  }
	  if (active)
	    atoms.push_back(i);
	}	  
      }

      std::vector<REAL> mmcharge;
      mmcharges(mmcharge,geom,reglist);
      
      for (int i:atoms){
	charges.push_back(mmcharge[i]);
	coords.push_back(geom.pos(i));
      }

      //core_shell_distance = _core_shell_distance;
      shell_model = find_core_shells(core_shell, geom, core_shell_distance);
    }

    coulomb_point_charges(const std::vector<std::vector<REAL> > &v,
			  REAL _too_close = (REAL)globals::too_close){
      length_units = "bohr";
      energy_units = "au";

      shell_model = false;
      
      for (const auto & l:v){
	charges.push_back(l[0]);
	coords.push_back({l[1],l[2],l[3]});
	core_shell.push_back(-1);
      }
      /*
      if (_core_shell.size() == coords.size()){
      shell_model = true;*/
      too_close = _too_close;
      setscales();
      core_shell_distance=globals::too_close;
	//      }
      
    }

    REAL Eint;
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1Eint, D2Eint;
    
    struct calc_thread{
      int n1, n2;
      bool finished;
      bool d1e, d2e;
      std::shared_ptr<coulomb_point_charges<REAL>> self, c2;

      REAL Epart;
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1Epart, D2Epart;

      calc_thread( std::shared_ptr<coulomb_point_charges<REAL>> _self,
		   std::shared_ptr<coulomb_point_charges<REAL>> _c2,
		   int _n1, int _n2, bool _d1e=false, bool _d2e = false
		 )
	: n1(_n1), n2(_n2), d1e(_d1e), d2e(_d2e),
	  self(_self), c2(_c2), finished(false){}

      
      void operator()(){
	int N = self->charges.size();
	REAL gscale = self->energy_scale/(self->length_scale*self->length_scale);
	REAL escale = self->energy_scale/self->length_scale;
	bool selfenergy = (c2 == self);
	std::cout << "coulomb thread " << n1 << " " << n2 << " is here\n";
	Epart = 0e0;
	D1Epart = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(N,3);
	D2Epart = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(N*3,N*3);
	
	for (int i=n1; i<=n2; i++){
	  vector3<REAL> g = 0e0;
	  for (int j = 0; j<c2->charges.size(); j++)
	    {
	      if (selfenergy && self->core_shell.size()>=N && self->core_shell[i] == j)
		continue;
	      vector3<REAL> r = self->coords[i] - c2->coords[j];
	      REAL R = r.norm();
	      if (R > self->too_close){
		if (d1e)
		  g -= self->charges[i]*c2->charges[j]*r/(R*R*R);
		Epart += escale*self->charges[i]*c2->charges[j]/R;
	      }
	    }
	  if (d1e)
	    D1Epart.template block<1,3>(i,0) =  gscale*g;
	  }
	finished = true;
      }
    };

    void calc_parallel(coulomb_point_charges<REAL> & c2, bool d1e=false, bool d2e = false){
      setscales();
      std::vector<calc_thread> thd;
      int N = charges.size();
      int ncores = globals::ncores;
      for (int n=0; n<ncores;n++)	
	thd.push_back(calc_thread(std::shared_ptr<SELF>(this,[](SELF*p){}),
				  std::shared_ptr<SELF>(&c2,[](SELF*p){}),
				  n*N/ncores, (n+1)*N/ncores-1,d1e,d2e));
      for (int i =0; i< ncores; i++) {
	std::thread th(std::ref(thd[i]));
	th.detach();
      }

      while(true){
	bool finished = true;
	for (int i =0; i< ncores; i++)
	  if (!thd[i].finished){
	    finished = false;
	    break;
	  }
	if (finished) break;
      	std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }

      Eint = thd[0].Epart;
      D1Eint = thd[0].D1Epart;
      D2Eint = thd[0].D2Epart;
      for (int i =1; i< ncores; i++) {
	Eint +=thd[i].Epart;
	D1Eint += thd[i].D1Epart;
	D2Eint += thd[i].D2Epart;
      }
    }

    void calculate_seq(coulomb_point_charges<REAL> & c2, bool d1e=false, bool d2e = false){
      setscales();
      REAL escale = energy_scale/(length_scale);
      REAL gscale = energy_scale/(length_scale*length_scale);
      bool selfenergy = (&c2 == this);
      Eint = 0e0;
      int N = charges.size();
      D1Eint = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(N,3);
      D2Eint = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(N*3,N*3);
      for (int i = 0; i<charges.size(); i++)
	{
	  vector3<REAL> g = 0e0;
	  for (int j = 0; j<c2.charges.size(); j++)
	    {
	      if (selfenergy && core_shell.size()>=charges.size() && core_shell[i] == j)
		continue;
	      vector3<REAL> r = coords[i] - c2.coords[j];
	      REAL R = r.norm();
	      if (R > too_close){
		if (d1e)
		  g -= charges[i]*c2.charges[j]*r/(R*R*R);
		Eint += escale*charges[i]*c2.charges[j]/R;
	      }
	    }
	  if (d1e)
	    D1Eint.template block<1,3>(i,0) =  gscale*g;
	}
    }

    void calculate( coulomb_point_charges<REAL> & c2, bool d1e=false, bool d2e = false){
      if (globals::ncores==1)
	calculate_seq(c2,d1e,d2e);
      else{
	std::cout << "Parallel coulomb interaction calc on ncores= " << globals::ncores<< "\n";
	calc_parallel(c2,d1e,d2e);
      }
    }
    //------------------------------------------------------------------------------------        
    // REAL interaction_energy_seq(const coulomb_point_charges<REAL> & c2)
    // {
    //   calculate_seq(c2);
    //   return Eint;
    // }

    // void interaction_gradients_seq(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E,
    // 			       const coulomb_point_charges<REAL> & c2)
    // {
    //   calculate_seq(c2, true);
    //   D1E = D1Eint;
    // }

    // void interaction_hessian_seq(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E,
    // 				 const coulomb_point_charges<REAL> & c2)
    // {
    //   calculate_seq(c2, true, true);
    //   D2E = D2Eint;
    // }

    //-------------------------------------------------------------------    
    REAL interaction_energy( coulomb_point_charges<REAL> & c2){
      Eint = 0e0;
      calculate(c2,false,false);
      return Eint;
    }

    void interaction_gradients(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E,
			        coulomb_point_charges<REAL> & c2)      
    {
      Eint = 0e0;
      D1Eint = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(charges.size(),3);
      calculate(c2,true,false);
      D1E=D1Eint;
    }

    void interaction_hessian(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E,
			      coulomb_point_charges<REAL> & c2){
      Eint = 0e0;
      D1Eint = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(charges.size(),3);
      D2Eint = Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>::Zero(charges.size()*3,charges.size()*3);
      calculate(c2,true,true);
      D2E = D2Eint;
    }

      
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    
    coulomb_point_charges(const py::list &l, REAL _too_close = (REAL)globals::too_close){
      length_units = "bohr";
      energy_units = "au";

      shell_model = false;
      too_close = _too_close;
      setscales();
      core_shell_distance=globals::too_close;
      
      for (const auto  &sl:l)
	if (!py::isinstance<py::list>(sl))
	  TypeError("Required here the list of the form [[q1,x1,y1,z1],[q2,x2,y2,z2],...]");
	else
	{
	  const auto  &rl = py::cast<py::list>(sl);
	  REAL q = py::cast<REAL>(rl[0]);
	  REAL x = py::cast<REAL>(rl[1]);
	  REAL y = py::cast<REAL>(rl[2]);
	  REAL z = py::cast<REAL>(rl[3]);
	  charges.push_back(q);
	  coords.push_back({x,y,z});	  
      }
      /*
      if (_core_shell.size() == coords.size()){
      shell_model = true;*/
            }
      
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>
    py_interaction_gradients( coulomb_point_charges<REAL> & c2){
      int n = charges.size();
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E(n,3);
      interaction_gradients(D1E,c2);
      return D1E;
    }
    
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>
    py_interaction_hessian( coulomb_point_charges<REAL> & c2){
      int n = charges.size();
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D2E(n*3,n*3);
      interaction_hessian(D2E,c2);
      return D2E;
    }
    
    static void py_export(py::module m, const char * pyname) {
      
      py::class_<coulomb_point_charges<REAL>, std::shared_ptr<coulomb_point_charges<REAL>> >(m,
					 (std::string("coulomb_point_charges_")+pyname).c_str())
	.def(py::init<>())
	.def(py::init< xgeometry<REAL > &,const std::vector<int>, REAL >(), py::arg("geom"),
	     py::arg("reglist") = std::vector<int>(), py::arg("too_close")=0e0)
	//.def(py::init< xgeometry<REAL> &,
	//   const std::vector<int> &, bool >(), py::arg("geom"),
	//     py::arg("reglist") = std::vector<int>(), py::arg("shell_model")=0e0)
	//.def(py::init< xgeometry<REAL > &,
	//     const std::vector<int> &, bool >(), py::arg("geom"),
	//     py::arg("reglist") = std::vector<int>(), py::arg("shell_model")=0e0)
	//.def(py::init< std::vector<std::vector<REAL> > , REAL>(),
	.def(py::init<const py::list &, REAL>(),
	     py::arg("charges"), py::arg("_too_close") = 0e0 )
	.def("interaction_energy", & coulomb_point_charges<REAL>::interaction_energy )
	.def("interaction_gradients", & coulomb_point_charges<REAL>::py_interaction_gradients )
	.def("interaction_hessian", & coulomb_point_charges<REAL>::py_interaction_hessian )
	.def("setscales", & coulomb_point_charges<REAL>::setscales )
	.def_readwrite("too_close", & coulomb_point_charges<REAL>::too_close)
	.def_readwrite("core_shell_distance", & coulomb_point_charges<REAL>::core_shell_distance)
	.def_readwrite("charges", & coulomb_point_charges<REAL>::charges )
	.def_readwrite("coords", & coulomb_point_charges<REAL>::coords )
	.def_readwrite("core_shell", & coulomb_point_charges<REAL>::core_shell )
	.def_readwrite("shell_model", & coulomb_point_charges<REAL>::shell_model )
	.def_readwrite("length_units", & coulomb_point_charges<REAL>::length_units )
	.def_readwrite("length_scale", & coulomb_point_charges<REAL>::length_scale )
	.def_readwrite("energy_units", & coulomb_point_charges<REAL>::energy_units )
	.def_readwrite("energy_scale", & coulomb_point_charges<REAL>::energy_scale )
	;
    }
    
#endif
    
  };
  /*
  template <class REAL>
  REAL coulomb_point_charges<REAL>::too_close = 1e-6;
  */
  // void coulomb_field( std::vector<REAL> & pot, std::vector<vector3<REAL> > & field, bool do_field,
  // 		      std::vector<vector3<REAL> > & where, const xgeometry<REAL> & geom)
  // {
  //   for (int i=0; i < where.size(); i++)
  //     for (int j=0; j < geom.size(); j++)
  // 	{
  // 	  vector3<REAL> Rji = where[i] - geom.pos(j);
  // 	  REAL q = geom.charge(j);
  // 	  REAL RRji = Rji.norm();
  // 	  pot[i] += q/RRji;
  // 	  if (do_field)
  // 	    field[i] += q*Rji/(RRji*RRji*RRji);
  // 	}
  // }

  // REAL coulomb_interaction_energy(const xgeometry<REAL> & geom1, const xgeometry<REAL> & geom2)
  // {
  //   std::vector<REAL> R;
  //   for (int i = 0; i < geom1.size(); i++)
  //     R.push_back(geom1.pos(i));
  //   std::vector<REAL> pot(geom1.size());
  //   std::vector<vector3<REAL> > empty_field;
  //   coulomb_field(pot,empty_field,false,R,geom2);
  //   REAL E=0;
  //   for (int i = 0; i < geom1.size(); i++)
  //     E += pot[i]*geom1.charge(i);
  //   return E;
  // }
  
};

#endif
