#include <pyqpp/pyqpp.hpp>
#include <classic/potentials.hpp>
#include <classic/pair_potentials.hpp>
#include <classic/potentials_3b.hpp>
#include <classic/coulomb.hpp>
#include <symm/gen_cell.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

namespace py = pybind11;
PYBIND11_MAKE_OPAQUE(const std::vector<std::vector<float> > &);
PYBIND11_MAKE_OPAQUE(const std::vector<std::vector<double> > &);
//py::implicitly_convertible<const std::vector<std::vector<float> > &, py::list>;
//py::implicitly_convertible<const std::vector<std::vector<double> > &, py::list>;

template<class REAL>
std::vector<int> py_find_core_shells(qpp::xgeometry<REAL> & geom, REAL maxdistance){
  std::vector<int> res;
  if (qpp::find_core_shells(res,geom,maxdistance))
    return res;
  else
    return std::vector<int>();
}


template<class REAL >
void py_potentials_export (py::module m, const char * pyname) {

  // generic potential
  qpp::classical_potential<REAL>::py_export(m,(std::string("classical_potential_") + pyname).c_str());
  qpp::pair_potential<REAL>::py_export(m,(std::string("pair_potential_") + pyname).c_str());
  qpp::potential_3body<REAL>::py_export(m,(std::string("potential_3body") + pyname).c_str());
  
  qpp::mm_calculator<REAL>::py_export(m,(std::string("mm_calculator_") + pyname).c_str());

  m.def("find_core_shells", &py_find_core_shells<REAL>);
  
  py::module pp = m.def_submodule("pp");
  
  qpp::buckingham_potential<REAL>::py_export(pp,(std::string("buckingham_") + pyname).c_str());
  qpp::buckingham4_potential<REAL>::py_export(pp,(std::string("buckingham4_") + pyname).c_str());
  qpp::morse_potential<REAL>::py_export(pp,(std::string("morse_") + pyname).c_str());
  qpp::cutcoulomb_potential<REAL>::py_export(pp,(std::string("cutcoulomb_") + pyname).c_str());
  qpp::spring_potential<REAL>::py_export(pp,(std::string("spring_") + pyname).c_str());
  qpp::three_harmonic<REAL>::py_export(pp,(std::string("three_harm_") + pyname).c_str());
}

void pyqpp_potentials_export (py::module m) {
  
  py_potentials_export<float>(m, "f");
  //  py_potentials_export<float, qpp::gen_cell<float, qpp::matrix3<float> > >(m, "pgf");
  //py_potentials_export<float, qpp::gen_cell<float, qpp::rotrans<float> > >(m, "cgf");

  qpp::coulomb_point_charges<float>::py_export(m,"f");
   
#ifdef PYTHON_EXP_EXT
  
  py_potentials_export<double>(m, "d");
  //  py_potentials_export<double, qpp::gen_cell<double, qpp::matrix3<double> > >(m, "pgd");
  //py_potentials_export<double, qpp::gen_cell<double, qpp::rotrans<double> > >(m, "cgd");

  qpp::coulomb_point_charges<double>::py_export(m,"d");
  
  
#endif
}
