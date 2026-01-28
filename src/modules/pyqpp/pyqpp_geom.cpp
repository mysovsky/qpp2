#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
//#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#pragma pop_macro("slots")
namespace py = pybind11;

#include <pyqpp/pyqpp.hpp>
#include <geom/geom.hpp>
#include <symm/gen_cell.hpp>
#include <geom/atom_vectors.hpp>

using namespace qpp;

template<class REAL>
void py_geom_export(py::module m, const char * pyname) {

  qpp::geometry<REAL>::py_props(m, pyname);
  py::class_<qpp::geometry<REAL>,
             std::shared_ptr<qpp::geometry<REAL>> >(m, pyname, py::dynamic_attr())
    .def(py::init<int,   const STRING_EX&>(),
	 py::arg("dim"), py::arg("__name") = "")
    //.def(py::init<qpp::periodic_cell<REAL>&, const STRING_EX&>(),
    //      py::arg("CELL"), py::arg("__name") = "")
    .def("add_observer",    & qpp::geometry<REAL>::add_observer )
    .def("remove_observer", & qpp::geometry<REAL>::remove_observer )
    .def("n_observers", & qpp::geometry<REAL>::n_observers )
    .def("build_types",     & qpp::geometry<REAL>::build_typetable)
    .def("set_typetable", & qpp::geometry<REAL>::set_typetable )
    .def("set_cell",   &qpp::geometry<REAL>::set_cell)
    .def("add",        &qpp::geometry<REAL>::py_add1)
    .def("add",        &qpp::geometry<REAL>::py_add2)
    .def("add",        &qpp::geometry<REAL>::py_add3)
    .def("erase",      &qpp::geometry<REAL>::py_erase)
    .def("erase",      &qpp::geometry<REAL>::py_erase2)
    .def("insert",     &qpp::geometry<REAL>::py_insert1)
    .def("insert",     &qpp::geometry<REAL>::py_insert2)
    .def("insert",     &qpp::geometry<REAL>::py_insert3)
    .def("clear",      &qpp::geometry<REAL>::clear)
    .def("sort",       &qpp::geometry<REAL>::sort)
    .def("replace",    &qpp::geometry<REAL>::replace)
    .def("pos",        &qpp::geometry<REAL>::py_pos1)
    .def("pos",        &qpp::geometry<REAL>::py_pos2)
    .def("pos",        &qpp::geometry<REAL>::py_pos3)
    .def("r",          &qpp::geometry<REAL>::py_pos1)
    .def("r",          &qpp::geometry<REAL>::py_pos2)
    .def("r",          &qpp::geometry<REAL>::py_pos3)
    .def("nat",        &qpp::geometry<REAL>::nat)
    .def("__len__",    &qpp::geometry<REAL>::nat)
    .def("__getitem__",&qpp::geometry<REAL>::py_getitem)
    .def("__setitem__",&qpp::geometry<REAL>::py_setitem)
    //.add_property("cell",
    //make_getter(& qpp::geometry<REAL,CELL>::cell, return_value_policy<reference_existing_object>()),
      //make_getter(& qpp::geometry<REAL,CELL>::cell, return_self<>()),
      //& qpp::geometry<REAL,CELL>::py_setcell)
      .def_readwrite("name",   &qpp::geometry<REAL>::name)
      .def_readwrite("cell",   &qpp::geometry<REAL>::cell)
      //.def_readonly( "dim",    &qpp::geometry<REAL,CELL>::DIM)
      .def( "dim",    &qpp::geometry<REAL>::DIM)
      .def_readwrite("atom",   &qpp::geometry<REAL>::py_atoms)
      .def_readwrite("coord",  &qpp::geometry<REAL>::py_coords)
      .def_readwrite("shadow", &qpp::geometry<REAL>::py_shadow)
      .def_readwrite("x",      &qpp::geometry<REAL>::py_x)
      .def_readwrite("y",      &qpp::geometry<REAL>::py_y)
      .def_readwrite("z",      &qpp::geometry<REAL>::py_z)
    .def_readwrite("observer",        & qpp::geometry<REAL>::py_observer)

      .def("typetable",     &qpp::geometry<REAL>::typetable)
    //  .def("atom_of_type", &qpp::geometry<REAL>::atom_of_type)
    //.def("type",         &qpp::geometry<REAL>::py_typeofatom)
    //  .def("n_types",      &qpp::geometry<REAL>::n_atom_types)
    /*
      .def_readwrite("auto_update_types",
                     &qpp::geometry<REAL>::auto_update_types,
                     "bool auto_update_types: Whether to update atomic type "
                     "table after any atom manipulations")
    */
      .def_readwrite("frac",
                     &qpp::geometry<REAL>::frac,
                     "bool frac: Whether to undestand the atomic"
                     " coordinates as fractional")

    .def_property("tol_geom",
		  [](const qpp::geometry<REAL>&self){return *self.tol_geom;},
		  [](qpp::geometry<REAL>&self, REAL v){*self.tol_geom=v;})
      ;
}


template<class REAL>
void py_observer_export(py::module m, const char * pyname) {
  py::class_<qpp::geometry_observer<REAL>,
      qpp::py_geometry_observer<REAL> >(m, pyname)
      .def(py::init<>())
      .def("added",    &qpp::geometry_observer<REAL>::added)
      .def("inserted", &qpp::geometry_observer<REAL>::inserted)
      .def("changed",  &qpp::geometry_observer<REAL>::changed)
      .def("erased",   &qpp::geometry_observer<REAL>::erased)
      .def("shaded",   &qpp::geometry_observer<REAL>::shaded)
      ;
}

void pyqpp_geom_export(py::module m) {
  geometry_type_table<float>::py_export(m,"type_table_f");
  py_geom_export<float>(m, "geometry_f");
  //  py_geom_export<float, qpp::gen_cell<float, qpp::matrix3<float> > >(m, "geometry_pgf");
  // py_geom_export<float, qpp::gen_cell<
  //  float, qpp::matrix3<float>   > >(m, "geometry_cgf");
  py_observer_export<float>(m, "gobserver_f");

#ifdef PYTHON_EXP_EXT
  geometry_type_table<double>::py_export(m,"type_table_d");
  py_geom_export<double>(m, "geometry_d");
  // py_geom_export<double,qpp::gen_cell<
  //    double, qpp::matrix3<double> > >(m, "geometry_pgd");
  // py_geom_export<double,qpp::gen_cell<
  //    double,qpp::rotrans<double>  > >(m, "geometry_cgd");
  py_observer_export<double>(m, "gobserver_d");
#endif

  py::enum_<qpp::before_after>(m, "geom_change_state")
      .value("before", qpp::before)
      .value("after",  qpp::after)
      ;
  qpp::geom_atom_vectors<float>::py_export(m,"atom_vector_f");
#ifdef PYTHON_EXP_EXT
  qpp::geom_atom_vectors<double>::py_export(m,"atom_vector_d");
#endif
  
}
