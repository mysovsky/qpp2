#include <pyqpp/pyqpp.hpp>
#include <symm/autosymm.hpp>
#include <symm/shoenflis.hpp>
#include <symm/point_groups.hpp>
#include <symm/permut.hpp>

template<class REAL>
void def_autosymm (py::module m, const char * pyname) {

  m.def("best_transform", &qpp::py_best_transform<REAL>);
  m.def("best_axes", &qpp::best_axes<REAL>);
  m.def("unitarize", &qpp::unitarize<REAL>);
  //m.def("rotate_transform", &qpp::rotate_pair<REAL>);
  m.def("analyze_transform", &qpp::py_analyze_transform<REAL>);
  m.def("find_point_symm", &qpp::find_point_symm<REAL>);
  m.def("bravais_point_group", &qpp::bravais_point_group<REAL>);
  //m.def("bravais_point_group1", &qpp::bravais_point_group1<REAL>);
  m.def("find_cryst_symm", &qpp::py_find_cryst_symm2<REAL>,
	py::arg("group"), py::arg("unit_cell"),
	py::arg("rtol") = qpp::geometry<REAL>::tol_geom_default);
  m.def("find_cryst_symm", &qpp::py_find_cryst_symm1<REAL>,
	py::arg("group"), py::arg("P"), py::arg("unit_cell"),
	py::arg("rtol") = qpp::geometry<REAL>::tol_geom_default);
  m.def("find_cryst_symm", &qpp::py_find_cryst_symm3<REAL>,
	py::arg("unit_cell"),
	py::arg("rtol") = qpp::geometry<REAL>::tol_geom_default);
  //m.def("find_point_subgroups", &qpp::py_find_point_subgroups<REAL,true>);
  //m.def("find_point_subgroups", &qpp::py_find_point_subgroups<REAL>);
  //m.def("find_point_subgroups", &qpp::py_find_point_subgroups2<REAL,true>);
  //m.def("find_point_subgroups", &qpp::py_find_point_subgroups2<REAL,false>);
  m.def("find_translations", &qpp::py_find_translations<REAL>,
	py::arg("g1"), py::arg("g2"), py::arg("cell"), py::arg("R"),
	py::arg("return_permut") = false);
  //m.def("finitize_point_group", &qpp::finitize_point_group<REAL>);
  m.def("reconstruct_point_group", &qpp::py_reconstruct_point_group<REAL>);
  m.def("reconstruct_cryst_group", &qpp::py_reconstruct_cryst_group<REAL>);
  //m.def("complete_group", &qpp::complete_group<qpp::matrix3<REAL> >);
  //m.def("complete_group", &qpp::complete_group<qpp::rotrans<REAL,true> >);
  m.def("point_group_symbol", &qpp::point_group_symbol<REAL>);
  //
  m.def("fix4_cryst_group", & qpp::py_fix4_cryst_group<REAL>);
  m.def("rotrans_sub", & qpp::rotrans_sub<REAL>);
  m.def("rotrans_diff", [] (const qpp::rotrans<REAL> & R1,
			    const qpp::rotrans<REAL> & R2){
    REAL rdiff, tdiff;
    qpp::rotrans_diff(rdiff, tdiff, R1, R2);
    return py::make_tuple(rdiff,tdiff);
  });
  //
  
  std::string shname = fmt::format("{0}_{1}","shnfl",pyname);
  
  py::class_<qpp::shnfl<REAL> >(m, shname.c_str())
    .def("Cs",  & qpp::shnfl<REAL>::Cs )
    .def("Ci",  & qpp::shnfl<REAL>::Ci )
    .def("Cn",  & qpp::shnfl<REAL>::Cn )
    .def("Cnv", & qpp::shnfl<REAL>::Cnv )
    .def("Cnh", & qpp::shnfl<REAL>::Cnh )
    .def("Dn",  & qpp::shnfl<REAL>::Dn )
    .def("Dnh", & qpp::shnfl<REAL>::Dnh )
    .def("Dnd", & qpp::shnfl<REAL>::Dnd )
    .def("S2n", & qpp::shnfl<REAL>::S2n )
    .def("T",   & qpp::shnfl<REAL>::T )
    .def("Td",  & qpp::shnfl<REAL>::Td )
    .def("Th",  & qpp::shnfl<REAL>::Th )
    .def("O",   & qpp::shnfl<REAL>::O )
    .def("Oh",  & qpp::shnfl<REAL>::Oh )
    .def("I",   & qpp::shnfl<REAL>::I )
    .def("Ih",  & qpp::shnfl<REAL>::Ih )
    .def("group", & qpp::shnfl<REAL>::group )
    .def("groups_by_order", & qpp::shnfl<REAL>::groups_by_order );
  std::string fpname = fmt::format("{0}_{1}",shname,"fingerprint");
  qpp::shnfl<REAL>::fingerprint::py_export(m,fpname.c_str());

  m.def("pg_approx_find",   &qpp::pg_approx_find<REAL>);
  m.def("pg_approx_mutab",  &qpp::pg_approx_multab<REAL>);
  m.def("pg_max_order",     &qpp::pg_max_order<REAL>,
	py::arg("G"), py::arg("angle_error") = 8*qpp::matrix3<REAL>::tol_equiv);

  //  std::string sbname = fmt::format("{0}_{1}","subspace3",pyname);  
  //qpp::subspace3<REAL>::py_export(m,sbname.c_str());
  //m.def("invariant_subspace",  &qpp::invariant_subspace<REAL>);

  std::string pgaxname = fmt::format("{0}_{1}","point_group_axes",pyname);
  py::class_<qpp::point_group_axes<REAL> >(m, pgaxname.c_str())
    .def(py::init<>())
    .def(py::init<const qpp::point_group_axes<REAL> &>())
    .def(py::init<const qpp::array_group< qpp::matrix3<REAL> > &>())
    .def_readwrite("axes",          &qpp::point_group_axes<REAL>::axes )
    .def_readwrite("planes",        &qpp::point_group_axes<REAL>::planes )
    .def_readwrite("orders",        &qpp::point_group_axes<REAL>::orders )
    .def_readwrite("rotoinversion", &qpp::point_group_axes<REAL>::rotoinversion )
    .def_readwrite("inversion",     &qpp::point_group_axes<REAL>::inversion )
    ;  
}

void pyqpp_autosymm_export (py::module m) {

  def_autosymm<float>(m, "f");

#ifdef PYTHON_EXP_EXT
  def_autosymm<double>(m, "d");
#endif

  qpp::permutation::py_export(m);

  py::class_<qpp::Bool>(m, "Bool")
    .def(py::init<>())
    .def(py::init<bool>())
    .def("__bool__",
	 [](const qpp::Bool & self) -> bool
	 {return bool(self);} )
    .def("__str__", [](const qpp::Bool & self) -> std::string
	 {if (self) return "true"; else return "false";})
    .def("__repr__", [](const qpp::Bool & self) -> std::string
	 {if (self) return "true"; else return "false";})
    ;
}
