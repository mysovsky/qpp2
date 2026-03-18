#include <pyqpp/pyqpp.hpp>
#include <grow/builders.hpp>
#include <grow/manyfold.hpp>
#include <geom/ngbr.hpp>
#include <symm/gen_cell.hpp>
//#include <pybind11/stl.h>
#include <pybind11/functional.h>

template<class REAL>
void def_builders4(py::module m){
  //  m.def("replicate",qpp::py_replicate1<REALDST,CELLDST, REALSRC,CELLSRC> );
  //  m.def("replicate",qpp::py_replicate2<REALDST,CELLDST, REALSRC,CELLSRC> );
  m.def("fill",qpp::py_fill1<REAL> );
  m.def("fill",qpp::py_fill2<REAL> );
  m.def("fill",qpp::py_fill3<REAL> );
}
/*
template<class REAL>
void def_builders3(py::module m){
  def_builders4< REAL>(m);
  def_builders4< REALDST, CELLDST, REALSRC,
      qpp::gen_cell<REALSRC,
      qpp::matrix3<REALSRC> > >(m);
  def_builders4< REALDST, CELLDST, REALSRC,
      qpp::gen_cell<REALSRC,
      qpp::rotrans<REALSRC> > >(m);
}

template<class REALDST, class CELLDST>
void def_builders2(py::module m){
  def_builders3< REALDST, CELLDST, float>(m);
  def_builders3< REALDST, CELLDST, double>(m);
}

template<class REALDST>
void def_builders1(py::module m){
  def_builders2< REALDST, qpp::periodic_cell<REALDST> >(m);
  def_builders2< REALDST, qpp::gen_cell<REALDST,
      qpp::matrix3<REALDST> > >(m);
  def_builders2< REALDST, qpp::gen_cell<REALDST,
      qpp::rotrans<REALDST> > >(m);
}
*/
template<class REAL, class SYMM>
void def_unique(py::module m, const std::string & kname){

  m.def("unique", qpp::py_unique1<REAL,SYMM>,
  py::arg("n_images"), py::arg("ugeom"), py::arg("nugeom"),
  py::arg("group"), py::arg("begin"), py::arg("end"),
  py::arg("key") =
  (std::function<REAL(const qpp::geometry<REAL> &, int)>)
  ([](const qpp::geometry<REAL> &g, int i) -> REAL
  { return std::sqrt(std::pow(g.pos(i)(1),2) + std::pow(g.pos(i)(2),2)); }),
  py::arg("eps") = qpp::geometry<REAL>::tol_geom_default);

  m.def("unique", qpp::py_unique2<REAL,SYMM>,
  py::arg("n_images"), py::arg("ugeom"), py::arg("nugeom"), py::arg("group"),
  py::arg("key") =
  (std::function<REAL(const qpp::geometry<REAL> &, int)>)
  ([](const qpp::geometry<REAL> &g, int i) -> REAL
  { return std::sqrt(std::pow(g.pos(i)(1),2) + std::pow(g.pos(i)(2),2)); }),
  py::arg("eps") = qpp::geometry<REAL>::tol_geom_default);

  m.def("unique", qpp::py_unique3<REAL>,
  py::arg("n_images"), py::arg("ugeom"), py::arg("nugeom"),
  py::arg("key") =
  (std::function<REAL(const qpp::geometry<REAL> &, int)>)
  ([](const qpp::geometry<REAL> &g, int i) -> REAL
  { return std::sqrt(std::pow(g.pos(i)(1),2) + std::pow(g.pos(i)(2),2)); }),
  py::arg("eps") = qpp::geometry<REAL>::tol_geom_default);

  m.def("unique", qpp::py_unique4<REAL,SYMM>,
  py::arg("ugeom"), py::arg("nugeom"), py::arg("group"),
  py::arg("key") =
  (std::function<REAL(const qpp::geometry<REAL> &, int)>)
  ([](const qpp::geometry<REAL> &g, int i) -> REAL
  { return std::sqrt(std::pow(g.pos(i)(1),2) + std::pow(g.pos(i)(2),2)); }),
  py::arg("eps") = qpp::geometry<REAL>::tol_geom_default);

  m.def("unique", qpp::py_unique5<REAL>,
  py::arg("ugeom"), py::arg("nugeom"), py::arg("key") =
  (std::function<REAL(const qpp::geometry<REAL> &, int)>)
  ([](const qpp::geometry<REAL> &g, int i) -> REAL
  { return std::sqrt(std::pow(g.pos(i)(1),2) + std::pow(g.pos(i)(2),2)); }),
  py::arg("eps") = qpp::geometry<REAL>::tol_geom_default);

}

template<class REAL>
void def_unique2(py::module m, const std::string & kname){

  /*
  std::cout << "registered " << "keyfunctype_"+kname << "\n";
  std::cout << typeid(REAL).name() << " " << typeid(NUCELL).name() << " " << typeid(NUCELL).name() << "\n";
  */
  py::class_<std::function<REAL(const qpp::geometry<REAL> &, int)> >(m,("keyfunctype_"+kname).c_str())
    .def(py::init<>());

  def_unique<REAL, qpp::periodic_cell<REAL>               >(m, kname+"t");
  def_unique<REAL, qpp::gen_cell<REAL, qpp::matrix3<REAL> >>(m, kname+"p");
  def_unique<REAL, qpp::gen_cell<REAL, qpp::rotrans<REAL> >>(m, kname+"c");
}

template<class REAL>
void def_unique3(py::module m, const std::string & kname){
  def_unique2<REAL, qpp::periodic_cell<REAL> >(m, kname+"t");
  def_unique2<REAL, qpp::gen_cell<REAL, qpp::matrix3<REAL> > >(m, kname+"p");
  def_unique2<REAL, qpp::gen_cell<REAL, qpp::rotrans<REAL> > >(m, kname+"c");
}

template<class REAL>
void def_paramsufr(py::module m){
  qpp::parametric_surface<float>:: py_export(m, "parametric_surface_f");
  qpp::orthoparametric_surface<float>:: py_export(m, "orthoparametric_surface_f");
  qpp::parametric_sphere<float>:: py_export(m, "parametric_sphere_f");
  qpp::parametric_torus<float>:: py_export(m, "parametric_torus_f");
}

void pyqpp_builders_export (py::module m) {

  m.def("treat_crowd", qpp::treat_crowd<float>);
  /*
  m.def("treat_crowd", qpp::treat_crowd<float,
                       qpp::gen_cell<float,
                       qpp::matrix3<float> > >);
  m.def("treat_crowd", qpp::treat_crowd<float,
                       qpp::gen_cell<float,
                       qpp::rotrans<float>  > >);
  */
  def_builders4<float>(m);
  def_unique2<float>(m, "f");

#ifdef PYTHON_EXP_EXT
  m.def("treat_crowd", qpp::treat_crowd<double >);
  /*
  m.def("treat_crowd", qpp::treat_crowd<double,
                       qpp::gen_cell<double,
                       qpp::matrix3<double> > >);
  m.def("treat_crowd", qpp::treat_crowd<double,
                       qpp::gen_cell<double,
                       qpp::rotrans<double> > >);
  */
  def_builders4<double>(m);
  def_unique2<double>(m, "d");
  qpp::parametric_surface<float>:: py_export(m, "parametric_surface_f");
  qpp::orthoparametric_surface<float>:: py_export(m, "orthoparametric_surface_f");
  qpp::parametric_sphere<float>:: py_export(m, "parametric_sphere_f");
  qpp::parametric_torus<float>:: py_export(m, "parametric_torus_f");

  qpp::parametric_surface<double>:: py_export(m, "parametric_surface_d");
  qpp::orthoparametric_surface<double>:: py_export(m, "orthoparametric_surface_d");
  qpp::parametric_sphere<double>:: py_export(m, "parametric_sphere_d");
  qpp::parametric_torus<double>:: py_export(m, "parametric_torus_d");
#endif

}
