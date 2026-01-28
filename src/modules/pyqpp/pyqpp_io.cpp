#include <pyqpp/pyqpp.hpp>

#include <io/geomio.hpp>
#include <symm/gen_cell.hpp>
#include <io/fdstream.hpp>
#include <io/vasp_io.hpp>

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

//vasp io start
template<class REAL>
//void py_read_vasp_poscar (int fd, qpp::geometry<REAL> & geom) {
void py_read_vasp_poscar (py::object & f, qpp::geometry<REAL> & geom) {

  boost::fdistream inp(qpp::fileno(f));
  qpp::read_vasp_poscar(inp,geom);

}

template< class REAL >
void py_write_vasp_poscar (py::object & f, qpp::geometry<REAL> & geom) {

  boost::fdostream out(qpp::fileno(f));
  qpp::write_vasp_poscar(out,geom);

}

//vasp io end
template<class REAL>
void py_read_xyz (py::object & f, qpp::geometry<REAL> & geom) {

  boost::fdistream inp(qpp::fileno(f));
  qpp::read_xyz(inp,geom);

}

template<class REAL>
void py_read_xyzq (py::object & f, qpp::xgeometry<REAL> & geom) {

  boost::fdistream inp(qpp::fileno(f));
  qpp::read_xyzq(inp,geom);

}

template< class REAL>
void py_write_xyz (py::object & f, const qpp::geometry<REAL> & geom) {

    boost::fdostream out(qpp::fileno(f));
    qpp::write_xyz(out,geom);

}

template< class REAL>
void py_write_xyzq (py::object & f, const qpp::xgeometry<REAL> & geom) {

    boost::fdostream out(qpp::fileno(f));
  qpp::write_xyzq(out,geom);

}

template< class REAL>
void py_export_ioxyz (py::module m) {

  m.def("read_xyz",  &py_read_xyz<REAL>);
  m.def("read_xyzq", &py_read_xyzq<REAL>);
  m.def("write_xyz", &py_write_xyz<REAL>);
  m.def("write_xyzq",&py_write_xyzq<REAL>);

}

template< class REAL >
void py_export_vasp_io (py::module m) {

  m.def("read_vasp_poscar", &py_read_vasp_poscar<REAL>);
  //m.def("read_vasp_outcar_md", &py_read_vasp_outcar_md<REAL, CELL> );
  m.def("write_vasp_poscar", &py_write_vasp_poscar<REAL>);

}

void pyqpp_io_export (py::module m) {

  py::module io = m.def_submodule("io");
  py_export_ioxyz<float>(io);
  py_export_vasp_io<float>(io);

#ifdef PYTHON_EXP_EXT
  py_export_ioxyz<double>(io);
  py_export_vasp_io<double>(io);
#endif

  /*
    py_export_ioxyz<float, qpp::generalized_cell<float,
        qpp::matrix3<float>  > >(m);

    py_export_ioxyz<double,
        qpp::generalized_cell<double, qpp::matrix3<double> > >(m);

    py_export_ioxyz<float,
        qpp::generalized_cell<float,  qpp::rotrans<float>   > >(m);

    py_export_ioxyz<double,
        qpp::generalized_cell<double, qpp::rotrans<double>  > >(m);
*/
}



#endif
