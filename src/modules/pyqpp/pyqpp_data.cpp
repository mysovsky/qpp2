#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/iostream.h>
#pragma pop_macro("slots")

namespace py = pybind11;

#include <pyqpp/pyqpp.hpp>
#include <geom/geom.hpp>
#include <symm/gen_cell.hpp>

using namespace qpp;

void pyqpp_data_export(py::module m) {  
  py::class_<qpp::globals>(m,"globals")
    .def_readwrite_static("ncores", &qpp::globals::ncores)
    .def_readwrite_static("too_close", &qpp::globals::too_close)
    ;
}
