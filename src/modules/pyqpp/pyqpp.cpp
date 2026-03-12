#include <pyqpp/pyqpp.hpp>
//#include <pybind11/embed.h> 

//#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
//PYBIND11_EMBEDDED_MODULE(pyqpp, m) {
//#else
PYBIND11_MODULE(pyqpp, m) {
  //#endif

  //std::cout << "qpp-rebuild version\n";
  try {
    pyqpp_linalg_export(m);
    pyqpp_data_export(m);
    pyqpp_opaque_types_export(m);

    pyqpp_math_export(m);
    pyqpp_cell_export(m);
    pyqpp_geom_export(m);
    pyqpp_xgeom_export(m);
    //pyqpp_basis_ecp_export(m);
    //pyqpp_molecule_export(m);
    pyqpp_shape_export(m);
    pyqpp_neighbours_export(m);
    pyqpp_potentials_export(m);
    pyqpp_builders_export(m);
    //pyqpp_potentials_export(m);
    //pyqpp_autosymm_export(m);
    pyqpp_io_export(m);
    //pyqpp_gmsio_export(m);
    pyqpp_ptable_export(m);
    //pyqpp_ccd_export(m);

    //pyqpp_experimental_export(m);
  } catch (const std::exception &e) {
    std::cerr << "FATAL: Exception during pyqpp module initialization: " << e.what() << std::endl;
    throw;
  } catch (...) {
    std::cerr << "FATAL: Unknown exception during pyqpp module initialization" << std::endl;
    throw;
  }
}
