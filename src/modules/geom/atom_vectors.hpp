#ifndef QPP_ATOM_VECTORS_H
#define QPP_ATOM_VECTORS_H

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
//#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif

#include<geom/xgeom.hpp>
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<qpp::vector3<float>> >);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<qpp::vector3<double>> >);
#endif
namespace qpp{

  template <class REAL>
  struct geom_atom_vectors{
    typedef geom_atom_vectors<REAL> SELF;
    
    geometry<REAL> * bound_geom;
    std::vector<vector3<REAL> > vectors;
    std::vector<vector3<REAL>> colors;

    geom_atom_vectors(geometry<REAL> &__bound_geom = nullptr){
      bound_geom = &__bound_geom;
      if (bound_geom){
	vectors.resize(bound_geom -> nat());
	std::fill(vectors.begin(), vectors.end(), vector3<REAL>(0e0));
	colors.resize(bound_geom -> nat());
	std::fill(colors.begin(), colors.end(), vector3<REAL>(0e0));
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
	py_vec.bind(this);
	py_crs.bind(this);
#endif
      }
    }

    geom_atom_vectors(const std::vector<vector3<REAL> > & __vectors,
		      geometry<REAL> & __bound_geom) :
      vectors(__vectors), bound_geom(&__bound_geom) {
    }

    bool is_valid() const{
      if (bound_geom)
	if (vectors.size() == bound_geom->nat() )
	  return true;
      return false;
    }

    vector3<REAL> start_pos(int i) const{
      return bound_geom->pos(i);
    }

    vector3<REAL> end_pos(int i, REAL scale = 1e0) const{
      return bound_geom->pos(i) + scale*vectors[i];
    }
    
    vector3<REAL> start_pos(int i, const index &I) const{
      return bound_geom->pos(i,I);
    }

    vector3<REAL> end_pos(int i, const index &I, REAL scale = 1e0) const{
      return bound_geom->pos(i,I) + scale*vectors[i];
    }
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    vector3<REAL> py_getvec(int i) {
      int n =  bound_geom->nat();
      if (i<0) i+= n;
      if (i<0 || i>=n)
	IndexError("atom_vector: index out of range");
      return vectors[i];
    }

    void py_setvec(int i, const vector3<REAL> &v)
    {
      int n =  bound_geom->nat();
      if (i<0) i+= n;
      if (i<0 || i>=n)
	IndexError("atom_vector: index out of range");
      vectors[i] = v;
    }
    
    vector3<REAL> py_getcrs(int i) {
      int n =  bound_geom->nat();
      if (i<0) i+= n;
      if (i<0 || i>=n)
	IndexError("atom_vector: index out of range");
      return colors[i];
    }

    void py_setcrs(int i, const vector3<REAL> &v)
    {
      int n =  bound_geom->nat();
      if (i<0) i+= n;
      if (i<0 || i>=n)
	IndexError("atom_vector: index out of range");
      colors[i] = v;
    }

    py_indexed_property< SELF, vector3<REAL>, int,&SELF::py_getvec, &SELF::py_setvec> py_vec;
    py_indexed_property< SELF, vector3<REAL>, int, &SELF::py_getcrs, &SELF::py_setcrs> py_crs;
    
    static void py_export(py::module m, const char * pyname){
      std::string  xvecname = fmt::format("{0}_{1}",pyname,"vec_indexed_property");
      std::string  xcolname = fmt::format("{0}_{1}",pyname,"clrs_indexed_property");
      py_indexed_property<  SELF, vector3<REAL>, int,
			    &SELF::py_getvec, &SELF::py_setvec>::py_export(m, xvecname.c_str());
      py_indexed_property<  SELF, vector3<REAL>, int,
      		    &SELF::py_getcrs, &SELF::py_setcrs>::py_export(m, xcolname.c_str());
      py::class_<SELF,std::shared_ptr<SELF>>(m,pyname)
	.def(py::init<geometry<REAL>&>())
	.def(py::init<const std::vector<vector3<REAL>> &,geometry<REAL>&>())
	.def_readwrite("vectors", &SELF::py_vec)
	.def_readwrite("colors", &SELF::py_crs);
    }
   
#endif
  };

}

#endif
