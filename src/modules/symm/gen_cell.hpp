#ifndef QPP_GCELL_H
#define QPP_GCELL_H

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif

#include <mathf/lace3d.hpp>
#include <symm/index.hpp>
#include <symm/transform.hpp>
#include <symm/groups.hpp>
#include <symm/cell.hpp>
#include <consts.hpp>
#include <stdexcept>
#include <cmath>

namespace qpp{

  /*
    The supercell concept generalization for the geometry class looks like:

  template< int DIM, class REAL>
  class CELL{
  public:

    vector3d<REAL> transform(const vector3d<REAL> & r, const index<DIM> & I) const;

    bool within(const vector3d<REAL> & r) const;

    vector3d<REAL> reduce(const vector3d<REAL> & r) const;

    // find high symmetry point within "radius" distance from given point "r"
    // makes sence for rotational symmetries
    vector3d<REAL> symmetrize(const vector3d<REAL> & r, REAL radius) const;

    // fractional to cartesian and vice versa transformation
    // makes sence only for periodic translational cells
    vector3d<REAL> frac2cart(const vector3d<REAL> & r) const;

    vector3d<REAL> cart2frac(const vector3d<REAL> & r) const;

  };

  Any class having these methods - can be used as a "supercell"
  for some geometry

  */

  template<class REAL,class TRANSFORM>
  struct gen_cell : public periodic_cell<REAL>{
    genform_group<TRANSFORM> G;
    //typedef rotrans<REAL> TRANSFORM;
  /*
    using genform_group<TRANSFORM>::_begin;
    using genform_group<TRANSFORM>::_end;
    using genform_group<TRANSFORM>::DIM;
    using genform_group<TRANSFORM>::generators;
    using genform_group<TRANSFORM>::generate;
  */
    typedef gen_cell<REAL,TRANSFORM> SELF;

    void init_default()
    {
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_gen.bind(this);
#endif
    }

    gen_cell(int dim=0) : genform_group<TRANSFORM>(dim)
    { init_default(); }

    gen_cell(const std::vector<TRANSFORM> & g,const index & __begin,
                     const index & __end) :
      G(g, __begin, __end)
    { init_default();}

    gen_cell(const std::vector<TRANSFORM> & g) :
      G(g)
    { init_default();}

    gen_cell(const genform_group<TRANSFORM> & _G) :
      G(_G)
    { init_default(); }

    gen_cell(const gen_cell<REAL,TRANSFORM> & _G) :
      G(_G.G)
    { init_default();}

    virtual vector3<REAL> transform(const vector3<REAL> & r,
                                    const index & I) const{
      return G(I)*r;
    }
    void auto_orders(int d=-1) {
      if (d==-1)
	G.auto_orders();
      else
	G.auto_order(d);
    }
    /*    
    void auto_orders() {
    }
    */
    int get_dim(){return G.get_dim();}

    void set_dim(int D){
      G.set_dim(D);
    }



    virtual bool within(const vector3<REAL> & r) const{
    // Answers the question whether r belongs to the unit cell
      // Does not make sense for gen_cell
      return true;
    }

    /*
    // Brings the point r into the volume of unit cell
    // by translations
    // unit cell is defined as parallelepiped with one vertex in
    // the coordinate origin
    // the others are pointed by v[0],v[1],v[2] vectors

    // Does not make sense for gen_cell */
    virtual vector3<REAL> reduce(const vector3<REAL> & r) const{
      return r;
    }

    // find high symmetry point within "radius" distance from given point "r"
    // makes sence for rotational symmetries
    virtual vector3<REAL> symmetrize(const vector3<REAL> & r, REAL radius) const{
      std::vector<vector3<REAL> > points;

      //std::cerr << "symmetrize for " << r << "\n";

      points.push_back(r);
      for (iterator I(G._begin,G._end); !I.end(); I++){
          if (I==index::D(G.DIM).all(0)) continue;

          vector3<REAL> r1 = transform(r,I);
          for (const vector3<REAL> & pt : points)
            if ((r1-pt).norm() < radius ){
	      std::cerr << r1 << "close\n";
                points.push_back(r1);
                break;
              }
        }
      vector3<REAL> res = vector3<REAL>::Zero();
      for (const vector3<REAL> & pt : points)
        res += pt;
      res /= points.size();
      return res;
    }

    // fractional to cartesian and vice versa transformation
    // makes sence only for periodic translational cells
    virtual vector3<REAL> frac2cart(const vector3<REAL> & r) const{
      // Works for DIM==3 and DIM==2
      // In the case DIM==2 the third coordinate ( z=r(2) ) is
      // orthogonal to both translation vectors

      // Does not make sense for gen_cell
      return r;
    }

    virtual vector3<REAL> cart2frac(const vector3<REAL> & r) const{
      // cartesian to fractional
      // works for DIM==3 and DIM==2
      // In the case DIM==2 the third coordinate ( z=r(2) ) is
      // orthogonal to both translation vectors

      // Does not make sense for gen_cell
      return r;
    }

    virtual void write(std::basic_ostream<CHAR_EX,TRAITS> &os, int offset=0) const{

      // fixme - implement this
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    gen_cell(const py::list & _G,
                     const index & __begin,
                     const index & __end):
                     genform_group<TRANSFORM>(){

      init_default();

      G.DIM = py::len(_G);
      for (int i=0; i<G.DIM; i++){
          /*  if (!py::cast<TRANSFORM>(G[i])())
            TypeError("gen_cell constructor: expected
            symmetry operation");*/
          G.generators.push_back(py::cast<TRANSFORM>(_G[i]));
        }
      G._begin = __begin;
      G._end   = __end;
    }

    gen_cell(const py::list & _G):
      G(){
      init_default();

      G.DIM = py::len(_G);
      for (int i=0; i<G.DIM; i++)
        {
          /*  if (!py::cast<TRANSFORM>(G[i])())
            TypeError("gen_cell constructor:
            expected symmetry operation");*/
          G.generators.push_back(py::cast<TRANSFORM>(_G[i]));
        }
      G._begin = G._end = index::D(G.DIM).all(0);
    }

    TRANSFORM py_getgen(int i) {
      if (i<0) i+=G.DIM;
      if (i<0 || i>=G.DIM) IndexError("gen_cell::index out of range");
      return G.generators[i];
    }

    void py_setgen(int i, const TRANSFORM & t) {
      if (i<0) i+=G.DIM;
      if (i<0 || i>=G.DIM) IndexError("gen_cell::index out of range");
      G.generators[i] = t;
    }

    py_indexed_property<SELF, TRANSFORM,
    int, &SELF::py_getgen, &SELF::py_setgen> py_gen;

    TRANSFORM py_call(index I){
      return G(I);
    }

    TRANSFORM py_call1(py::tuple I){
      //std::cout << "() from tuple called\n";
      return G(index(I));
    }

    py::list py_generate(){
      std::vector<TRANSFORM> v;
      G.generate(v);
      py::list l;
      for (const auto & g : v)
        l.append(g);
      return l;
    }

    virtual index begin() const{
      return G.begin();
    }
    virtual index end() const{
      return G.end();
    }

    static void py_export(py::module m, const char * pyname){
      std::string indexPropName = fmt::format("{0}_{1}",pyname,"indexed_property");
      py_indexed_property<SELF, TRANSFORM, int, & SELF::py_getgen,
          &SELF::py_setgen>::py_export(m, indexPropName.c_str());
      py::class_<gen_cell<REAL,TRANSFORM>, periodic_cell<REAL>,
		 std::shared_ptr<gen_cell<REAL,TRANSFORM> >>(m, pyname)
	//.def(py::init<const gen_cell<REAL,TRANSFORM> >())
	//.def(py::init<const py::list &, const index &, const index&>())
	.def(py::init<const py::list &>())
	.def("begin",  & SELF::begin )
	.def("end",    & SELF::end )
	.def_readwrite("generator", &SELF::py_gen)
	.def_property("dim", &SELF::get_dim, &SELF::set_dim)
	.def("transform",   &SELF::transform)
	.def("within",      &SELF::within)
	.def("cart2frac",   &SELF::cart2frac)
	.def("frac2cart",   &SELF::frac2cart)
	.def("auto_orders", &SELF::auto_orders, py::arg("d")=-1)
	//	.def("auto_orders", &SELF::auto_orders)
	.def("generate",    &SELF::py_generate)
	.def("__call__",    &SELF::py_call)
	.def("__call__",    &SELF::py_call1)
	.def_readwrite("G",&SELF::G)
	;
    }

#endif

  };

}

#endif
