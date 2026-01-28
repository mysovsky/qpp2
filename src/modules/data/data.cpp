#include "data.hpp"
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif
namespace qpp{

  int globals::ncores = 1;

  double globals::too_close = 1e-6;
  /*
  template<> int datum::get<int>() const {
    if (tp!=type_int) TypeError("Type mismatch in datum conversion to int");
    if (ptr) return *pi;
    else return i;
  }

  template<> float datum::get<float>() const {
    if (tp==type_float){
        if (ptr) return *pf;
        else return f;
      }
    else if (tp==type_double) {
        if (ptr) return *pd;
        else return d;
      }
    else
      throw std::invalid_argument("Type mismatch in datum conversion to float");
  }

  template<> double datum::get<double>() const {
    if (tp==type_float){
        if (ptr) return *pf;
        else return f;
      }
    else if (tp==type_double){
        if (ptr) return *pd;
        else return d;
      }
    else
      throw std::invalid_argument("Type mismatch in datum conversion to double");
  }

  template<> bool datum::get<bool>() const {
    if (tp!=type_bool) throw std::invalid_argument("Type mismatch in datum conversion to bool");
    if (ptr) return *pb;
    else return b;
  }

  template<> STRING_EX datum::get<STRING_EX>() const {
    if (tp!=type_string) throw std::invalid_argument("Type mismatch in datum conversion to STRING");
    if (ptr) return *ps;
    else return s;
  }

  // -------------------------------------------------------------------

  template<> int& datum::ref<int>() {
    if (tp!=type_int) throw std::invalid_argument("Type mismatch in datum conversion to int");
    if (ptr) return *pi;
    else return i;
  }

  template<> float& datum::ref<float>() {
    if (tp!=type_float) throw std::invalid_argument("Type mismatch in datum conversion to float");
    if (ptr) return *pf;
    else return f;
  }

  template<> double& datum::ref<double>() {
    if (tp!=type_double) throw std::invalid_argument("Type mismatch in datum conversion to double");
    if (ptr) return *pd;
    else return d;
  }

  template<> bool& datum::ref<bool>() {
    if (tp!=type_bool) throw std::invalid_argument("Type mismatch in datum conversion to bool");
    if (ptr) return *pb;
    else return b;
  }

  template<> STRING_EX& datum::ref<STRING_EX>() {
    if (tp!=type_string) throw std::invalid_argument("Type mismatch in datum conversion to STRING");
    if (ptr) return *ps;
    else return s;
  }

  datum::datum() {
    tp = type_bool; ptr = false; b = false;
  }

  datum::datum( int    _i){set(_i);}
  datum::datum( float  _f){set(_f);}
  datum::datum( double _d){set(_d);}
  datum::datum( bool   _b){set(_b);}
  datum::datum( const STRING_EX & _s){set(_s);}
  datum::datum( const char * _s){set(_s);}

  datum::datum( int    *_i){set(_i);}
  datum::datum( float  *_f){set(_f);}
  datum::datum( double *_d){set(_d);}
  datum::datum( bool   *_b){set(_b);}
  datum::datum( STRING_EX *_s){set(_s);}

  datum::datum(const datum & v) {
    ptr = v.ptr;
    tp = v.tp;
    if (!ptr)
      switch(tp){
        case type_bool:   b = v.b; break;
        case type_int:    i = v.i; break;
        case type_float:  f = v.f; break;
        case type_double: d = v.d; break;
        case type_string: s = v.s; break;
        default: break;
        }
    else
      switch(tp){
        case type_bool:   pb = v.pb; break;
        case type_int:    pi = v.pi; break;
        case type_float:  pf = v.pf; break;
        case type_double: pd = v.pd; break;
        case type_string: ps = v.ps; break;
        default: break;
        }
  }
  */
  Boolarray:: Boolarray(){}

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
  Boolarray::Boolarray(const py::list &l){
    for (auto it:l)
      {
	bool b = py::cast<bool>(it);
	arr.push_back(b);
      }
  }
#endif
  bool Boolarray::getitem(int i) const{return arr[i];}

  void Boolarray::setitem(int i, bool b){arr[i] = b;}
 
  void Boolarray::append(bool b){arr.push_back(b);}

  std::string Boolarray::str(){
    std::string s="[";
    bool first = true;    
    for (bool i:arr){
      if (first)
	first = false;
      else
	s += ", ";
      if (i)
	s += "True";
      else
	s += "False";
    }
    s += "]";
    return s;
  }
  void Boolarray::clear(){
    arr.clear();
  }
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
  
  void Boolarray::py_export(py::module m, const char * pyname){
    py::class_<Boolarray>(m,pyname)
      .def(py::init<>())
      .def_readwrite("arr", & Boolarray::arr)
      .def("append", &Boolarray::append)
      .def("insert", [](Boolarray & self, int i, bool b){self.arr.insert(self.arr.begin()+i, b);})
      .def("erase", [](Boolarray & self, int i){self.arr.erase(self.arr.begin()+i);})
      .def("__getitem__", &Boolarray::getitem)
      .def("__setitem__", &Boolarray::setitem)
      .def("__len__", [](Boolarray &self){return self.arr.size();} )
      .def("__str__", &Boolarray::str)
      .def("__repr__", &Boolarray::str);
    
  }
#endif
};
