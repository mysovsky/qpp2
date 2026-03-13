#ifndef QPP_GEOM_EXTRAS_H
#define QPP_GEOM_EXTRAS_H

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
#include <cassert>
#include <data/types.hpp>
#include <data/errors.hpp>
#include <data/data.hpp>
#include <geom/geom.hpp>
#include <io/strfun.hpp>
#include <io/simplefun.hpp>
#include <vector>
#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <variant>
#include <typeinfo>
#include <data/types.hpp>
#include <data/errors.hpp>
#include <data/data.hpp>
#include <geom/geom.hpp>
#include <io/strfun.hpp>
//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float> >);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double> >);  
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int> >);  
//PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Bool> >);

PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<float> > >);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<double> > >);  
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<int> > >);  
//PYBIND11_MAKE_OPAQUE(std::vector<unsigned char>);
//PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<Bool> > >);
#endif

namespace qpp {
  
  template<class T>
  void vecinsert(std::vector<T> & v, int i, const T& elem){
    v.insert(v.begin()+i,elem);
  }

  template<class T>
  void vecerase(std::vector<T> & v, int i){
    v.erase(v.begin()+i);
  }

  template<class T>
  struct XFTYPE{
    typedef T xftype;
  };
  
  template<>
  struct XFTYPE<bool>{
    typedef unsigned char xftype;
  };

  template<>
  struct XFTYPE<const char*>{
    typedef STRING_EX xftype;
  };
  
  template<class REAL>
  ///
  /// \brief  Geometry with extra fields - for storing any additional data
  ///
  class xgeometry : public geometry<REAL> {
  public:
    using geometry<REAL>::_atm;
    using geometry<REAL>::_crd;
    using geometry<REAL>::_shadow;
    //using geometry<REAL>::reorder_types;
    using geometry<REAL>::size;
    using geometry<REAL>::has_observers;
    using geometry<REAL>::observers;
    using geometry<REAL>::shadow;
    using geometry<REAL>::nat;
    using geometry<REAL>::atom;
    using geometry<REAL>::coord;
    using geometry<REAL>::erase;
    using geometry<REAL>::_typetable;
    //using geometry<REAL>::fieldtypes;
    typedef std::variant<REAL, int, bool, STRING_EX, std::vector<REAL>,
    			 std::vector<int>, std::vector<STRING_EX>, Boolarray> fieldtypes;
    typedef periodic_cell<REAL> CELL;
    //private:

   
  public:
    typedef std::variant<std::vector<REAL>, std::vector<int>, std::vector<unsigned char>,
			 std::vector<STRING_EX>,
			 std::vector<std::vector<REAL>>, std::vector<std::vector<int>>,
			 std::vector<std::vector<STRING_EX>>,std::vector<Boolarray>> xtrafieldtype;
    std::vector<unsigned char> _additive;
    std::vector<STRING_EX> field_names;
    std::vector<int> field_types;
    //tmp  
  public:
    xtrafieldtype& newcolumn(int t){
      //std::cout << "newcolumn " << t << "\n";
      if (t&type_array){
	if (t&type_real){
	  return *new xtrafieldtype(std::vector<std::vector<REAL>>());
	}
	else if (t&type_int)
	  //std::cout << "vector(int)\n";
	  return *new xtrafieldtype(std::vector<std::vector<int>>());	
	else if (t&type_bool)
	  return *new xtrafieldtype(std::vector<Boolarray>());
	else if (t&type_string)
	  return *new xtrafieldtype(std::vector<std::vector<STRING_EX>>());	
      }
      else{
	//	std::cout << "newline not array " << t <<"\n";
	if (t&type_real!=0){
	  //std::cout << "still real " << t<< "\n";
	  return *new xtrafieldtype(std::vector<REAL>());
	}
	else if (t==type_int)
	  return *new xtrafieldtype(std::vector<int>());
	else if (t==type_bool){
	  auto nn= new xtrafieldtype(std::vector<unsigned char>());
	  //std::cout << "new bool "<< nn-> index() << "\n";
	  return *nn;}
	else if (t==type_string)
	  return *new xtrafieldtype(std::vector<STRING_EX>());
       }
      return  *new xtrafieldtype;
    }

  public:
    typedef xgeometry<REAL> SELF;
    
    std::vector<xtrafieldtype> _xfields;

    int xftype(const STRING_EX& s){
      auto fld = split(s,"()");
      if (fld.size()==1){
	if (s=="s" || s=="str")
	  return type_string;
	else if (s=="i" || s=="int")
	  return type_int;
	else if (s=="bool" || s=="b")
	  return type_bool;
	else if (s=="d" || s=="f" || s=="r" || s=="double" || s=="float" || s=="real")
	  return type_real;
	else ValueError(("Incorrect value type specifier "+s).c_str());
      }
      else if (fld[0]=="list" || fld[0]=="l" || fld[0]=="array"){
	if (fld[1]=="s" || fld[1]=="str")
	  return type_string + type_array;
	else if (fld[1]=="i" || fld[1]=="int")
	  return type_int + type_array;
	else if (fld[1]=="bool" || fld[1]=="b")
	  return type_bool + type_array;
	else if (fld[1]=="d" || fld[1]=="f" || fld[1]=="r"
		 || fld[1]=="double" || fld[1]=="float" || fld[1]=="real")
	  return type_real + type_array;
	else ValueError(("Incorrect value type specifier "+s).c_str());
      }
      else ValueError(("Incorrect value type specifier "+s).c_str());
      return 0;
    }      

    void create_table(const std::vector<STRING_EX> & fnames, const std::vector<int> & ftypes,
		      const std::vector<unsigned char> _add){
      if (fnames.size()!=ftypes.size())
	IndexError("number of parameters does not match number o their type specifiers!");
      field_names = fnames;
      field_types = ftypes;
      _additive.clear();
      for (int i=0; i< ftypes.size(); i++){
	auto nc = &newcolumn(field_types[i]);
	//std::cout << "index = " << nc->index() << "\n";
	_xfields.push_back(*nc);
	if (_add.size()>i)
	  _additive.push_back(_add[i]);
	else
	  _additive.push_back(false);
      }
    }
    
    void create_table(const std::vector<STRING_EX> & fnames, const std::vector<STRING_EX> & ftypes){
      if (fnames.size()!=ftypes.size())
	IndexError("number of parameters does not match number o their type specifiers!");
      field_names = fnames;
      field_types.clear();
      std::vector<int> ftint;
      for (int i=0; i< ftypes.size(); i++)
	ftint.push_back(xftype(ftypes[i]));
      create_table(fnames,ftint,{});
      /*      for (int i=0; i< ftypes.size(); i++){
	auto nc = &newcolumn(field_types[i]);
	//std::cout << "index = " << nc->index() << "\n";
	_xfields.push_back(*nc);
	_additive.push_back(false);
	} */     
    }
    
    void clear_table(){
      _xfields.clear();
      _additive.clear();
    }
    
    void newline(){
      for (int i=0; i< field_types.size(); i++){
	int t = field_types[i];
	if (t&type_array){
	  if (t&type_real)
	    std::get<std::vector<std::vector<REAL>>>(_xfields[i]).push_back(std::vector<REAL>());
	  else if (t&type_int)
	    std::get<std::vector<std::vector<int>>>(_xfields[i]).push_back(std::vector<int>());
	  else if (t&type_bool)
	    std::get<std::vector<Boolarray>>(_xfields[i]).push_back(Boolarray());
	  else
	    std::get<std::vector<std::vector<STRING_EX>>>(_xfields[i]).push_back(std::vector<STRING_EX>());
	}
	else{
	  if (t&type_real)
	    std::get<std::vector<REAL>>(_xfields[i]).push_back(0e0);
	  else if (t==type_int)
	    std::get<std::vector<int>>(_xfields[i]).push_back(0);
	  else if (t==type_bool)
	    std::get<std::vector<unsigned char>>(_xfields[i]).push_back(false);
	  else
	    std::get<std::vector<STRING_EX>>(_xfields[i]).push_back("");
	}
      }
    }

    void insertline(int at){
      //std::cout << "insertline " << at << "\n";
      for (int i=0; i< field_types.size(); i++){
	int t = field_types[i];
	if (t&type_array){
	  if (t&type_real)
	    vecinsert(std::get<std::vector<std::vector<REAL>>>(_xfields[i]),at, std::vector<REAL>());
	  else if (t&type_int)
	    vecinsert(std::get<std::vector<std::vector<int>>>(_xfields[i]),at,std::vector<int>());
	  else if (t&type_bool)
	    vecinsert(std::get<std::vector<Boolarray>>(_xfields[i]),at,Boolarray());
	  else
	    vecinsert(std::get<std::vector<std::vector<STRING_EX>>>(_xfields[i]),at,std::vector<STRING_EX>());
	}
	else{
	  if (t&type_real)
	    vecinsert(std::get<std::vector<REAL>>(_xfields[i]),at,REAL(0));
	  else if (t==type_int)
	    vecinsert(std::get<std::vector<int>>(_xfields[i]),at,0);
	  else if (t==type_bool)
	    vecinsert(std::get<std::vector<unsigned char>>(_xfields[i]),at,(unsigned char)(false));
	  else
	    vecinsert(std::get<std::vector<STRING_EX>>(_xfields[i]),at,STRING_EX(""));
	}
      }

    }

    void eraseline(int at){
      for (int i=0; i< field_types.size(); i++){
	int t = field_types[i];
	if (t&type_array){
	  if (t&type_real)
	    vecerase(std::get<std::vector<std::vector<REAL>>>(_xfields[i]),at);
	  else if (t&type_int)
	    vecerase(std::get<std::vector<std::vector<int>>>(_xfields[i]),at);
	  else if (t&type_bool)
	    vecerase(std::get<std::vector<Boolarray>>(_xfields[i]),at);
	  else
	    vecerase(std::get<std::vector<std::vector<STRING_EX>>>(_xfields[i]),at);
	}
	else{
	  if (t&type_real)
	    vecerase(std::get<std::vector<REAL>>(_xfields[i]),at);
	  else if (t==type_int)
	    vecerase(std::get<std::vector<int>>(_xfields[i]),at);
	  else if (t==type_bool)
	    vecerase(std::get<std::vector<unsigned char>>(_xfields[i]),at);
	  else
	    vecerase(std::get<std::vector<STRING_EX>>(_xfields[i]),at);
	}
      }

    }

    int findex(const STRING_EX& f) const{
      for (int i=0; i< field_names.size(); i++)
	if (f==field_names[i])
	  return i;
      return -1;
    }
    
    void clone(xgeometry<REAL> &dst, const bool copy_data = true) {
      dst = xgeometry<REAL>(*this);
    }
    
    void get_format (std::vector<STRING_EX> & fn, std::vector<int> & ft) const
    {
      fn = field_names;
      ft = field_types;
    }

    void set_format(const std::vector<STRING_EX> & fn,
		    const std::vector<int> & ft) {
      clear_table();
      create_table(fn,ft,{});
      for (int i=0; i< nat(); i++)
	insertline(0);      
    }

    int field_type(int f) const{return field_types[f];}

    STRING_EX field_name(int f) const{return field_names[f];}

    template<class T>
    T xfield(int f, int i) const{
      try{
	return (T)(std::get<std::vector<typename XFTYPE<T>::xftype>>(_xfields[f])[i]);
      }
      catch (const std::bad_variant_access& e){
	TypeError(STRING_EX("Bad type for extra field in xeometry")+e.what());
	return T();}
    }
      
    template<typename T>
    T& xfield(int f, int i)    {
      //XFTYPE<T>::xftype val;
      try{
	return reinterpret_cast<T&>(std::get<std::vector<typename XFTYPE<T>::xftype> >(_xfields[f])[i]);
      }
      catch (const std::bad_variant_access& e){
	TypeError(STRING_EX("Bad type for extra feild in xgeometry")+e.what());
	return *(T*)nullptr;}
    }

    REAL charge(int j) const{
      int f = findex("charge");
      if (f==-1){
	KeyError("This xgeom does not have charge xfield");
	return 0e0;
      }
      else
	return xfield<REAL>(f,j);
    }
    
    REAL & charge(int j){
      int f = findex("charge");
      if (f==-1){
	KeyError("This xgeom does not have charge xfield");
	return *(REAL*)nullptr;
      }     
      else
	return xfield<REAL>(f,j);
    }
    
    template<typename T>
    int setxfield(int at, int f, const T& x)
    {
      //std::cout << at << " " << f << typeid(T).name() << "\n";
      try{
	std::get<std::vector<typename XFTYPE<T>::xftype> >(_xfields[f])[at]  = typename XFTYPE<T>::xftype(x);
      }
      catch (const std::bad_variant_access& e){
	return -1;
      }
      return 0;
    }

    
    
    template <typename... Types>
    void fillxfields(int i,  Types... ffields){
      //std::cout << i << ","<< f << " " << x << " " << typeid(x).name() << " " << "idx="<<_xfields[f].index() <<"\n";
      int f =0;
      auto x = {0,(setxfield(i,f++,ffields))... };
    }
    virtual void get_fields (int j, std::vector<fieldtypes> & v) const override{

      geometry<REAL>::get_fields(j,v);
      for (int f=0; f< field_types.size(); f++){
	int t = field_types[f];
	if (t&type_array){
	  if (t&type_real)
	    v.push_back(xfield<std::vector<REAL>>(f,j));
	  else if (t&type_int)
	    v.push_back(xfield<std::vector<int>>(f,j));
	  else if (t&type_bool)
	    v.push_back(xfield<Boolarray>(f,j));
	  else
	    v.push_back(xfield<std::vector<STRING_EX>>(f,j));
	}
	else{
	  if (t&type_real)
	    v.push_back(xfield<REAL>(f,j));
	  else if (t==type_int)
	    v.push_back(xfield<int>(f,j));
	  else if (t==type_bool)
	    v.push_back(xfield<bool>(f,j));
	  else
	    v.push_back(xfield<STRING_EX>(f,j));
	}
      }
    }
    
    virtual void set_fields (int j, const std::vector<fieldtypes> & v) override{
      std::vector<fieldtypes> v1 = {&v[0],&v[4]};
      std::vector<fieldtypes> v2(v.begin()+4, v.end());
      //std::cout << "xgeom::set_fields " << v1.size() << "\n";
      geometry<REAL>::set_fields(j,v1);
      if (v2.size() > field_types.size())
	IndexError("Number of provided fields greater than number of fields in xgeometry!!");
      for (int f=0; f< v2.size(); f++){
	int t = field_types[f];
	if (t&type_array){
	  if (t&type_real)
	    xfield<std::vector<REAL>>(f,j) = std::get<std::vector<REAL>>(v2[f]);
	  else if (t&type_int)
	    xfield<std::vector<int>>(f,j) = std::get<std::vector<int>>(v2[f]);
	  else if (t&type_bool)
	    xfield<Boolarray>(f,j) = std::get<Boolarray>(v2[f]);
	  else
	    xfield<std::vector<STRING_EX>>(f,j) = std::get<std::vector<STRING_EX>>(v2[f]);
	}
	else{
	  if (t&type_real)
	    xfield<REAL>(f,j) = std::get<REAL>(v2[f]);
	  else if (t==type_int)
	    xfield<int>(f,j) = std::get<int>(v2[f]);
	  else if (t==type_bool)
	    xfield<bool>(f,j) = std::get<bool>(v2[f]);
	  else
	    xfield<STRING_EX>(f,j) = std::get<STRING_EX>(v2[f]);
	}
      }
    }
    
    virtual void add (const STRING_EX & a, const vector3<REAL> & r1){
      geometry<REAL>::add(a,r1);
      newline();
    }

    virtual void add (const STRING_EX &a, const REAL _x, const REAL _y, const REAL _z){
      geometry<REAL>::add(a,_x,_y,_z);
      newline();
    }
    
    template <typename... Types>
    void add(const STRING_EX &a, const REAL _x, const REAL _y, const REAL _z, Types... ffields)
    {
      geometry<REAL>::add(a,_x,_y,_z);
      newline();
      fillxfields(nat()-1, ffields...);
    }
    
    virtual void add_fields( const std::vector<fieldtypes> & v){
      // all fields in one vector v including atom and x,y,z
      STRING_EX at = std::get<STRING_EX>(v[0]);
      REAL x = std::get<REAL>(v[1]);
      REAL y = std::get<REAL>(v[2]);
      REAL z = std::get<REAL>(v[3]);
      std::vector<fieldtypes> xv(v.begin()+4,v.end());
      add(at,x,y,z);
      set_fields(nat()-1,xv);
    }


    virtual void insert (int at, const STRING_EX & a, const vector3<REAL> & r1){
      //std::cout << "insert(a,r) " << at << " " << a << " " << r1 << "\n";
      geometry<REAL>::_insert(at,a,r1);      
      insertline(at);
    }

    
    template <typename... Types>
    void insert (int at, const STRING_EX & a,  REAL _x, REAL _y, REAL _z,Types... ffields){
      //std::cout << "insert: " << at << a << " " << _x << " " << _y << " " << _z<< "\n";
      geometry<REAL>::_insert(at,a,{_x,_y,_z});
      insertline(at);
      fillxfields(at, ffields...);
    }
    
    virtual void erase (int at)
    {
      //std::cout << "xgeom erase " << at << "\n";
      geometry<REAL>::_erase(at);
      eraseline(at);
    }

    virtual void reorder (const std::vector<int> & ord) override{
      
      for (int i=0; i<observers.size(); i++)
	observers[i]->reordered(ord, before);
      // fixme - might be inefficient for large molecules
      
      std::vector<STRING_EX> __atm(_atm);
      std::vector<vector3<REAL> > __crd(_crd);
      std::vector<Bool> __shadow(_shadow);
      
      std::vector<std::vector<fieldtypes>>  v(nat());
      
      for (int i=0; i<size(); i++) {

	get_fields(ord[i],v[i]);
	_shadow[i] = __shadow[ord[i]];	
      }
      for (int i=0; i<size(); i++)
	set_fields(i,v[i]);
      
	if (_typetable)
	  _typetable->reorder_types(ord);

        for (int i=0; i<observers.size(); i++)
          observers[i]->reordered(ord, after);

      }
    
    virtual bool is_xgeometry () const
    {return true;}
    
    int nfields()const{return _xfields.size();}
    
    //bool additive(int f) const(){return _additive[f];}

    bool& additive(int f){return *(bool*)(&_additive[f]);}

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    
    py::dict pyheader(){
      py::dict hd;
      for (int i=0; i<field_names.size(); i++){ 
	STRING_EX fn = field_names[i];
	STRING_EX ft = type_data::type_name[field_types[i]];
	hd[fn.c_str()] = py::cast(ft.c_str());
      }
      return hd;
    }

    void setfromlist(int j, py::list l){
      //py::print("setfromlist",l);
      if (j<0)
	j+= nat();
      if (j<0 || j>= nat())
	IndexError("atom number out of range");
      for (int i=0; i< nfields(); i++){
	int t = field_types[i];
	if (t&type_array){
	  if (t&type_real)
	    xfield<std::vector<REAL>>(i,j) = py::cast<std::vector<REAL>>(l[i]);
	  else if (t&type_int){
	    py::print(l[i],py::cast<std::vector<int>>(l[i]), j );
	    xfield<std::vector<int> >(i,j) = py::cast<std::vector<int>>(l[i]);
	  }
	  else if (t&type_bool)
	    xfield<Boolarray>(i,j) = Boolarray(py::cast<py::list>(l[i]));
	  else
	    xfield<std::vector<STRING_EX>>(i,j) = py::cast<std::vector<STRING_EX>>(l[i]);
	}
	else{
	  if (t&type_real)
	    xfield<REAL>(i,j) = py::cast<py::float_>(l[i]);
	  else if (t==type_int)
	    xfield<int>(i,j) = py::cast<int>(l[i]);
	  else if (t==type_bool)
	    xfield<bool>(i,j) = py::cast<bool>(l[i]);
	  else
	    xfield<STRING_EX>(i,j) =  py::cast<py::str>(l[i]);
	}
      }
    }
    
    void py_add(const py::list &l){
      STRING_EX at = py::cast<py::str>(l.attr("pop")(0));
      REAL x = py::cast<py::float_>(   l.attr("pop")(0));
      REAL y = py::cast<py::float_>(   l.attr("pop")(0));
      REAL z = py::cast<py::float_>(   l.attr("pop")(0));
      geometry<REAL>::add(at,x,y,z);
      newline();
      setfromlist(-1,l);
      //std::cout << "nat= " << nat()<< " nf=" << nfields() << "\n";
      
    }

    void py_insert(int i, const py::list &l){
      STRING_EX at = py::cast<py::str>(l.attr("pop")(0));
      REAL x = py::cast<py::float_>(   l.attr("pop")(0));
      REAL y = py::cast<py::float_>(   l.attr("pop")(0));
      REAL z = py::cast<py::float_>(   l.attr("pop")(0));
      geometry<REAL>::_insert(i,at,{x,y,z});
      insertline(i);
      setfromlist(i,l);
    }
    
    virtual py::list py_getitem(int j) const {
      if (j<0)
	j+= nat();
      if (j<0 || j>= nat())
	IndexError("atom number out of range");
      py::list res;
      res.append(atom(j));
      res.append(coord(j)[0]);
      res.append(coord(j)[1]);
      res.append(coord(j)[2]);
      for (int i=0; i< nfields(); i++){
	int t = field_types[i];
	if (t&type_array){
	  if (t&type_real)
	    res.append(py::cast(xfield<std::vector<REAL>>(i,j)));
	  else if (t&type_int)
	    res.append(py::cast(xfield<std::vector<int> >(i,j)));
	  else if (t&type_bool)
	    res.append(py::cast(xfield<Boolarray>(i,j)));
	  else
	    res.append(py::cast(xfield<std::vector<STRING_EX>>(i,j)));
	}
	else{
	  if (t&type_real)
	    res.append(py::cast(xfield<REAL>(i,j)));
	  else if (t==type_int)
	    res.append(py::cast(xfield<int>(i,j)));
	  else if (t==type_bool)
	    res.append(py::cast(xfield<bool>(i,j)));
	  else
	    res.append(py::cast(xfield<STRING_EX>(i,j)));
	}
      }
      return res;
    }
    
    virtual void py_setitem(int j, py::list  l) {
      STRING_EX at = py::cast<py::str>(l.attr("pop")(0));
      REAL x = py::cast<py::float_>(   l.attr("pop")(0));
      REAL y = py::cast<py::float_>(   l.attr("pop")(0));
      REAL z = py::cast<py::float_>(   l.attr("pop")(0));
      _atm[j]  = at;
      _crd[j][0] = x;
      _crd[j][1] = y;
      _crd[j][2] = z;
      setfromlist(j,l);
    }
    
    py::object py_getfield1(int i) {
      TypeError("xgeometry::field accepts 2 indicies");
      return py::none();
    }
    
    void py_setfield1(int i, const py::object & o) {
      TypeError("xgeometry::field accepts 2 indicies");
    }

    py::object py_getfield(int i, int j)
    {
      /*
      std::cout << " nf=" << nfields()<< "\n";
      std::cout << " nat=" << nat()<< "\n";
      std::cout << "xfield"<< " " << i << " " << j <<  "\n";
      */
      if (j<0)
	j += nat();
      std::cout << i<< " "<< j << " " << nat() << " " << nfields() << "\n";
      if (j<0|| j>=nat())
	IndexError("out or range atom number");
      if (i<0)
	i+= nfields();
      if (i<0 || i>=field_types.size() )
	IndexError("xfield index out of range");
      int t = field_types[i];
      //std::cout << "type" << t << " nat=" << nat()<< "\n""\n";
     
      if (t&type_array){
	if (t&type_real)
	  return py::cast(xfield<std::vector<REAL> >(i,j));
	else if (t&type_int)
	  return py::cast(xfield<std::vector<int>>(i,j));
	else if (t&type_bool)
	  return py::cast(xfield<Boolarray>(i,j));
	else
	  return py::cast(xfield<std::vector<STRING_EX>>(i,j));
	}
      else{
	if (t&type_real)
	  return py::cast(xfield<REAL>(i,j));
	else if (t==type_int)
	  return py::cast(xfield<int>(i,j));
	else if (t==type_bool)
	  return py::cast(xfield<bool>(i,j));
	else
	  return py::cast(xfield<STRING_EX>(i,j));
      }
    }
		       
    void py_setfield(int i, int j, const py::object & o)
    {
            if (j<0)
	j += nat();
      if (j<0|| j>=nat())
	IndexError("out or range atom number");
      if (i<0)
	i+= nfields();
      if (i<0 || i>=field_types.size() )
	IndexError("xfield index out of range");
      int t = field_types[i];
      if (t&type_array){
	if (t&type_real)
	  xfield<std::vector<REAL> >(i,j) = py::cast<std::vector<REAL>>(o);
	else if (t&type_int)
	  xfield<std::vector<int>>(i,j)= py::cast<std::vector<int>>(o);
	else if (t&type_bool)
	  xfield<Boolarray>(i,j) = Boolarray(py::cast<py::list>(o));
	else
	  xfield<std::vector<STRING_EX>>(i,j)= py::cast<std::vector<STRING_EX>>(o);
      }
      else{
	if (t&type_real)
	  xfield<REAL>(i,j) = py::cast<REAL>(o);
	else if (t==type_int)
	  xfield<int>(i,j) = py::cast<int>(o);
	else if (t==type_bool)
	  xfield<bool>(i,j) = py::cast<bool>(o);
	else
	  xfield<STRING_EX>(i,j) = py::cast<STRING_EX>(o);
	}
    }
    
    py_2indexed_property<SELF,py::object, py::object, int,
			 &SELF::py_getfield1, &SELF::py_setfield1,
			 &SELF::py_getfield, &SELF::py_setfield > py_xfield;
    
    
    bool py_getadd (int f) {
        if (f<0) f+=nfields();
        if (f<0 || f>=nfields()) IndexError("geometry::atom::field out of range");
        return additive(f);
      }

    void py_setadd (int f, const bool & a) {
        if (f<0) f+=nfields();
        if (f<0 || f>=nfields()) IndexError("geometry::atom::field out of range");
	additive(f) = a;
       }

    py_indexed_property< SELF, bool, int,
			 &SELF::py_getadd, &SELF::py_setadd> py_additive;

#endif
    
    xgeometry(CELL & __cell, const py::dict& header ):
      geometry<REAL>(__cell)
    {      
      std::vector<STRING_EX> fnames;
      std::vector<STRING_EX> ftypes;
      for (auto it : header){
	fnames.push_back(py::cast<STRING_EX>(it.first));
	ftypes.push_back(py::cast<STRING_EX>(it.second));
      }
      create_table(fnames,ftypes);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif
    }
    
    xgeometry(int dim, const py::dict& header):geometry<REAL>(dim)
    {
      std::vector<STRING_EX> fnames;
      std::vector<STRING_EX> ftypes;
      for (auto it : header){
	fnames.push_back(py::cast<STRING_EX>(it.first));
	ftypes.push_back(py::cast<STRING_EX>(it.second));
      }
      create_table(fnames,ftypes);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif
    }
    
    xgeometry(const std::vector<STRING_EX> & fnames, const std::vector<STRING_EX> & ftypes,
	      CELL & __cell = periodic_cell<REAL>(0)): geometry<REAL>(__cell){
      create_table(fnames,ftypes);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif

    }
    
    xgeometry(const std::vector<STRING_EX> & fnames, const std::vector<STRING_EX> & ftypes,
	      int dim=0): geometry<REAL>(dim){
      create_table(fnames,ftypes);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif

    }

    xgeometry(int dim, const STRING_EX & _name = ""): geometry<REAL>(dim,_name){
      create_table({},{});
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif
    }

    xgeometry(const xgeometry<REAL> & src): geometry<REAL>(src),
					    field_names(src.field_names), field_types(src.field_types),
					    _additive(src._additive), _xfields(src._xfields){
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
      py_xfield.bind(this);
      py_additive.bind(this);
#endif
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    void py_erase1(int at){
      erase(at);
    }

    void py_erase2(const std::vector<int> &ats){
      erase(ats);
    }
    static void py_export(py::module m, const char * pyname) {
      STRING_EX  xfldname = fmt::format("{0}_{1}",pyname,"fld_indexed_property");
      STRING_EX  xaddname = fmt::format("{0}_{1}",pyname,"add_indexed_property");
      py_2indexed_property<SELF,py::object, py::object, int,  &SELF::py_getfield1, &SELF::py_setfield1,
			   &SELF::py_getfield, &SELF::py_setfield >::py_2export(m, xfldname.c_str());
      py_indexed_property< SELF, bool, int,
			   &SELF::py_getadd, &SELF::py_setadd>::py_export(m, xaddname.c_str());
      py::class_<xgeometry<REAL>, geometry<REAL>, std::shared_ptr<xgeometry<REAL>>>(m,pyname)
	.def(py::init<int,const py::dict& >(), py::arg("dim")=0, py::arg("header"))
	.def(py::init<CELL&,const py::dict& >(),
	     py::arg("_cell")=periodic_cell<REAL>(0),
	     py::arg("header"))
	.def("header", &xgeometry<REAL>::pyheader )
	.def("add",&SELF::py_add)
	.def("insert",&SELF::py_insert)
	//	.def("erase",&SELF::py_erase1)
	//.def("erase",&SELF::py_erase2)
	.def_readwrite("field", &SELF::py_xfield)
	.def_readwrite("additive", &SELF::py_additive)
	.def("__getitem__", &SELF::py_getitem)
	.def("__setitem__", &SELF::py_setitem)
	.def("nfields", &SELF::nfields)
	.def_readonly("fnames", &SELF::field_names)
	.def_readonly("ftypes", &SELF::field_types)
	.def("newline", &SELF::newline)
	//.def_readwrite("additive", &SELF::_additive)
	;
    }
    
#endif
  };

};
#endif
