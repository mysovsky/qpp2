#ifndef QPP_INDEX_H
#define QPP_INDEX_H

#//if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
//#endif


#include <ostream>
#include <sstream>
#include <signal.h>
#include <initializer_list>
#include <stdexcept>
#include <mathf/lace3d.hpp>
#include <data/errors.hpp>
#include <io/strfun.hpp>

namespace qpp {
  
  class index;
  
  index atom_index(int at, const index & I);

  // -------------------------------------------------------------

  bool compare_atindex(const index & at1,
                       const index & at2);


  /// \brief index is a general purpose complex index, having DIM integer components
  class index {

    std::shared_ptr<int> idx;

    public:

    // \brief DIM(int) The index dimension - number of components
    int DIM;

      //! Typecast operator to int type returns the 0th component of the index
    inline operator int () const {return idx.get()[0];}

      //! The d-th component of the index
      inline int operator () (int d) const {
	if (d<0)
	  d+=DIM;
	if (d<0 || d>=DIM)
	  IndexError("qpp::index: out of range");
	return idx.get()[d];
      }

      //! The d-th component of the index
      inline int& operator () (int d) {
	if (d<0)
	  d+=DIM;
	if (d<0 || d>=DIM)
	  IndexError("qpp::index: out of range");
	return idx.get()[d];
      }

      //\brief Typecast from integer: the 0th component of index is set to i, the rest is set to 0
      inline index& operator= (int i) {

        idx.get()[0] = i;
        for (int d=1; d<DIM; d++)
          idx.get()[d] = 0;
        return *this;

      }

      //\brief Assigment operator
      inline index& operator= (const index & I) {

        if (DIM != I.DIM){
            DIM = I.DIM;
            idx = std::shared_ptr<int>(new int[DIM], [](int *p) { delete [] p;});
          }

        for (int d=0; d<DIM; d++)
          idx.get()[d] = I(d);
        return *this;

      }

      //\brief The same as assigment operator in the form of explicitly called method
      inline void set (const index & I) {
        *this = I;
      }

      //! Componentwise addition of two indicies
      inline index operator+(const index & I) const {
	if (DIM!=I.DIM)
	  IndexError("DIfferent dimensions in index addition");

        index res = D(DIM);
        for (int d=0; d<DIM; d++)
          res(d) = idx.get()[d] + I(d);
        return res;

      }

      //\brief Componentwise subtraction of two indicies
      inline index operator-(const index & I) const {
	if (DIM!=I.DIM)
	  IndexError("DIfferent dimensions in index subtraction");

        index res = D(DIM);
        for (int d=0; d<DIM; d++)
          res(d) = idx.get()[d] - I(d);
        return res;

      }

      //\brief Index I is added to this index componentwise
      inline index& operator+= (const index & I) {
	if (DIM!=I.DIM)
	  IndexError("DIfferent dimensions in index addition");
        for (int d=0; d<DIM; d++)
          idx.get()[d] = idx.get()[d] + I(d);
        return *this;

      }

      //\brief Index I is subtracted from this index componentwise
      inline index& operator-= (const index & I){
	if (DIM!=I.DIM)
	  IndexError("DIfferent dimensions in index subtraction");
        for (int d=0; d<DIM; d++)
          idx.get()[d] = idx.get()[d] - I(d);
        return *this;

      }

      /*! \brief Using std::initializer_list to set the components of this index. Example:
      qpp::index I({1,2,3,4,5});
      std::cout << I << std::endl; // (1,2,3,4,5)
      I.set({5,4,3,2,1});
      std::cout << I << std::endl; // (5,4,3,2,1)
      I.set({1,2,3}); // IndexError: Wrong number of index components
     */
      inline void set (const std::initializer_list<int> &li) {

        if (li.size() != DIM) IndexError("Wrong number of index components");
        int d=0;
        for (int i : li)
          idx.get()[d++] = i;
      }

      
    index(int dim)
    {
      DIM = dim;
      idx = std::shared_ptr<int>(new int[DIM], [](int *p) { delete [] p;});
    }
    

    index () {
      DIM = 0;
      idx = nullptr;
    }

      /*
    index(int i)
    {
      DIM=1;
      idx = new int;
      del = true;
      *idx = i;
    }
    */

    index (const index & I) {
      
      DIM = I.DIM;
      idx = std::shared_ptr<int>(new int[DIM], [](int *p){delete [] p;});
      for (int d=0; d<DIM; d++)
	idx.get()[d] = I(d);
    }

    index (const index & I, int d1, int d2 = -1) {
      //std::cout << "d1d2 case " << I << " " << I.DIM << " "<< d1 << " " << d2 << "\n";
      if (d1<0) d1+=I.DIM;
      if (d2<0) d2+=I.DIM;

      //std::cout << "d1d2 case " << d1 << " " << d2 << "\n";

      if (d1<0 || d1>=I.DIM) //raise(SIGSEGV);
	IndexError("d1 d2 out of range");
      if (d2<0 || d2>=I.DIM) //raise(SIGSEGV);
	IndexError("d1 d2 out of range");
      DIM = d2-d1+1;
      idx = std::shared_ptr<int>(new int[DIM], [](int *p) { delete [] p;});     
      for (int i=0; i< DIM; i++)
	idx.get()[i] = I(i+d1);
    }

    index (const std::initializer_list<int> &li) {
      DIM = li.size();
      idx = std::shared_ptr<int>(new int[DIM], [](int *p){delete [] p;});
      set(li);      
    }
    
    ~index () {
      
      idx = nullptr;
      
    }
    // \brief subindex, i.e. I[d1:d2] part of thid index
    inline index cut (int d1, int d2 = -1) const {
      return index(*this,d1,d2);
    }
    
    inline bool operator== (const index &I) const {
      
      bool res = DIM == I.DIM;
      if (res)
	for (int d=0; d<DIM; d++)
	  if (idx.get()[d]!=I(d)){
	    res = false;
	    break;
	  }
      return res;
    }

    bool is_zero() const {
      for (int d=0; d<DIM; d++)
	if ((*this)(d)!=0)
	  return false;
      return true;
    }

    inline bool operator!= (const index &I) const {
      return ! ((*this) == I);
    }

    inline index head(int d)const{
      return cut(0,d);
    }
    
    inline index tail(int d) const{
      if (DIM<=d)
	return index();
      else
	return cut(d);
    }

    inline index cat(const index & I)const{
      index res(DIM+I.DIM);
      for (int d = 0; d< DIM; d++)
	res(d) = (*this)(d);
      for (int d = 0; d< I.DIM; d++)
	res(DIM+d) = I(d);
      return res;
    }
    
    struct factory {
      
      int DIM;
      
      factory (int dim) {DIM=dim;}
      
      index all (int a) {
	
	index t(DIM);
	
	for (int d=0; d<DIM; d++)
	  t(d) = a;
	return t;
	
      }
      /*  
      index atom (int a) {
	
	index t = all(0);
	t(0)=a;
	return t;
	
      }
      */
      operator index () {return all(0);}
      
    };
    
    static factory D(int dim) {
      
      return factory(dim);
    }
    
    //#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      // --------------- PYTHON -------------------------------
    index(py::args args){
      DIM = py::len(args);
      idx = std::shared_ptr<int>(new int[DIM], [](int *p) { delete [] p;});
      for (int d=0; d < DIM; d++)
	idx.get()[d] = py::cast<int>(args[d]);
    }

      int py_getitem(int d) const{
	if (d<0)
	  d+=DIM;
	if (d<0 || d>=DIM)
	  IndexError("qpp::index: out of range");
	return idx.get()[d];
      }

      void py_setitem(int d, int v){
	if (d<0)
	  d+=DIM;
	if (d<0 || d>=DIM)
	  IndexError("qpp::index: out of range");
        idx.get()[d] = v;
      }

      index(const py::list &l){

        // Assuming that type checks have been performed inside python code
        DIM = py::len(l);
        idx =  std::shared_ptr<int>(new int[DIM], [](int *p){delete [] p;});
        for (int d=0; d<DIM; d++)
          idx.get()[d] =py::cast<int>(l[d]);
      }

      index(const py::tuple &l){

        // Assuming that type checks have been performed inside python code
        DIM = py::len(l);
        idx = std::shared_ptr<int>(new int[DIM], [](int *p){delete [] p;});
        for (int d=0; d<DIM; d++)
          idx.get()[d] = py::cast<int>(l[d]);
      }

      std::string print(){
        std::string _tmp = "{";

        _tmp = "idx(";// + t2s(idx.get()[0]);

        for (int d=0; d<DIM; d++){
          _tmp += t2s(idx.get()[d]);
	  if (d<DIM-1)
	    _tmp += ",";
	}
        _tmp += ")";

        return _tmp;

      }

      index py_add(const index &I2) const
      { return (*this)+I2; }

      index py_sub(const index &I2) const
      { return (*this)-I2; }

    static void py_export( py::module m, const char * pyname){

      py::class_<index >(m, pyname)
	.def(py::init<>())
	.def(py::init<py::list&>())
	.def(py::init<py::tuple&>())
	.def(py::init<index const&>())
	//TODO: last parameter is optional
	.def(py::init<index const&, int, int >())
	.def(py::init<py::args>())
	//.def(py::init<int>())
	.def("__getitem__",&index::py_getitem)
	.def("__setitem__",&index::py_setitem)
	.def("cut",  &index::cut, py::arg("d1"), py::arg("d2")=-1)
	.def("head",  &index::head)
	.def("tail", &index::tail)
	.def("cat", &index::cat)
	//	.def(py::str(py::self))
	//	.def(py::repr(py::self))
	.def(py::self + py::self)
	.def(py::self - py::self)
	//.def("__add__", & index::py_add )
	//.def("__sub__", & index::py_sub )
	.def(py::self == py::self)
	.def(py::self!= py::self)
	.def("__str__", &index::print)
	.def("__repr__", &index::print)
	.def_readonly("DIM", &index::DIM)
	.def("__rand__", [](const index & self, int at)->index{ return atom_index(at,self);})
	;
    }

    //#endif

  };

  // ------------------------------------------------------

  template<typename _CharT, class _Traits>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT,
             _Traits>& __os,
             const index &I){

    std::basic_ostringstream<_CharT, _Traits> __s;
    __s.flags(__os.flags());
    __s.imbue(__os.getloc());
    __s.precision(__os.precision());

    if (I.DIM > 0){

        __s  << "(" << I(0);
        for (int d=1; d<I.DIM; d++)
          __s << "," << I(d);
        __s << ")";
      }
    else

      __s << "(0)";
    return __os << __s.str();
  }

  // ------------------- iterator class --------------------
  // Iterator allows you run through all (or some) atoms of this cell
  // and of surrounding cells

  class iterator : public index {

    protected:
      index a, b;
      bool _end;
      // a - from
      // b - to

      //using index::idx;

      inline void inc(){

        if (DIM==0){

            _end = true;
            return;
          }

        int d = 0;

        while (++(*this)(d) > b(d)){

            (*this)(d)=a(d);
            if (++d >= DIM) {

                _end = true;
                break;
              }
          }
      }

    public:

    iterator(const index & _a, const index & _b) :
      index(_a), a(_a), b(_b){
      
      _end = false;
    }

    iterator(const py::tuple &_a,const py::tuple &_b):
      a(_a), b(_b){
      _end = false;
      set(a);
    }
    
     iterator(const py::list &_a,const py::list &_b):
      a(_a), b(_b){
       set(a);
      _end = false;
      }
    
    iterator(const std::initializer_list<int> &li1,
	     const std::initializer_list<int> &li2):
      a(li1),b(li2){
      set(a);
      //std::cout << "a=" << a << " b="<< b << *this << "\n";
      _end = false;
    }
    inline void reset(){
      
      _end = false;
      set(a);
    }
    
    iterator & operator++(int){
      
      inc();
      return *this;
    }
    
    inline bool end(){return _end;}
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    
    // --------------- PYTHON -------------------------------
    
    index py_next(){
      
      if (end())
	StopIter();
      index res = *this;
      inc();
      return res;
    }
    
      static void py_export(py::module m, const char * pyname){

        py::class_<iterator>(m, pyname)
	  .def(py::init<const index&, const index&>())
	  .def(py::init<py::list,py::list>())
	  .def(py::init<py::tuple,py::tuple>())
	  .def("next", & iterator::py_next)
	  .def("__next__", & iterator::py_next)
	  ;
      }
    
#endif

  };


  class index_range{

      index a,b;

    public:

    index_range(const std::initializer_list<int> &li1,
		const std::initializer_list<int> &li2):
      a(li1),b(li2){
      
    }

    index_range(const index & _a,
		const index & _b) :
      a(_a), b(_b) {}
    
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    index_range(const py::list & l1,
		const py::list & l2) :
      a(l1), b(l2){}
    
    index_range(const py::tuple & l1,
		const py::tuple & l2) :
      a(l1), b(l2) {}
#endif    
    iterator  __iter__()
    { return iterator(a,b); }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)    
    static void py_export(py::module m,
			  const char * pyname){
      py::class_<index_range>(m, pyname)
	.def(py::init<const index&, const index&>())
	.def(py::init<const py::list&, const py::list&>())
	.def(py::init<const py::tuple&, const py::tuple&>())
	.def("__iter__", & index_range::__iter__)
	;
    }
#endif
  };



}

#endif
