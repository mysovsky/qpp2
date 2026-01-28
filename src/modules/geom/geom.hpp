#ifndef QPP_GEOM_H
#define QPP_GEOM_H

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

#include <algorithm>
#include <vector>
#include <cmath>
#include <functional>
#include <variant>
#include <mathf/lace3d.hpp>   
#include <data/data.hpp>
#include <symm/index.hpp>
#include <symm/cell.hpp>

namespace qpp {

  const uint32_t GEOM_DEFAULT_RESERVE_AMOUNT = 128;

  template< class REAL = double>
  class geometry;

  //--------------------------------------------------------------//

  /*!

The supercell concept generalization for the geometry class looks like:
    template< int DIM, class REAL>
    class CELL{
    public:
    index<DIM> begin();
    index<DIM> end();
    vector3d<REAL> transform(const vector3d<REAL> & r, const index<DIM> & I);
    bool within(const vector3d<REAL> & r);
    vector3d<REAL> reduce(const vector3d<REAL> & r);
    // find high symmetry point within "radius" distance from given point "r"
    // makes sence for rotational symmetries
    vector3d<REAL> symmetrize(const vector3d<REAL> & r, REAL radius);
    // fractional to cartesian and vice versa transformation
    // makes sence only for periodic translational cells
    vector3d<REAL> frac2cart(const vector3d<REAL> & r);
    vector3d<REAL> cart2frac(const vector3d<REAL> & r);
  };

        Any class inheriting from periodic_cell<REAL and >having these methods - can be used as a "supercell" for some geometry

    template <int DIM>
    struct atom_index {
    int atom;
    index<DIM> cell;

    atom_index(int _atom, const index<DIM> & _cell = index<DIM>::all(0))
    {
      atom = _atom;
      cell = _cell;
    }

    inline operator int() const {return atom;}

    inline bool operator==(const atom_index<DIM> & ai)
    {
      return atom == ai.atom && cell == ai.cell;
    }

    inline bool operator!=(const atom_index<DIM> & ai)
    {
      return atom != ai.atom || cell != ai.cell;
    }

    };*/

  enum before_after {before = 0, after = 1};

  /*! \class geometry_observer
    \brief Geometry updated objects
    One geometry can maintain arbitrary number of "observers", i.e.
    objects which need to know about the changes made to geometry.
    Geometry will inform them all when atoms are added, inserted or removed
  */

  template <class REAL>
  struct geometry_observer {

      virtual void added (before_after, const STRING_EX &,
                          const vector3<REAL> &) =0;
      virtual void inserted (int at, before_after,
                             const STRING_EX &, const vector3<REAL> &) = 0;
      virtual void changed (int at, before_after,
                            const STRING_EX &, const vector3<REAL> &) = 0;
      virtual void erased (int at, before_after) =0;
      virtual void shaded (int at, before_after, bool) =0;
      virtual void reordered (const std::vector<int> &, before_after) = 0;

  };

  /*! \class geometry
    \brief geometry basically can store atomic symbols and coordinates for a molecule or crystal.
    However, the functionality of geometry class is much more than that.
    First, besides atomic symbols and coordinates any other additional data of real,
    integer, boolean or string type can be stored for each atom (see xgeometry).
    Second, geometry can handle any type of space symmetry, periodic or nonperiodic.
  */

  template <class REAL>
  struct geometry_type_table: public geometry_observer<REAL>{
  private:
    std::vector<STRING_EX> _atm_types;
    std::vector<int> _type_table;
    std::vector<REAL> _symm_rad;
    geometry<REAL> & geom;
    
  public:
    bool auto_update;
    bool auto_symmetrize;
    REAL default_symmetrize_radius;
    typedef geometry_type_table<REAL> SELF;

    geometry_type_table(geometry<REAL>& g, bool aupd = false):geom(g){      
      auto_update = aupd;
      auto_symmetrize = false;
      default_symmetrize_radius = 0e0;
    }

    //! Number of atomic type for atom named at
    inline int type_of_atom (const STRING_EX & at) const {

      for (int t=0; t < _atm_types.size(); t++)
	if ( _atm_types[t] == at )
	  return t;
      return -1;
    }
    //! Number of different atomic types in molecule
    inline int n_atom_types () const {
      return _atm_types.size();
    }

    //! Then synonim for n_atom_types()
    inline int n_types () const
    { return n_atom_types();}

    /*! \brief The string name of atom of type number t
      @param t - the atomic type number
    */
    inline STRING_EX atomic_type (int t) const {
      if (t<0)
	t+= _atm_types.size();
      if (t<0|| t>= _atm_types.size())
	IndexError("atomic type out of range");
      return _atm_types[t];
    }

    inline STRING_EX atom_name (int t) const {
      if (t<0)
	t+= _atm_types.size();
      if (t<0|| t>= _atm_types.size())
	IndexError("atomic type out of range");
      return _atm_types[_type_table[t]];
    }

    inline int type(int i){
      if (i<0)
	i+= geom.nat();
      if (i<0|| i>=geom.nat())
	IndexError("atom out of range");
      return _type_table[i];
    }

    inline int operator()(int i){
      if (i<0)
	i+= geom.nat();
      if (i<0|| i>=geom.nat())
	IndexError("atom out of range");
      return _type_table[i];
    }

    //! Synonim for type_of_atom(const STRING_EX & at)
    //inline int type (const STRING_EX & at) const
    // { return type_of_atom(at); }

    int define_type (const STRING_EX & at) {
      int t = type_of_atom(at);
      if (t==-1){
	t = _atm_types.size();
	_atm_types.push_back(at);
	_symm_rad.push_back(default_symmetrize_radius);
      }
      return t;

    }
    //inifficent - needs to be cached
    inline int atom_count_by_type (int atype) {      
      int count = 0;      
      for(int i = 0; i < _type_table.size(); i++)
	if(_type_table[i] == atype) count+=1;
      return count;      
    }
    
    void clear() {
      _atm_types.clear();
      _symm_rad.clear();
      _type_table.resize(geom.size());
    }
    // --------------- Symmetrization radius for each atomic type ----

    // whether to be aware of high symmetry points when placing atoms
    // (if atom is placed close to high symmetry point, it would be
    // also very close to its image)

    REAL symmetrize_radius (int t) const {
      if (t<0)
	t+= _atm_types.size();
      if (t<0|| t>= _atm_types.size())
	IndexError("atomic type out of range");
      return _symm_rad[t];
    }

    void set_symmetrize_radius (int t, REAL rad) {
      if (t<0)
	t+= _atm_types.size();
      if (t<0|| t>= _atm_types.size())
	IndexError("atomic type out of range");
      _symm_rad[t] = rad;
    }

    REAL symmetrize_radius_by_label (const STRING_EX & a) const {

      int t =type_of_atom(a);

      //std::cerr << "t= " << t << "\n";

      if (t==-1)
	return default_symmetrize_radius;
      else
	return symmetrize_radius (t);

    }

    void set_symmetrize_radius_by_label(const STRING_EX & a, REAL rad) {

      int t = type_of_atom(a);

      if (t == -1){
	t = _atm_types.size();
	_atm_types.push_back(a);
	_symm_rad.push_back(rad);
      }
      else
	_symm_rad[t] = rad;

    }

    virtual void build() {
      
      _atm_types.clear();
      _symm_rad.clear();
      _type_table.resize(geom.size());
      
      for (int i=0; i<geom.size(); i++) {
	
	int t = type_of_atom(geom.atom(i));
	
	if (t == -1) {
	  
	  t = _atm_types.size();
	  _atm_types.push_back(geom.atom(i));
	  _symm_rad.push_back(default_symmetrize_radius);
	  
	}
	
	_type_table[i] = t;
	
      }
      
    }
    
      void reorder_types(const std::vector<int> & ord) {

        if (_type_table.size() != geom.size()) return;

        std::vector<int> __type_table(_type_table);
        for (int i=0; i<geom.size(); i++)
          _type_table[i] = __type_table[ord[i]];

      }

    virtual void added (before_after ba, const STRING_EX &at,
			const vector3<REAL> &r)
    {
      //      std::cout << "added "<< at << " "<< ba << "\n";
      if (!auto_update || ba == before) return;
      _type_table.push_back(define_type(at));
      _symm_rad.push_back(default_symmetrize_radius);
    }
    
    virtual void inserted (int i, before_after ba,
			   const STRING_EX &at, const vector3<REAL> &){
      if (!auto_update || ba == before) return;
      _type_table.insert(_type_table.begin()+i,define_type(at));
      _symm_rad.insert(_symm_rad.begin()+i,default_symmetrize_radius);
    }
    
    virtual void changed (int at, before_after,
			  const STRING_EX &, const vector3<REAL> &){}
    virtual void erased (int at, before_after){}
    virtual void shaded (int at, before_after, bool){}
    virtual void reordered (const std::vector<int> &, before_after) {}
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    py::dict py_symmrad2dict () {
      
      py::dict d;
      d["default"] = default_symmetrize_radius;
      
      for (int i=0; i<n_atom_types(); i++)
	d[(atomic_type(i)).c_str()] = symmetrize_radius(i);
      
      return d;
      
    }
    
    static void py_export(py::module m, const char * pyname){
      //std::cout << "Exporting typetable " << pyname << "\n";
      py::class_<SELF,std::shared_ptr<SELF>>(m,pyname)
	.def(py::init<geometry<REAL>& ,bool>(), py::arg("geom"), py::arg("auto_update")=false)
	.def( "type_of_atom", &SELF::type_of_atom)
	.def("n_atom_types",  &SELF::n_atom_types)
	.def("atomic_type",   &SELF::atomic_type)
	.def("define_type",   &SELF::define_type)
	.def("atom_count_by_type",  &SELF:: atom_count_by_type)
	.def("symmetrize_radius", &SELF::symmetrize_radius)
	.def("set_symmetrize_radius", &SELF::set_symmetrize_radius)
	.def("set_symmetrize_radius", &SELF::set_symmetrize_radius_by_label)
	.def("symmetrize_radius", &SELF::symmetrize_radius_by_label)
	.def("build",  &SELF::build)
	.def("reorder_types",&SELF:: reorder_types)
	.def_readwrite("auto_update", &SELF::auto_update)
        .def_readwrite("auto_symmetrize",
		       &SELF::auto_symmetrize,
		       "bool auto_symmetrize: How to treat atoms added near \
                     'high symmetry points\'.\nIf auto_symmetrize==true,"
		       " put the new atom in high symmetry point instead of "
		       "given coordinates, when the distance\n between atom and "
		       "its image less then symmetrization radius (symmrad)")
	
	.def_readwrite("default_symmetrize_radius",
		       &SELF::default_symmetrize_radius,
		       "real default_symmrad: symmetrization "
		       "radius for all atoms")
	/*
	  .def_readwrite("symmrad",
	  &qpp::geometry<REAL>::py_symmrad,
	  "Usage: symmrad[at] with string at. "
	  "Symmetrization radius for atom at")
	*/
	.def("symmrad_to_dict",
	     &SELF::py_symmrad2dict,
	     "Outputs the symmrad values in dictionary form")
	;
      
    }
#endif
  };
  
  // ----------------------------------------------------------------------------
  // ----- The main class here - geometry --------------------------------------
  //----------------------------------------------------------------------------
  
  template< class REAL>
  class geometry {
    
  public:
    
    // Storage of atoms
    std::vector<STRING_EX> _atm;
    
    // Storage of coordinates
    std::vector<vector3<REAL> > _crd;
    
    //! Special logical array allows to "hide" or "shadow" some atoms
    std::vector<Bool> _shadow;

    typedef std::variant<REAL, int, bool, STRING_EX, std::vector<REAL>,
			 std::vector<int>, std::vector<STRING_EX>, Boolarray> fieldtypes;
    
  protected:
    // The dependent objects array
    std::vector<std::shared_ptr<geometry_observer<REAL>>> observers;
    bool has_observers;
    
  public:    
    typedef geometry_observer<REAL> DEP;
    typedef geometry<REAL> SELF;
    typedef periodic_cell<REAL> CELL;
    //! geometry tolerance - maximal distance between two atoms
    //! where they are considered at the same point
    // ! is the reference of cell->tol_transl
    std::shared_ptr<REAL> tol_geom;

      //! default value for tol_geom
      static REAL tol_geom_default;

      //! cell for periodic boundary conditions and on-the-fly transformations
    std::shared_ptr<CELL> cell;

    //! whether to update typetable on-the-fly
    //bool auto_update_types;
    
    //! whether fractional coordinates are used
    bool frac;
    
    //! the name of this geometry
    STRING_EX name;
    
    virtual bool is_xgeometry () const
    {return false;}
    
    //! \brief Dimension of periodicity
    inline int DIM() const
    {return cell->DIM;}

    void set_cell(CELL &cl){
      cell = std::shared_ptr<CELL>(&cl,[](CELL*p){});
    }
    
    //! \brief Number of atoms
    inline int size () const
    {return _crd.size();}

    //! \brief Number of atoms
    inline int nat () const
    {return _atm.size();}

    //! \brief The name of i-th atom
    inline STRING_EX & atom (int at)
    {return _atm[at];}

    inline STRING_EX atom (int at) const
    {return _atm[at];}

    /// \brief Gives the coordinates of an atom in the geometry
    ///  @param at - the number of atom in the geometry
    inline vector3<REAL> coord (int at) const {return _crd[at];}
    
    inline vector3<REAL>& coord (int at) {return _crd[at];}
    
    /// \brief real-space position of atom number at
    /// \param at - the number of the atom in the geometry
    inline vector3<REAL> r (int at) const{
      vector3<REAL> r1 = _crd[at];
      if (frac)
	r1 = cell->frac2cart(r1);  
      return cell->transform(r1, index::D(DIM()).all(0));
    }

    /// \brief real-space position of an atom
    /// \param at  - the number of atom in the geometry
    /// \param I   - the cell or symmetry indicies
    inline vector3<REAL> r (int at, const index & I) const {
      
      vector3<REAL> r1 = _crd[at];
      if (frac)
	r1 = cell->frac2cart(r1);
      return cell->transform(r1, I);
      
    }
    
    /// \brief r_frac
    /// \param at
    /// \param I
    /// \return
    inline vector3<REAL> r_frac (int at, const index & I) const {

      vector3<REAL> r1 = _crd[at];
      r1 = cell.cart2frac(r1);
      return cell.transform(r1, I);
      
    }

      /// \brief r_frac
      /// \param at
      /// \param I
      /// \return
      inline vector3<REAL> r_frac (int at) const {

        vector3<REAL> r1 = _crd[at];
        r1 = cell.cart2frac(r1);
        return cell.transform(r1, index::D(DIM()).all(0));

      }


      /// \brief real-space position of an atom
      /// \param ai  - complex index; ai[0] is the number of atom in the geometry,
      ///   ai[1:DIM] are the cell or symmetry indicies
      inline vector3<REAL> r (const index & ai) const {
        return cell->transform(_crd[ai(0)], ai.tail(1) );
      }

      //! \brief The synonym for r(at)
      inline vector3<REAL> pos (int at) const
      { return r(at); }

      //! \brief The synonym for r(at,I)
      inline vector3<REAL> pos (int at, const index & I ) const
      { return r(at,I); }

      //! \brief The synonym for r(ai)
      inline vector3<REAL> pos (const index & ai) const {
        return cell->transform(_crd[ai(0)], ai.tail(1) );
      }

      inline bool shadow (int at) const {
        return _shadow[at];
      }

      // ------------------- Typetable for atoms ------------------------------
    std::shared_ptr<geometry_type_table<REAL>> _typetable;

    std::shared_ptr<geometry_type_table<REAL>> typetable(){
      if (_typetable)
	return _typetable;
      else
	ValueError("Initialize typetable first!");
      return nullptr;
    }

    void build_typetable()
    {
      _typetable  = std::make_shared<geometry_type_table<REAL>>(geometry_type_table<REAL>(*this));
      _typetable->build();
      add_observer(_typetable);
    }

    void set_typetable(std::shared_ptr<geometry_type_table<REAL>> tt){
      _typetable = tt;
      add_observer(_typetable);
    }



      // --------------- Constructors & destructors --------------------

    private:

      void init_default () {

        //auto_update_types = true;
        frac = false;
        //auto_symmetrize = false;
        has_observers = false;
        //default_symmetrize_radius = 0e0;
	//if (cell == nullptr)
	tol_geom = cell -> tol_transl;
        _atm.reserve(GEOM_DEFAULT_RESERVE_AMOUNT);
        _crd.reserve(GEOM_DEFAULT_RESERVE_AMOUNT);
        _shadow.reserve(GEOM_DEFAULT_RESERVE_AMOUNT);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
        py_atoms.bind(this);
        py_shadow.bind(this);
        py_x.bind(this);
        py_y.bind(this);
        py_z.bind(this);
        py_coords.bind(this);
        py_symmrad.bind(this);
	py_observer.bind(this);
#endif

      }

    public:

    geometry (CELL & __cell, const STRING_EX & __name = ""){
      cell = std::shared_ptr<CELL>(&__cell, [](CELL*p){});
      
      init_default();
      //DIM = cell.DIM;
      name = __name;
      
    }


    geometry (int dim = 0, const STRING_EX & __name = "") : cell(std::make_shared<periodic_cell<REAL>>(periodic_cell<REAL>(dim))) {

        init_default();
        //DIM = dim;
        name = __name;

      }

    geometry (const geometry<REAL> & g) :
      cell(g.cell), _atm(g._atm), _crd(g._crd),
      _shadow(g._shadow), _typetable(g._typetable),
      observers(g.observers) {
      
      init_default();
      
      //DIM = g.DIM;
      
      // options
      //auto_update_types = g.auto_update_types;
      frac = g.frac;
      if (_typetable !=nullptr)
	{
	  _typetable->auto_symmetrize = g._typetable->auto_symmetrize;	 
	  _typetable->default_symmetrize_radius = g._typetable->default_symmetrize_radius;
	}
      tol_geom = cell->tol_transl;
      
      // observers
      has_observers = g.has_observers;
      
      
    }
    
    ~geometry () {
      //fixme
    }
    
      /*
    void copy(const geometry<DIM, REAL> &G)
    {
      // fixme - what is necessary to copy?
      clear();
      trnf = G.trnf;
      atm = G.atm;
      _crd = G._crd;
      _shadow = G._shadow;
      geomtol = G.geomtol;
      frac = G.frac;
    }
    */
      // ----------------------- Managing observers -----------------------

    void add_observer (const std::shared_ptr<DEP>  &d) {

        auto it = std::find(observers.begin(), observers.end(), d);

        if (it != observers.end()) {
            has_observers = true;
            return;
          }

        observers.push_back(d);
        has_observers = true;

      }

    void remove_observer (const std::shared_ptr<DEP> &d) {

        auto i = observers.begin();

        while(i != observers.end()) {
            if (*i == d) {
                observers.erase(i);
                break;
              }
            i++;
          }

        if (observers.size()==0)
          has_observers = false;

      }

    std::shared_ptr<DEP> get_observer(int i){
      return observers[i];
    }

    int n_observers()const{return observers.size();}
    
    void  set_observer(int i,const  std::shared_ptr<DEP>  &d){
      observers[i] = d;
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)    
    py_indexed_property<SELF,std::shared_ptr<DEP>,int, &SELF::get_observer, &SELF::set_observer> py_observer;

#endif
      // ----------------------- Manipulations with atoms -----------------------

    protected:

      inline void _add (const STRING_EX & a, const vector3<REAL> & r1) {

        //std::cerr << "geometry::add entry\n";

        vector3<REAL> r2 = r1;
        if (_typetable!=nullptr && _typetable->auto_symmetrize){
            REAL rad = _typetable->symmetrize_radius_by_label(a);
            r2 = cell->symmetrize(r1,rad);
          }

        if (has_observers)
          for (int i=0; i<observers.size(); i++)
            observers[i]->added(before, a, r2);

        _atm.push_back(a);
        _crd.push_back(r2);
        _shadow.push_back(false);

        //if (auto_update_types)
        //  _type_table.push_back(define_type(a));

        if (has_observers)
          for (int i=0; i<observers.size(); i++)
            observers[i]->added(after, a, r2);

      }

      inline void _erase (int at) {

        if (has_observers)
          for (int j=0; j<observers.size(); j++)
            observers[j]->erased(at,before);

        _atm.erase(_atm.begin()+at);
        _crd.erase(_crd.begin()+at);
        _shadow.erase(_shadow.begin()+at);
        //if (auto_update_types)
        //  _type_table.erase(_type_table.begin()+at);

        if (has_observers)
          for (int j=0; j<observers.size(); j++)
            observers[j]->erased(at,after);

      }

      inline void _insert (int at, const STRING_EX & a, const vector3<REAL> & r1) {

        vector3<REAL> r2 = r1;
        if (_typetable!=nullptr && _typetable->auto_symmetrize)
          r2 = cell->symmetrize(r1,_typetable->symmetrize_radius_by_label(a));

        if (has_observers)
          for (int j=0; j<observers.size(); j++)
            observers[j]->inserted(at,before, a, r2);

        _atm.insert(_atm.begin()+at,a);
        _crd.insert(_crd.begin()+at,r2);
        _shadow.insert(_shadow.begin()+at,false);

	// if (auto_update_types)
        //  _type_table.insert(_type_table.begin()+at,define_type(a));

        if (has_observers)
          for (int j=0; j<observers.size(); j++)
            observers[j]->inserted(at,after, a, r2);

      }

    inline void _change (int at, const STRING_EX & a1, const vector3<REAL> & r1) {
      //std::cout << "_change\n";
      if (has_observers)
	for (int i = 0; i < observers.size(); i++)
	  observers[i]->changed(at, before, a1, r1);
      //std::cout << "observers before done\n";
      _crd[at] = r1;
      _atm[at] = a1;
      
      //if (auto_update_types && at < _type_table.size())
      //  _type_table[at] = define_type(a1);
      //std::cout << "atom changed\n";

      if (has_observers)
	for (int i = 0; i < observers.size(); i++)
	  observers[i]->changed(at, after, a1, r1);
      //std::cout << "observers after done\n";
    }

  public:
    
    virtual void add (const STRING_EX & a, const vector3<REAL> & r1)
    { _add(a,r1); }
    
    virtual void add (const STRING_EX &a, const REAL _x, const REAL _y, const REAL _z)
    {  add(a,{_x,_y,_z}); }
    
    virtual void insert (int at, const STRING_EX & a, const vector3<REAL> & r1)
    { _insert(at,a,r1); }
    
    void insert (int at, const STRING_EX & a, const REAL _x, const REAL _y, const REAL _z)
    { insert(at,a,{_x,_y,_z}); }
    
    virtual void erase (int at)
    {  _erase(at); }
    
    void erase(std::vector<int>  ats){
      std::sort(ats.begin(), ats.end(), std::greater<int>());
      for (int at:ats)
	erase(at);
    }
    
    inline void shadow (int at, bool sh) {
      
      for (int i=0; i<observers.size(); i++)
	observers[i]->shaded(at, before, sh);
      _shadow[at] = sh;
      for (int i=0; i<observers.size(); i++)
	observers[i]->shaded(at, after, sh);
      
    }

    virtual void change (int at, const STRING_EX & a1, const vector3<REAL> & r1)
    { _change(at,a1,r1); }
    
    virtual void change_pos (int at, const vector3<REAL> & r1)
    { _change(at, _atm[at], r1); }
    
    void replace(const STRING_EX& at1, const STRING_EX& at2,
		 const std::function<bool(const geometry<REAL> &, int)> & where){
    for (int i=0; i< nat(); i++)
      if (atom(i)==at1 && where(*this,i))
	change(i,at2,coord(i));
    }
    
      virtual void reorder (const std::vector<int> & ord) {

        for (int i=0; i<observers.size(); i++)
          observers[i]->reordered(ord, before);
        // fixme - might be inefficient for large molecules

        std::vector<STRING_EX> __atm(_atm);
        std::vector<vector3<REAL> > __crd(_crd);
        std::vector<Bool> __shadow(_shadow);

        //bool reorder_types = (_type_table.size() == size());

        for (int i=0; i<size(); i++) {

            _atm[i] = __atm[ord[i]];
            _crd[i] = __crd[ord[i]];
            _shadow[i] = __shadow[ord[i]];

          }
	if (_typetable)
	  _typetable->reorder_types(ord);

        for (int i=0; i<observers.size(); i++)
          observers[i]->reordered(ord, after);

      }

      void sort (const std::function<REAL(const geometry<REAL> &, int)> & key) {

        std::vector<int> ord(size());
        for (int i = 0; i < size(); i++)
          ord[i] = i;
        std::sort(ord.begin(), ord.end(),
                  [&](const int& a, const int& b){
            return (key(*this,a) < key(*this,b));
          });
        reorder(ord);

      }

      virtual void clear () {

        _crd.clear();
        _atm.clear();
        _shadow.clear();
	if (_typetable)
	  _typetable->clear();

      }

      virtual void get_fields (int j, std::vector<fieldtypes> & v) const {

        if (j<0) j+=nat();
        if (j<0 || j>= nat())
          IndexError("xgeometry::get_fields: index out of range");

        v.clear();
        v.push_back(atom(j));
        vector3<REAL> rl = coord(j);
        v.push_back(rl(0));
        v.push_back(rl(1));
        v.push_back(rl(2));

      }

      virtual void set_fields (int j, const std::vector<fieldtypes> & v) {

        if (j<0) j+=nat();

        if (j<0 || j>= nat())
          IndexError("xgeometry::set_fields: index out of range");

        if (v.size()!=4)
          IndexError("geometry::set_fields: wrong number of fields, must be 4");

        change(j, std::get<STRING_EX>(v[0]), qpp::vector3<REAL>(std::get<REAL>(v[1]),
								std::get<REAL>(v[2]), std::get<REAL>(v[3])));

      }


      std::vector<fieldtypes> operator[](int j) {
        std::vector<fieldtypes> v;
        get_fields(j,v);
        return v;

      }

      virtual void write (std::basic_ostream<CHAR_EX,TRAITS> &os, int offset=0) const {

        for (int k=0; k<offset; k++) os << " ";
        os << "geometry";
        if (name != "")
          os << " " << name;
        os << "(" << DIM() << "d,atom,x,y,z){\n";

        /*
      os << "(";
      for (int i=0; i<n_param(); i++)
        {
          param(i)->write(os);
          if (i<n_param()-1)
            os << ",";
        }
      os << "){\n";
      */

        cell->write(os,offset+4);

        for (int i=0; i<size(); i++){
            for (int k=0; k<offset+2; k++) os << " ";
            //TODO Migrate to fmt
            // os << _atm[i] << boost::format(" %11.6f %11.6f %11.6f\n") %
            //_crd[i].x() % _crd[i].y() % _crd[i].z();

          }

        for (int k=0; k<offset; k++) os << " ";
        os << "}\n";

      }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      /*
    inline CELL & py_getcell()
    { return *cell; }

    inline void py_setcell(CELL & cl)
    {
      cell = & cl;
      DIM = cell -> DIM;
    }
    */
      // --------------------------------------------------------

      inline vector3<REAL> py_pos1 (int at) const
      { return r(at); }

      inline vector3<REAL> py_pos2 (int at, const index & I ) const
      { return r(at,I); }

      inline vector3<REAL> py_pos3 (const index & AI ) const
      { return r(AI); }


      // --------------------------------------------------------------------
      // Implementing some sexual perversions to make nice&simple properties
      // (atom, coord, x, y, z, shadow)

      STRING_EX py_getatom (int i) {
        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::atom::Index out of range");
        return atom(i);
      }

      void py_setatom (int i, const STRING_EX & a) {
        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::atom::Index out of range");
        change(i,a,coord(i));
      }

      py_indexed_property< SELF, STRING_EX, int,
      &SELF::py_getatom, &SELF::py_setatom> py_atoms;

      bool py_getshadow (int i) {
        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::shadow::Index out of range");
        return shadow(i);
      }

      void py_setshadow (int i, const bool & b) {
        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::shadow::Index out of range");
        shadow(i,b);
      }

      py_indexed_property< SELF, bool, int,
      &SELF::py_getshadow, &SELF::py_setshadow> py_shadow;

      //    py_array_vector<REAL> py_coords;

      REAL py_getcoord (int i, int d) {
        if (i<0) i+=nat();
        if (i<0 || i>=nat() || d<0 || d>2 ) IndexError("geometry::coord::Index out of range");
        return coord(i)(d);
      }

      void py_setcoord (int i, int d, const REAL & c) {

        if (i<0) i+=nat();
        if (i<0 || i>=nat() || d<0 || d>2 ) IndexError("geometry::coord::Index out of range");
        vector3<REAL> r1 = coord(i);
        r1(d) = c;
        change(i,atom(i),r1);

      }

      vector3<REAL> py_getvec (int i) {

        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::coord::Index out of range");
        return coord(i);

      }

      void py_setvec (int i, const vector3<REAL> & v) {

        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::coord::Index out of range");
        change(i,atom(i),v);

      }

      py_2indexed_property<SELF,vector3<REAL>,REAL, int, &SELF::py_getvec, &SELF::py_setvec,
      &SELF::py_getcoord, &SELF::py_setcoord > py_coords;

      template<int d>
      REAL py_getxyz (int i) {

        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::x,y,z:: Index out of range");
        return py_getcoord(i,d);

      }

      template<int d>
      void py_setxyz (int i, const REAL & c) {

        if (i<0) i+=nat();
        if (i<0 || i>=nat()) IndexError("geometry::x,y,z:: Index out of range");
        py_setcoord(i,d,c);

      }

      py_indexed_property<SELF, REAL, int,
      &SELF::py_getxyz<0>, &SELF::py_setxyz<0> > py_x;
      py_indexed_property<SELF, REAL, int,
      &SELF::py_getxyz<1>, &SELF::py_setxyz<1> > py_y;
      py_indexed_property<SELF, REAL, int,
      &SELF::py_getxyz<2>, &SELF::py_setxyz<2> > py_z;

      REAL py_getsymmrad (const char * at) {

        int t = _typetable->type_of_atom(at);

        if (t==-1)
          KeyError("geometry::symmrad::get_symmrad: Atom not specified");
        else
          return _typetable->symmetrize_radius(t);
	return 0e0;
      }

      void py_setsymmrad (const char *at, const REAL &r)
      { _typetable->set_symmetrize_radius_by_label(at,r);  }

      py_indexed_property<SELF,REAL,const char *,
      &SELF::py_getsymmrad, &SELF::py_setsymmrad> py_symmrad;


      // --------------------------------------------------------------------

      inline int py_typeofatom (const STRING_EX & at)
      {  return _typetable->type_of_atom(at); }

      // --------------------------------------------------------------------

      static void py_props (py::module m, const char * pyname) {
        std::string sPropNameAtom =
            fmt::format("{0}_{1}",pyname,"idx_prop_atom");
        py_indexed_property<
            SELF,
            STRING_EX,
            int,
            &SELF::py_getatom,
            &SELF::py_setatom>::py_export(m, sPropNameAtom.c_str());

        std::string sPropNameShadow =
            fmt::format("{0}_{1}",pyname,"idx_prop_shadow");
        py_indexed_property<
            SELF,
            bool,
            int,
            &SELF::py_getshadow,
            &SELF::py_setshadow>::py_export(m, sPropNameShadow.c_str());


        std::string sPropNameCoord =
            fmt::format("{0}_{1}",pyname,"2idx_prop_coord");
        py_2indexed_property<
            SELF,
            vector3<REAL>,
            REAL,
            int,
            &SELF::py_getvec,
            &SELF::py_setvec,
            &SELF::py_getcoord,
            &SELF::py_setcoord >::py_2export(m, sPropNameCoord.c_str());

        std::string sPropNameXYZ0 =
            fmt::format("{0}_{1}",pyname,"idx_prop_xyz0");
        py_indexed_property<
            SELF,
            REAL,
            int,
            &SELF::py_getxyz<0>,
            &SELF::py_setxyz<0> >::py_export(m, sPropNameXYZ0.c_str());

        std::string sPropNameXYZ1 =
            fmt::format("{0}_{1}",pyname,"idx_prop_xyz1");
        py_indexed_property<
            SELF,
            REAL,
            int,
            &SELF::py_getxyz<1>,
            &SELF::py_setxyz<1> >::py_export(m, sPropNameXYZ1.c_str());

        std::string sPropNameXYZ2 =
            fmt::format("{0}_{1}",pyname,"idx_prop_xyz2");
        py_indexed_property<
            SELF,
            REAL,
            int,
            &SELF::py_getxyz<2>,
            &SELF::py_setxyz<2> >::py_export(m, sPropNameXYZ2.c_str());

        std::string sPropNameSymRad =
            fmt::format("{0}_{1}",pyname,"idx_prop_symrad");
        py_indexed_property<
            SELF,
            REAL,
            const char *,
            &SELF::py_getsymmrad,
            &SELF::py_setsymmrad>::py_export(m, sPropNameSymRad.c_str());

	std::string sPropObsrv=fmt::format("{0}_{1}",pyname,"idx_prop_obsrv");
	py_indexed_property<SELF,std::shared_ptr<DEP>,int, &SELF::get_observer, &SELF::set_observer>::
	  py_export(m,sPropObsrv.c_str());
      }

      // --------------------------------------------------------------------

      void py_add1 (const STRING_EX & a, const vector3<REAL> & r1)
      { add(a,r1); }

      void py_add2 (STRING_EX a, const REAL _x, const REAL _y, const REAL _z)
      { add(a,_x,_y,_z); }

      bool py_check_axyz (const py::list & l) {

        return
            py::len(l) == 4 &&
            py::isinstance<py::str>(l[0]) &&
            py::isinstance<py::float_>(l[1]) &&
            py::isinstance<py::float_>(l[2]) &&
            py::isinstance<py::float_>(l[3]);
      }

      virtual void py_add_list (const py::list &l) {

        if (!py_check_axyz(l))
          TypeError("geometry::Invalid list. List must be [atom,x,y,z]");

        add(py::cast<py::str>(l[0]),
            py::cast<py::float_>(l[1]),
            py::cast<py::float_>(l[2]),
            py::cast<py::float_>(l[3]));

      }

      void py_add3 (const py::list & l) {

        if (py::len(l)==0) return;
        if (py::isinstance<py::list>(l[0]))
          // list of lists
          for (int i=0; i<py::len(l); i++)
            py_add_list(py::cast<py::list>(l[i]) );
        else
          py_add_list(l);

      }

      void py_insert1 (int i, const STRING_EX & a, const vector3<REAL> & r1) {

        if (i<0) i+=nat();
        if (i<0 || i>nat()) IndexError("geometry::Index out of range");

        insert(i,a,r1);

      }

      void py_insert2 (int i, const STRING_EX & a,
                       const REAL _x,
                       const REAL _y,
                       const REAL _z) {

        if (i<0) i+=nat();
        if (i<0 || i>nat()) IndexError("geometry::Index out of range");
        insert(i,a,_x,_y,_z);

      }

      virtual void py_insert_list (int i, const py::list & l) {

        if (!py_check_axyz(l))
          TypeError("geometry::Invalid list. List must be [atom,x,y,z]");

        insert(i,
               py::cast<STRING_EX>(l[0]),
            py::cast<REAL>(l[1]),
            py::cast<REAL>(l[2]),
            py::cast<REAL>(l[3]));

      }

      void py_insert3 (int at, const py::list & l) {
        if (at<0) at+=nat();
        if (at<0 || at>nat()) IndexError("geometry::Index out of range");

        if (py::len(l)==0) return;
        if (py::isinstance<py::list>(l[0]))
          // list of lists
          for (int i=0; i<py::len(l); i++)
            py_insert_list(at+i, py::cast<py::list>(l[i]) );
        else
          py_insert_list(at, l);

      }

      void py_erase (int at) {

        if (at<0) at+=nat();
        if (at<0 || at>=nat()) IndexError("geometry::Index out of range");

        erase(at);
      }

    void py_erase2(const std::vector<int> ats){
      erase(ats);
    }

      virtual py::list py_getitem (int i) const {

        if (i<0)
          i+=nat();

        if (i<0 || i>=nat())
          IndexError("geometry::Index out of range");

        py::list l;
        l.append(atom(i));
        l.append(coord(i).x());
        l.append(coord(i).y());
        l.append(coord(i).z());
        return l;

      }

      virtual void py_setitem (int i, const py::list & l) {

        if (i<0)
          i+=nat();

        if (i<0 || i>=nat())
          IndexError("geometry::Index out of range");

        if (!py_check_axyz(l))
          TypeError("geometry::Invalid list. List must be [atom,x,y,z]");

        STRING_EX a1 = py::cast<STRING_EX>(l[0]);
        REAL x1 = py::cast<REAL>(l[1]),
            y1 = py::cast<REAL>(l[2]),
            z1 = py::cast<REAL>(l[3]);

        change(i,a1,{x1,y1,z1});

      }

#endif

  };

  template< class REAL>
  REAL geometry<REAL>::tol_geom_default = 1e-6;

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

  template<class REAL>
  struct py_geometry_observer : geometry_observer<REAL> {

      using geometry_observer<REAL>::geometry_observer;
      void added (before_after s, const STRING_EX & a, const vector3<REAL> & v)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              added,
              s,a,v);
      }

      void inserted (int at, before_after s, const STRING_EX & a,
                     const vector3<REAL> & v)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              inserted,
              at,s,a,v);
      }

      void changed (int at, before_after s, const STRING_EX & a,
                    const vector3<REAL> & v)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              changed,
              at,s,a,v);
      }

      void erased (int at, before_after s)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              erased,
              at,s);
      }

      void shaded (int at, before_after s, bool h)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              shaded,
              at,s,h);
      }

      void reordered (const std::vector<int> & ord, before_after s)
      override {
        PYBIND11_OVERLOAD_PURE(
              void,
              geometry_observer<REAL>,
              reordered,
              ord, s);
      }

  };
#endif

} // namespace qpp

#endif

