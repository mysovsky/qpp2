#ifndef QPP_NGBR_H
#define QPP_NGBR_H

#include <typeinfo>
#include <utility>
#include <set>
#include <algorithm>
#include <cmath>
#include <geom/geom.hpp>
#include <string>
#include <chrono>
#include <thread>

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
//#include <pybind11/stl.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif

namespace qpp {

  // ------------------------------------------------------------


  template <class REAL>
  class bonding_table {

      struct _covrad_record {
          STRING_EX at;
          REAL d;

          _covrad_record(const STRING_EX & _at, REAL _d){
            at = _at;
            d = _d;
          }

      };

      std::vector<_covrad_record> _covrads;

      inline bool _covrad_match(const STRING_EX & at, int i) const {
        return at == _covrads[i].at;
      }

      inline int _find_covrad(const STRING_EX & at) const {
        for (int i=0; i<_covrads.size(); i++)
          if (_covrad_match(at,i))
            return i;
        return -1;
      }

      // -------------------------------------------------

      struct _ngbr_record {
          STRING_EX at1, at2;
          REAL d;

          _ngbr_record(const STRING_EX & _at1, const STRING_EX & _at2, REAL _d){
            at1 = _at1;
            at2 = _at2;
            d = _d;
          }
      };

      std::vector<_ngbr_record> _records;

      inline bool _record_match(const STRING_EX & at1,
                                const STRING_EX & at2, int i) const{
        return ( at1 == _records[i].at1 && at2 == _records[i].at2 ) ||
            ( at2 == _records[i].at1 && at1 == _records[i].at2 );
      }

      inline int _find_record(const STRING_EX & at1, const STRING_EX & at2) const{
        for (int i=0; i<_records.size(); i++)
          if (_record_match(at1,at2,i))
            return i;
        return -1;
      }

    public:

      typedef bonding_table<REAL> SELF;

      STRING_EX name;

      REAL default_distance;

      REAL covrad(const STRING_EX & at){
        int i = _find_covrad(at);
        return i>-1 ? _covrads[i].d : 0e0;
      }

      void set_covrad(const STRING_EX  & at, const REAL & d){
        int i = _find_covrad(at);
        if (i>-1)
          _covrads[i].d = d;
        else
          _covrads.push_back(_covrad_record(at,d));
      }

      REAL pair(const STRING_EX   at1, const STRING_EX   at2){
        int i = _find_record(at1,at2);
        return i>-1 ? _records[i].d : 0e0;
      }

      void set_pair(const STRING_EX   at1, const STRING_EX   at2, const REAL &d){
        int i = _find_record(at1,at2);
        if (i>-1)
          _records[i].d = d;
        else
          _records.push_back(_ngbr_record(at1,at2,d));
      }

      REAL distance(const STRING_EX & at1, const STRING_EX   at2){
        int i = _find_record(at1,at2);
        if (i>-1)
          return _records[i].d;

        i = _find_covrad(at1);
        int j = _find_covrad(at2);

        if ( i>-1 && j>-1 )
          return _covrads[i].d + _covrads[j].d;

        return default_distance;
      }


      void clear(){
        _records.clear();
        _covrads.clear();
        default_distance = 0e0;
      }

      virtual void write(std::basic_ostream<CHAR_EX,TRAITS> &os, int offset=0) const{
        for (int k=0; k<offset; k++) os << " ";
        os << "bonding_table";
        if (name != "")
          os << " " << name;
        os << "{\n";

        for (int k=0; k<offset+2; k++) os << " ";
        os << "default = " << default_distance << ";\n";

        for (const auto & c : _covrads){
            for (int k=0; k<offset+2; k++) os << " ";
            os << "covrad(" << c.at << ") = " << c.d << ";\n";
          }

        for (const auto & p : _records){
            for (int k=0; k<offset+2; k++) os << " ";
            os << "pair(" << p.at1 << "," << p.at2 << ") = " << p.d << ";\n";
          }

        for (int k=0; k<offset; k++) os << " ";
        os << "}\n";
      }

      void merge(const bonding_table<REAL> & bt){
        default_distance = std::max(default_distance, bt.default_distance);
        for (int i=0; i<bt._records.size(); i++){
            int j = _find_record( bt._records[i].at1, bt._records[i].at2);
            if (j==-1)
              _records.push_back(bt._records[i]);
            else
              _records[j].d = std::max( _records[j].d, bt._records[i].d);
          }

        for (int i=0; i<bt._covrads.size(); i++){
            int j = _find_covrad( bt._covrads[i].at );
            if (j==-1)
              _covrads.push_back(bt._covrads[i]);
            else
              _covrads[j].d = std::max( _covrads[j].d, bt._covrads[i].d);
          }
      }


#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      // --------------- PYTHON -------------------------------

      py_indexed_property<SELF, REAL, const STRING_EX &,
      &SELF::covrad, &SELF::set_covrad> py_covrad;

      REAL py_getpair(const STRING_EX){ return -1;}
      void py_setpair(const STRING_EX, const REAL &){}

      REAL py_getpair2(const STRING_EX at1, const STRING_EX  at2){
        return pair(at1,at2);
      }

      void py_setpair2(const STRING_EX at1, const STRING_EX  at2, const REAL & d){
        set_pair(at1,at2,d);
      }

      py_2indexed_property<SELF, REAL, REAL, const STRING_EX,
      &SELF::py_getpair, &SELF::py_setpair,
      &SELF::py_getpair2, &SELF::py_setpair2 > py_pair;

      py::dict to_dict(){
        py::dict d;
        d["default"] = default_distance;

        for (auto i = _covrads.begin(); i!=_covrads.end(); i++)
          d[i->at.c_str()] = i->d;

        for (auto i = _records.begin(); i!=_records.end(); i++)
          d[py::make_tuple(i->at1, i->at2)] = i->d;

        return d;
      }

      void from_dict( py::dict & d){

        for(auto p : d){
            //std::vector p = py::cast<std::vector>(pnc);
            REAL val = this->default_distance;
            if (!py::isinstance<py::float_>(p.second))
              TypeError("bonding_table::from_dict: "
                        "Invalid dictionary, values must be of real type");
            else val = py::cast<REAL>(p.second);

            //bp::extract<STRING> ks(p[0]);
            if (py::isinstance<py::str>(p.first)){
                STRING_EX s = py::cast<STRING_EX>(p.first);
                if (s=="default")
                  default_distance = val;
                else
                  set_covrad(s,val);
              }

            else if (py::cast<py::tuple>(p.first)){
                py::tuple t = py::cast<py::tuple>(p.first);
                if (py::len(t)!=2)
                  TypeError("bonding_table::from_dict: Invalid dictionary, "
                            "tuple key can contain a pair of atoms only");

                if (py::isinstance<py::str>(t[0]) &&
                    py::isinstance<py::str>(t[1]))
                  set_pair(py::cast<STRING_EX>(t[0]),
                      py::cast<STRING_EX>(t[1]), val);
                else
                  TypeError("bonding_table::from_dict: Invalid dictionary,"
                            " tuple key can contain a pair of atoms only");
              }
            else
              TypeError("bonding_table::from_dict: Invalid dictionary, "
                        "key can be \'default\', an atom, or a pair of atoms");
          }
      }

      static void py_export(py::module m, const char * pyname){
        std::string sPropNameCovRad =
            fmt::format("{0}_{1}",pyname,"idx_prop_covrad");
        py_indexed_property<SELF, REAL, const STRING_EX & , &SELF::covrad,
            &SELF::set_covrad>::py_export(m, sPropNameCovRad.c_str());

        std::string sPropNamePair =
            fmt::format("{0}_{1}",pyname,"idx_prop_pair");
        py_2indexed_property<SELF,REAL,REAL, const STRING_EX , &SELF::py_getpair,
            &SELF::py_setpair,&SELF::py_getpair2, &SELF::py_setpair2 >
            ::py_2export(m, sPropNamePair.c_str(),true);

        py::class_<bonding_table<REAL> >(m, pyname)
            .def( py::init<>())
            .def_readwrite("default", &SELF::default_distance,
                           "The default bonding distance between any atoms (r/w)" )
            .def_readwrite("covrad",  &SELF::py_covrad,
                           "Covalent radius for atom. "
                           "Usage: covrad[at] with string at(r/w)")
            .def_readwrite("pair",    &SELF::py_pair,
                           "Pair specific bonding distance. "
                           "Usage: pair[at1,at2] "
                           "with string at1,at2 (r/w)")
            .def("distance", &SELF::distance,
                 "Bonding distance between two atoms. "
                 "Usage: distance(at1,at2)"
                 " with string at1,at2")
            .def("to_dict", & SELF::to_dict,
                 "Output the content of bonding"
                 " table as dictionary. Usage: to_dict()")
            .def("from_dict", & SELF::from_dict,
                 "Fill the bonding table from dictionary. "
                 "Usage: from_dict(dict)")
            .def("clear", & SELF::clear, "Erase all data records")
            .def("merge",  & SELF::merge, "Merge another bonding table"
                                          " to this one");

      }

#endif

      bonding_table(const STRING_EX & __name = ""){
        name = __name;
        default_distance = 0e0;
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
        py_covrad.bind(this);
        py_pair.bind(this);
#endif
      }


  };

  // ------------------- Neighbours table -----------------------

  template <class REAL>
  class neighbours_table : public geometry_observer<REAL> {

      typedef neighbours_table<REAL> SELF;

      // The periodicity the geometry
      int DIM;

      // The geometry for which neighbours table is
      geometry<REAL> * geom;

      // The bonding distances data
      bonding_table<REAL> *btbl;

      // Internal bonding distances table
    REAL* _disttable{nullptr};

      // Number of atomic types
      int ntp;

      // Bonding distance between i-th and j-th atomic types
      inline REAL distance(int i, int j) {return _disttable[i*ntp+j];}

      inline void resize_disttable(){
        //      if (_disttable != nullptr)
        //        delete []_disttable;
        // delete []_disttable;
	if (_disttable != nullptr)
	  delete []_disttable;
	ntp = geom->typetable()->n_types();
	if (ntp < 1)
	  geom-> build_typetable();
	//IndexError("Number of atomic types is zero. Maybe, you should initialize typetable?");
	_disttable = new REAL[ntp*ntp];
      }

      void build_disttable(){
        resize_disttable();
        for (int i=0; i<ntp; i++)
          for (int j=0; j<=i; j++) {
              _disttable[ntp*i+j] = btbl->distance(geom->typetable()->atomic_type(i),geom->typetable()->atomic_type(j));
              _disttable[ntp*j+i] = btbl->distance(geom->typetable()->atomic_type(i),geom->typetable()->atomic_type(j));

            }
      }
      // ------------------------------------------

      std::vector<std::vector<index> > _table;

    public:
      // Number of neighbours of i-th atom
      inline int n(int i) const {return _table[i].size();}

      // j-th neighbour of i-th atom
      inline index table(int i, int j) const {return _table[i][j];}

      // Synonym
      inline index operator()(int i, int j) const{
        return table(i,j);
      }

      // The same but more general - i-th atom specified with full index
      inline index operator()(const index & i, int j) const{
        index res = table(i,j);
        for (int d=0; d<DIM; d++)
          res(d+1) += i(d+1);
        return res;
      }

    private:

      inline void _add_ngbr(int i, const index & j){
        bool found = false;
        for (int k = 0; k < _table[i].size(); k++ )
          if (_table[i][k]==j){
              found = true;
              break;
            }
        if (!found)
          _table[i].push_back(j);
      }

      // ------------------------------------------

      bool auto_update;

    public:
      bool auto_grainsize;

      void set_auto_update(bool au){
        if (au && !auto_update){
	  geom->add_observer(std::shared_ptr<SELF>(this));
            auto_update = true;
          }
        else if (!au && auto_update){
	  geom->remove_observer(std::shared_ptr<SELF>(this));
            auto_update = false;
          }
      }

      bool get_auto_update() const {return auto_update;}

    private:
      REAL grainsize;

      vector3<REAL> Rmin, Rmax;
      index ngrain;
      std::vector<std::vector<index> > grains;

      void translational_grain_setup(){
        // find largest neighbourung distance
        //std::cout << "trnsl grain setup\n";
	auto t1 = std::chrono::high_resolution_clock::now();
        bool ndef = true;

        for (int i=0; i<geom->size(); i++){
            vector3<REAL> r = geom->pos(i);
            if (ndef){
                Rmin = Rmax = r;
                ndef = false;
              }
            for (int j=0; j<3; j++){
                if ( Rmin(j) > r(j) )
                  Rmin(j) = r(j);
                if ( Rmax(j) < r(j) )
                  Rmax(j) = r(j);
              }
          }
        //std::cout <<  "Rmin = " << Rmin << " Rmax= " << Rmax << "\n";

        if (DIM!=0){
            //std::cout << "Extension of Rmin, Rmax\n";
            Rmin -= vector3<REAL>(1,1,1)*(grainsize+1e-5);
            Rmax += vector3<REAL>(1,1,1)*(grainsize+1e-5);
          }

        ngrain = index::D(3);
        for (int i=0; i<3; i++)
          ngrain(i) = int( (Rmax(i)-Rmin(i))/grainsize ) + 1;

        //std::cout << "grain setup finished\n" << "Rmin = " << Rmin
        //<< " Rmax= " << Rmax << " gs= " << grainsize <<
        //" ngrain = " << ngrain << "\n";
	auto t2 = std::chrono::high_resolution_clock::now();
	int t=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
	std::cout << "Grain setup finished in "<< t << " msec\n";

      }

      REAL optimal_grainsize(){
        const REAL alpha = 1, beta = 0;
        REAL Vol = (Rmax(0)-Rmin(0))*(Rmax(1)-Rmin(1))*(Rmax(2)-Rmin(2));
        return std::pow(Vol,1./3)*std::pow(alpha/geom->nat(), 1./(3*(beta+1)));

      }

      void grain_setup(){

        if (auto_grainsize){
            grainsize = REAL(0);
            int ntp = geom->typetable()->n_atom_types();
            for (int i=0; i<ntp*ntp; i++)
              if (_disttable[i] > grainsize )
                grainsize = _disttable[i];
          }

        translational_grain_setup();

        REAL opt_gs = optimal_grainsize();

        //std::cout << "gs= " << grainsize << " opt_gs= " << opt_gs << "\n";

        if (auto_grainsize && grainsize < opt_gs ){
            grainsize = opt_gs;
            for (int i=0; i<3; i++)
              ngrain(i) = int( (Rmax(i)-Rmin(i))/grainsize ) + 1;
          }
	//std::cout << "leaving grain setup\n";

      }

      inline int gidx(int i, int j, int k){
        return i + ngrain(0)*j + ngrain(0)*ngrain(1)*k;
      }

      inline int gidx(const index & I){
        return I(0) + ngrain(0)*I(1) + ngrain(0)*ngrain(1)*I(2);
      }

      inline std::vector<index> & grain(int i, int j, int k){
        return grains[gidx(i,j,k)];
      }

      inline std::vector<index> & grain(const index & I){return grains[gidx(I)];}

      inline index igrain(const index & at){
        vector3<REAL> r = geom->r(at);
        index I = index::D(3);
	//std::cout << "igrain "<< I<< "\n";
        for (int j=0; j<3; j++)
          I(j) = std::floor( (r(j)-Rmin(j))/grainsize );
	//std::cout << "igrain "<< I<< "\n";
        return I;
      }

      inline void to_grain(const index & at){
        index I=igrain(at);
	//std::cout << "to_grain"<< I<< "\n";
        if ( I(0)>=0 && I(0)<ngrain(0) && I(1)>=0 &&
             I(1)<ngrain(1) && I(2)>=0 && I(2)<ngrain(2) )
          grain(I).push_back(at);
      }

      void graining(){
        for (int i=0; i<grains.size(); i++)
          grains[i].clear();
        grains.clear();
        grains.resize(ngrain(0)*ngrain(1)*ngrain(2));
	//std::cout << "graining middle\n";

        for (int at=0; at<geom->nat(); at++)
          if (! geom->shadow(at) )
            for ( iterator I(geom->cell->begin(),geom->cell->end()); !I.end(); I++){
	      //std::cout << I.print() << " " << I.DIM <<" at="<< at<< "\n";
	      //std::cout << index({at}).cat(I)<< "\n";
              to_grain(index({at}).cat(I));
	    }
        //      debug
        
	//std::cout << "graining finished\n";
      /*
      for (iterator I({0,0,0}, ngrain-index({1,1,1})); !I.end(); I++)
        {
          std::cout << "grain" <<  I << " " << gidx(I) << "\n";

          for (int l=0; l<grain(I).size(); l++)
            std::cout << " " << grain(I)[l] << " ";
          std::cout << "\n";

        }
      */

      }

      // -------------------------------------------------------------------


    public:

      bool reference_mode;
      bool transl_mode;
    int nparr;

      neighbours_table( geometry<REAL> & g, bonding_table<REAL> & t) :
        ngrain(index::D(3)){
        btbl = &t;
        geom = & g;
        DIM = geom -> DIM();
        //_disttable.clear();
        build_disttable();
        reference_mode = false;
        transl_mode = true;
	nparr = globals::ncores;
        auto_grainsize = true;
        auto_update = false;
	if (geom->_typetable==nullptr)
	  geom->build_typetable();
      }

      ~neighbours_table() {
        //        if (_disttable != nullptr)
        //          delete [] _disttable;
        delete [] _disttable;
      }

      REAL get_grain_size(){return grainsize;}

      void set_grain_size(REAL gs){auto_grainsize = false;grainsize = gs;}

      // -------------------------------------------------------------------

      void reference_build(){

        //std::cout << "reference build\n";

        for (int i=0; i<geom->nat(); i++)
          if (!geom->shadow(i))
            for (int k=0; k<geom->nat(); k++)
              if (!geom->shadow(k))
                for (iterator I(geom->cell->begin(), geom->cell->end()); !I.end(); I++)
                  if ( (geom->pos(i) - geom->pos(k,I)).norm() <
                       distance(geom->typetable()->type(i), geom->typetable()->type(k))
                       && !( i==k && I==index::D(DIM).all(0)) )
                    _table[i].push_back(atom_index(k,I));
      }

      void translational_build(){
        std::cout << "translational build\n";
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t10 = std::chrono::high_resolution_clock::now();

	std::vector<int> atypes(geom-> nat());
	for (int i=0; i<geom-> nat(); i++ )
	  atypes[i] = geom->typetable()->type(i);
	
        index dirray[14];
        int n=0;
        for (int i=-1; i<=0; i++)
          for (int j=-1; j<=-i; j++)
            for (int k=-1; k<=-i || k<=-j; k++)
              dirray[n++] = {i,j,k};


        // debug
        /*
      std::cout << "after dirray alive\n";
      for (const index & DI : dirray)
        std::cout << DI;
      std::cout << "\n";
      */


        //index shift1 = {1,1,1}, shift2={2,2,2};
        //index shift1 = {1,1,1}, shift2={1,1,1};
        index shift1 = {0,0,0}, shift2={1,1,1};
        /*
      if (geom->DIM==0)
        {
          //shift1 = {0,0,0};
          //shift2 = {1,1,1};
          shift1 = {1,1,1};
          shift2 = {1,1,1};
        }
      */
        int Nprint=200;
        int nprint=0;
	int iftimens = 0, prevt =0;
	auto t0 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

        for (iterator I(shift1, ngrain-shift2); !I.end(); I++){
            int g1 = gidx(I);
            if (grains[g1].size()>0)
              for (const index & DI : dirray){
                  index J = I + DI;
                  if ( J(0)<0 || J(0)>=ngrain(0) || J(1)<0 || J(1)>=ngrain(1) ||
                       J(2)<0 || J(2)>=ngrain(2))
                    continue;
                  int g2 = gidx(J);

                  if (++nprint == Nprint){
                      nprint =0;
                      //std::cout << I << " = " << g1 << " + " <<
                      //DI << " = " << J << " = " << g2 << " grains.size= "
                      //<< grains.size() << "\n";
                    }

                  for (int c2 = 0; c2 < grains[g2].size(); c2++)
                    for (int c1 = 0; c1 < ( g1==g2? c2 : grains[g1].size()); c1++){
		      t0 = std::chrono::high_resolution_clock::now();
                        index at1 = grains[g1][c1];
                        index at2 = grains[g2][c2];
                        REAL r = (geom->pos(at1) - geom->pos(at2)).norm();

                        /*
                      if (at1==31|| at2==31)
                        {
                          std::cout << "at1 = " << at1 << " grain " << I <<
 " at2= "<< at2 << " grain " << J <<  " DI =" << DI <<"\n";
                          std::cout << "type1= " << geom->type(at1) <<
  " type2= " << geom->type(at2) << " d= "
<< distance(geom->type(at1), geom->type(at2)) << " r= " << r << "\n";
                        }
                      */
			t1 = std::chrono::high_resolution_clock::now();
                        if ( r <= distance(atypes[at1], atypes[at2])){
                            if ( at1.tail(1) == index::D(DIM).all(0) ){
                                _add_ngbr(at1,at2);
                                index at= at1;
                                at.tail(1) -= at2.tail(1);
                                /*
                              for (int dd=0; dd<DIM; dd++)
                                at.setcell(dd,-at2.cell(dd));
                              */
                                _add_ngbr(at2,at);

                                /*
                              if (at1==31 || at1==41)
                                std:: cout << "added:" << int(at1) <<
                                at2 << " added " << int(at2) << at << "\n";
                              */

                              }
                            else if ( at2.tail(1) == index::D(DIM).all(0) ){
                                _add_ngbr(at2,at1);
                                index at= at2;
                                at.tail(1) -= at1.tail(1);
                                _add_ngbr(at1,at);

			    }
			    t2 = std::chrono::high_resolution_clock::now();
			    int dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
			    int dt0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
			    iftimens +=dt;
			    prevt+=dt0;
			  

                          }
                      }
                }
          }
	t2 = std::chrono::high_resolution_clock::now();
	int t=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t10).count();
	std::cout << "Translational build finished in "<< t << " msec\n";
	std::cout << "central if takes " << iftimens << " ns , the rest "<< prevt << " ns\n";
      }

    struct ngbfinder_thread{
      int n1,n2;
      bool imfinished;
      std::vector<std::pair<int,index>> my_table;
      std::vector<std::pair<int,int> > &pairlist;
      neighbours_table<REAL> &self;
      ngbfinder_thread(int _n1, int _n2,
		       std::vector<std::pair<int,int>> & _pl,
		       neighbours_table<REAL> &_self): n1(_n1),n2(_n2),
						       imfinished(false), pairlist(_pl),self(_self)
      {}
      /*
	inline void my_add_ngbr(int i, const index & j ){	  
	bool found = false;
	for (int k = 0; k < _table[i].size(); k++ )
	if (_table[i-n1][k]==j){
	found = true;
	break;
	}
	if (!found)
	_table[i-n1].push_back(j);
	}*/
      void printme(){
	std::cout << "ngbs building thread from "<< n1 <<" to " << n2 <<"\n";
	std::cout << "pairlist " << (&pairlist) << " of len "<< pairlist.size()<< "\n";
      }
      void operator()(){
	for (int i=n1; i<=n2; i++){
	  auto  gp = pairlist[i];
	  int g1 = gp.first, g2 = gp.second;
	  for (int c1 = 0; c1 < self.grains[g1].size(); c1++)
	    for (int c2 = 0; c2 < (g1==g2 ? c1 : self.grains[g2].size()); c2++){
	      index at1 = self.grains[g1][c1];
	      index at2 = self.grains[g2][c2];

	      bool skip = self.geom->shadow(at1) || self.geom->shadow(at2);
	      if (!skip)
		skip = at1.tail(1) != index::D(self.DIM).all(0)
		  && at2.tail(1) != index::D(self.DIM).all(0);
	      if (skip)
		continue;

	      REAL r = (self.geom->pos(at1) - self.geom->pos(at2)).norm();
		
	      if ( r <= self.distance(self.geom->typetable()->type(at1),
				      self.geom->typetable()->type(at2))){
		if ( at1.tail(1) == index::D(self.DIM).all(0) )
		  my_table.push_back(std::pair<int,index>(at1, at2));
		if ( at2.tail(1) == index::D(self.DIM).all(0) )
		  my_table.push_back(std::pair<int,index>(at2, at1));
	      }
		
	    }
	}
	imfinished=true;
	printme();
	//std::cout << "I am finished "<< imfinished << " size " << my_table.size() << "\n";
      }
    };
    void build_parallel(){
      // graining is fast - no paralellization
      grain_setup();
      graining();
      std::set<std::pair<int,int> > gpairs;
      std::vector<char> marked(ngrain(0)*ngrain(1)*ngrain(2), false);

      for (int at=0; at<geom->nat(); at++)
	if (!geom->shadow(at)){
	  index I = igrain(index({at}).cat(index::D(DIM).all(0)));
	  int g1 = gidx(I);

	  if (marked[g1])
	    continue;
	  else
	    marked[g1] = true;

	  for (iterator DI({-1,-1,-1},{1,1,1}); !DI.end(); DI++){
	    auto J = I + DI;
	    if (J(0)>=0 && J(0)<ngrain(0) && J(1)>=0 &&
		J(1)<ngrain(1) && J(2)>=0 && J(2)<ngrain(2) ){
	      int g2 = gidx(J);

	      if (g1>=g2)
		gpairs.insert(std::pair<int,int>(g1,g2));
	      else
		gpairs.insert(std::pair<int,int>(g2,g1));
	    }
	  }
	}
      std::vector<std::pair<int,int> > pairlist(gpairs.begin(),gpairs.end());
      
      std::vector<ngbfinder_thread> threads;
      for (int i =0; i< nparr; i++)
	threads.push_back(ngbfinder_thread(i*pairlist.size()/nparr, (i+1)*pairlist.size()/nparr -1,
					   pairlist, *this));
      for (int i =0; i< nparr; i++) {
	std::thread th(std::ref(threads[i]));
	th.detach();
      }

      while(true){
	bool finished = true;
	//	for (int i =0; i< nparr; i++) std::cout << " " << i << " " << threads[i].imfinished
	//					<< " size " << threads[i].my_table.size()  ;
	//std::cout << "\n";
	for (int i =0; i< nparr; i++)
	  if (!threads[i].imfinished){
	    finished = false;
	    break;
	  }
	if(finished) break;
	std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
      for (int i =0; i< nparr; i++)
	for (const auto x: threads[i].my_table)
	  _add_ngbr(x.first, x.second);
    }

    void build(){
      for (int i=0; i<_table.size(); i++)
	_table[i].clear();
      _table.resize(geom->nat());
      
      if (nparr>1){
	build_parallel();
	return;
      }
      
      if (reference_mode){
	reference_build();
	return;
      }
      
      //std::cout << "grain setup\n";
      grain_setup();
      graining();
      //std::cout << "finished grain setup\n";
      /*
	
	std::cout << typeid(CELL).name() << "\n";
	std::cout <<typeid(periodic_cell<DIM,REAL>).name() << "\n";
	std::cout << typeid(CELL).hash_code() << "\n";
	std::cout <<typeid(periodic_cell<DIM,REAL>).hash_code() << "\n";
      */
      
      if (transl_mode && typeid(*geom->cell) == typeid(periodic_cell<REAL>)){
	translational_build();
	return;
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      
      std::set<std::pair<int,int> > gpairs;
      std::vector<char> marked(ngrain(0)*ngrain(1)*ngrain(2), false);

        /*
      std::cout << "alive1 " << marked.size()<< "\n";
      std::cout << "ngrain= " << ngrain << " Rmin= " <<
      Rmin << " Rmax= " << Rmax << " gs= " << grainsize << "\n";
      */

      for (int at=0; at<geom->nat(); at++)
	if (!geom->shadow(at)){
	  index I = igrain(index({at}).cat(index::D(DIM).all(0)));
	  int g1 = gidx(I);

	  if (marked[g1])
	    continue;
	  else
	    marked[g1] = true;

	  //std::cout << "g1= " << g1 << " I= " << I << "\n";

	  for (iterator DI({-1,-1,-1},{1,1,1}); !DI.end(); DI++){
	    auto J = I + DI;
	    if (J(0)>=0 && J(0)<ngrain(0) && J(1)>=0 &&
		J(1)<ngrain(1) && J(2)>=0 && J(2)<ngrain(2) ){
	      int g2 = gidx(J);

	      //std::cout << g1 << I << g2 << J << "\n";

	      if (g1>=g2)
		gpairs.insert(std::pair<int,int>(g1,g2));
	      else
		gpairs.insert(std::pair<int,int>(g2,g1));
	    }
	  }
	}
      auto t2 = std::chrono::high_resolution_clock::now();
      int t=std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
      std::cout << "Pair list formed  in "<< t << " msec\n";
    /*
      std::sort(gpairs.begin(),gpairs.end());
      auto last = std::unique(gpairs.begin(),gpairs.end());
      gpairs.resize( std::distance(gpairs.begin(),last) );
      */

      for (const auto & gp : gpairs){
	int g1 = gp.first, g2 = gp.second;
	for (int c1 = 0; c1 < grains[g1].size(); c1++)
	  for (int c2 = 0; c2 < (g1==g2 ? c1 : grains[g2].size()); c2++){
	    index at1 = grains[g1][c1];
	    index at2 = grains[g2][c2];

	    bool skip = geom->shadow(at1) || geom->shadow(at2);
	    if (!skip)
	      skip = at1.tail(1) != index::D(DIM).all(0)
		&& at2.tail(1) != index::D(DIM).all(0);
	    if (skip)
	      continue;

	    REAL r = (geom->pos(at1) - geom->pos(at2)).norm();

	    //std::cout << "at1 = " << at1.atom << at1.cell
	    //<< " grain " << I << " at2= "<< at2.atom <<
	    //at2.cell << " grain " << I+DI <<  " DI =" << DI <<"\n";

	    if ( r <= distance(geom->typetable()->type(at1),
			       geom->typetable()->type(at2))){
	      if ( at1.tail(1) == index::D(DIM).all(0) )
		_add_ngbr(at1, at2);
	      if ( at2.tail(1) == index::D(DIM).all(0) )
		_add_ngbr(at2, at1);
	    }

	  }
      }

      //      for (auto gp:gpairs)
      //	std::cout << gp.first << " " << gp.second << "\n";
    }

      typedef  geometry_observer<REAL> DEP;

      void ref_inserted(int at,
                        before_after st,
                        const STRING_EX & a,
                        const vector3<REAL> & r){
        if (st == DEP::after){
            _table.insert(_table.begin()+at, std::vector<index>());
            for (int i=0; i<geom->nat(); i++)
              if (! geom->shadow(i))
                for (iterator j(geom->cell.begin(),
                                geom->cell.end()); j.end(); j++)
                  if ( !(i==at && j==index::D(DIM).all(0)) &&
                       norm(geom->r(at) - geom->r(i, j)) <
                       distance(geom->typetable()->type(at), geom->typetable()->type(i))
                       ){
                      _table[at].push_back(index({i,j}));
                      index iat({at});
                      iat.tail(1) -= j;
                      _table[i].push_back(iat);
                    }
          }
      }

      void ref_added(before_after st,
                     const STRING_EX & a,
                     const vector3<REAL> & r){
        ref_inserted(geom->nat()-1,st,a,r);
      }

      void ref_moved(int at,
                     before_after st,
                     const vector3<REAL> & r){

      }

      void ref_erased(int at,
                      before_after st){

        if (st==DEP::before){
            _table.erase(_table.begin()+at);

            for (int i=geom->nat()-2; i>=0; i--)
              for (int j=n(i)-1; j>=0; j--)
                if (_table[i][j](0)==at){
                    //std::cerr << "erase " << i << " " << _table[i][j] << "\n";
                    _table[i].erase(_table[i].begin()+j);
                  }
                else if (_table[i][j](0)>at){
                    //std::cerr << _table[i][j] << "->";
                    _table[i][j](0)--;
                    //std::cerr << _table[i][j] << "\n";
                  }
          }
      }

      void ref_shaded(  int at,
                        before_after st,
                        bool sh){}

      virtual void added( before_after st,
                          const STRING_EX & a,
                          const vector3<REAL> & r)
      {}

      virtual void inserted(int at,
                            before_after st,
                            const STRING_EX & a,
                            const vector3<REAL> & r)
      {}

      virtual void changed(int at,
                           before_after st,
                           const STRING_EX & a,
                           const vector3<REAL> & r)
      {}

      virtual void erased(int at,
                          before_after st)
      {}

      virtual void shaded(int at,
                          before_after st,
                          bool sh)
      {}

      virtual void reordered(const std::vector<int> &,
                             before_after)
      {}

      //------------------------------------------------------------
      bool operator==(const neighbours_table & t) const{
        if (_table.size() != t._table.size())
          return false;
        for (int i=0; i<_table.size(); i++){
            if (_table[i].size() != t._table[i].size())
              return false;
            for (int j=0; j<_table[i].size(); j++){
                bool found = false;
                for (int k=0; k<_table[i].size(); k++)
                  if (_table[i][j]==t._table[i][k]){
                      found=true;
                      break;
                    }
                if (!found)
                  return false;
              }
          }
        return true;
      }

      bool operator!=(const neighbours_table & t) const {
        return ! (*this == t);
      }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      py::list py_getitem(int i) const{
        if (i<0)
          i += geom -> nat();
        if (i<0 || i>= geom->nat())
          IndexError("Index out of range in ngbr_table::getitem");
        py::list l;
        for(int j=0; j<n(i); j++)
          l.append(table(i, j));
        return l;
      }

      index py_getitem2(py::tuple I) const{
        if (py::len(I)!=2)
          TypeError("neigbour_table::getitem. Expected 2 integers");
        if ( ! py::isinstance<py::int_>(I[0]) || ! py::isinstance<py::int_>(I[1]))
          TypeError("neigbour_table::getitem. Expected 2 integers");

        int i = py::cast<int>(I[0]);
        int j = py::cast<int>(I[1]);

        if (i<0)
          i += geom -> nat();
        if (i<0 || i>= geom->nat())
          IndexError("Index out of range in ngbr_table::getitem");

        if (j<0)
          j += n(i);
        if (j<0 || j>= n(i))
          IndexError("Index out of range in ngbr_table::getitem");

        return table(i,j);
      }

      static void py_export(py::module m, const char * pyname){
        py::class_<neighbours_table<REAL> >(m, pyname)
            .def(py::init< geometry<REAL> &, bonding_table<REAL> & >())
            .def_property("auto_update",
                          &SELF::get_auto_update,
                          &SELF::set_auto_update,
                          "bool auto_update (default: False). If auto_update==True, "
                          "any modifications made to atoms are\nautomatically "
                          "reflected by neighbours_table")

            .def_readwrite("reference_mode",
                           &SELF::reference_mode,
                           "bool reference_mode (default: False). "
                           "If reference_mode==True, simple but very inefficient\n"
                           "algorythm is used for building neighbours_table. For "
                           "small molecules (less then 100 atoms) and for"
                           " debugging purpose")

	  .def_readwrite("nparr",
			 &SELF::nparr,
			 "int nparr (default: 1 - no parallelization)"
			 "number of parallel threads for building ngbr table"
			 )

            .def_property("grain_size",
                          &SELF::get_grain_size,
                          &SELF::set_grain_size,
                          "real grain_size (default: auto). Building and updating"
                          " neighbour_table is done by dividing the\nmolecule "
                          "space into cubic grains of grain_size. The value of "
                          "grain_size\nis normally selected automatically, "
                          "however, you can set it manually")

            .def_readwrite("auto_grainsize",
                           &SELF::auto_grainsize,
                           "bool auto_grainsize (default: True). Whether to "
                           "choose grain_size value automatically")

            .def_readwrite("transl_mode",
                           &SELF::transl_mode)
            .def("build", &SELF::build, "Build the neighbours_table")
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("len", &SELF::n,
                 "Usage: len(int i) - number of neighbours of i-th atom")
            .def("__getitem__", &SELF::py_getitem)
            .def("__getitem__", &SELF::py_getitem2)
            ;
      }

#endif

  };

}

#endif
