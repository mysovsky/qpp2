#ifndef QPP_TRANSFORM_H
#define QPP_TRANSFORM_H
#include <mathf/lace3d.hpp>
#include <symm/cell.hpp>
#include <consts.hpp>
#include <vector>
#include <memory>

namespace qpp{

  /*! \brief rotrans is the general rotational and translational symmetry operation turning the vector r into R*r+T, where R is rotation (and probably reflection) matrix and T is a translation vector. Thus, rotrans class implements all possible operations of space groups.
    @param REAL either float or double depending of what precision is necessary
    @param BOUND boolean parameter specifying if the rotrans should be "bound" to certain periodic cell. The need for "bound rotrans" is the following: a symmetry group consisting of rotrans operations is, generally speaking, infinite. "Bound rotrans" makes it finite by reducing the symmetry operation result (a vector, a point) to 3d-periodic cell with given translation vectors. Such periodic cell must be the same for all ""bound rotranses" comprising the group.
   */
  template <typename ELEM> using matrix4 = Eigen::Matrix<ELEM, 4, 4>;

  template<class REAL>
  struct rotrans {
    typedef periodic_cell<REAL>  BOUNDARY;
    
    //!\brief the threshhold to consider two translation vectors equal
    //   altref<REAL> tol_transl;
    
    //!\brief the threshhold to consider two rotation matricies equal
    //altref<REAL> tol_rot;
    
    static rotrans<REAL> unity()
    {
      return rotrans<REAL>(vector3<REAL>::Zero(), matrix3<REAL>::unity());
    }
    
    //!\brief the translation vector
    vector3<REAL> T;
    
    //!\brief the rotation matrix
    matrix3<REAL> R;
    
    //!\brief the periodic cell to which the bound rotrans is bound
    std::shared_ptr<BOUNDARY> cell{nullptr};

    //!\brief empty constructor creates unity operation
    rotrans():      cell (std::make_shared<periodic_cell<REAL>>(0)){
		    //tol_transl(cell->tol_transl),
		    //tol_rot(cell->tol_rot){
      T = vector3<REAL>(0);
      R = matrix3<REAL>::unity();

    }

      //!\brief copy constructor
    rotrans(const rotrans<REAL> & a): cell(a.cell){
					//tol_transl(cell->tol_transl),
					//tol_rot(cell->tol_rot){
      T = a.T;
      R = a.R;
      cell = a.cell;
    }

    //!\brief creates translation operation with _T vector
    rotrans(const vector3<REAL> &_T, std::shared_ptr<BOUNDARY>  _cell = nullptr):
      cell(_cell){//,tol_transl(_cell->tol_transl)
      //tol_rot(_cell->tol_rot){
      T=_T;
      R = matrix3<REAL>::unity();
      if (_cell==nullptr)
	cell = std::make_shared<periodic_cell<REAL>>(0);
      else
	cell = _cell;
    }

    //!\brief creates rotation operation with _R matrix
    rotrans(const matrix3<REAL> & _R, std::shared_ptr<BOUNDARY>  _cell = nullptr):
      cell(_cell==nullptr ? std::make_shared<periodic_cell<REAL>>(periodic_cell<REAL>(0)): _cell)
      //tol_transl(_cell->tol_transl),
      //tol_rot(_cell->tol_rot)
    {
      T = vector3<REAL>(0);
      R = _R;/*
      if (_cell==nullptr)
	cell = std::make_shared<periodic_cell<REAL>>(periodic_cell<REAL>(0));
      else
      cell = _cell;*/
    }
    
    //!\brief creates rotation-translation operation with vector _T and matrix _R
    rotrans(const vector3<REAL> &_T, const matrix3<REAL> &_R,
	    std::shared_ptr<BOUNDARY> _cell = nullptr)/*
      cell(_cell==nullptr ? std::make_shared<periodic_cell<REAL>>(periodic_cell<REAL>(0)): _cell),
      tol_transl(cell->tol_transl),
      tol_rot(cell->tol_rot)*/
    {
      T = _T;
      R = _R;
      if (_cell==nullptr)
	{
	  cell = std::make_shared<periodic_cell<REAL>>(periodic_cell<REAL>(0));
	  //std::cout << "new cellcreated\n";
	}
      else
	{
	  cell = _cell;
	  //std::cout << "cell :" << cell->v[0] << cell->v[1] << cell->v[2] << "\n";
	}
      //tol_transl.bind(cell->tol_transl);
      //tol_rot.bind(cell->tol_rot);
      //std::cout << tol_transl() << " " << tol_rot<< "\n";
	     
    }
    void set_cell( std::shared_ptr<BOUNDARY> _cell)
    {
      cell = _cell;
    }

    //!\brief Multiplication of two rotrans operations
    inline rotrans<REAL> operator*(const rotrans<REAL> & b) const {
      vector3<REAL> t = T + R*b.T;
      vector3<REAL> f =  cell -> cart2frac(t);
      for (int d=0; d<cell->DIM; d++) {
	f(d) -= floor(f(d));
	if (std::abs(f(d)-REAL(1)) < *cell->tol_transl)
	  f(d)=REAL(0);
      }
      t = cell -> frac2cart(f);
      return rotrans<REAL>(t, R*b.R, cell);
    }
  

      /*!\brief Comparison of two rotrans operations. They are considered equal
      if their translations differ by less than tol_trans and their
      rotation matricies differ by less than tol_rot.
    */
    inline bool operator==(const rotrans<REAL> & b) const {
      std::cout << "rotrans== tol:"<< cell->tol_transl<< "DT:"
		<< (T - b.T).norm()<< "DR:"<< (R - b.R).norm()<< "\n";
      return (T - b.T).norm()<= *cell->tol_transl &&
	(R - b.R).norm() <= *cell->tol_rot;
    }
  /*
        else {
            if ( (R - b.R).norm() > tol_rot) return false;
            vector3<REAL> f =  cell -> cart2frac(T - b.T);

            //debug
            //std::cout << "rotrans::== f= " << f;

            for (int d=0; d<cell->DIM; d++) f(d) -= floor(f(d)+.5);
            f = cell->frac2cart(f);

            //std::cout << " t= " << f << "\n";

            return (f).norm() <= tol_trans;
          }
      }
  */

      //!\brief Inequality operator. Simply !(a==b).
      inline bool operator!=(const rotrans<REAL> & b) const {
        return !(*this == b);
      }

      //!\brief Rotrans times vector multiplication means that rotrans
      //!  acts on this vector
      inline vector3<REAL> operator*(const vector3<REAL> & v) const {
        vector3<REAL> res = T+R*v;
	res = cell -> reduce(res);
        return res;
      }

      inline rotrans<REAL> pow(REAL fn) const {
        int n = floor(fn);
        // fixme - inefficient
        rotrans<REAL> A = rotrans<REAL>::unity();
        rotrans<REAL> C = *this;
        A.cell = C.cell;
        if (n>0){
            while (n-- > 0)
              A = (C)*A;
          }
        else if (n<0){
            rotrans<REAL>  B = C.inverse();
            while (n++ < 0)
              A = (B)*A;
          }
        return A;
      }

      inline rotrans<REAL> inverse() const {
        matrix3<REAL> A = ((*this).R).inverse();
        vector3<REAL> t = - A*(*this).T;
        return rotrans<REAL>(t, A, (*this).cell);
      }

      //!\ Rotrans oputput in qpp format
      virtual void write(std::basic_ostream<CHAR_EX,TRAITS> &os, int offset=0) const{
        for (int k=0; k<offset; k++) os << " ";
	os << "rotrans(";
        os << T << "," << R << ")";
      }


#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

      inline vector3<REAL> py_mulv(const vector3<REAL> & v) const
      {return (*this)*v; }
      inline rotrans<REAL> py_mulr(const rotrans<REAL> & b) const {return (*this)*b; }

    REAL get_tol_transl(){
      std::cout << "get_tol_transl\n";
      //std::cout << "isnull=" << tol_transl.isnull()<< "\n";
      if (cell==nullptr)
	{
	  std::cout << "nullptr\n";
	return 0e0;
	}
      else{
	std::cout << cell->tol_transl<< "\n";
	return *cell->tol_transl;}
    }
    
    REAL get_tol_rot(){
      if (cell==nullptr)
	return 0e0;
      else
	return *(cell->tol_rot);}
    
    void set_tol_transl(REAL v){
      if (cell==nullptr)
	return;
      else
	*(cell->tol_transl)=v;}

    void set_tol_rot(REAL v){
      if (cell==nullptr)
	return;
      else
	*(cell->tol_rot)=v;}
#endif

  };

  template<typename _CharT, class _Traits, class VALTYPE>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT, _Traits>& __os,
             const rotrans<VALTYPE> & r){
    r.write(__os);
    return __os;
  }

/*

  template<class REAL, bool BOUND>
  REAL rotrans<REAL>::tol_trans = vector3<REAL>::tol_equiv;

  template<class REAL, bool BOUND>
  REAL rotrans<REAL>::tol_rot = matrix3<REAL>::tol_equiv;
*/
  /* 
  template<class REAL>
  rotrans<REAL> rotrans<REAL>::unity= rotrans<REAL>(vector3<REAL>::Zero(), matrix3<REAL>::unity);
  }
  */
  template<class REAL>
  matrix4<REAL> rotrans4d(const rotrans<REAL> & R)
  {
    matrix4<REAL> res = matrix4<REAL>::Identity();
    res.block(0,0,3,3) = R.R;
    res.block(0,3,3,1) = R.T;
    return res;
  }
  
  template <class REAL>
  rotrans<REAL> rotrans_shift(const rotrans<REAL> & r, const index & I)
  {
    auto & cl = *r.cell;
    return rotrans<REAL>(vector3<REAL>(r.T + cl(0)*I(0) + cl(1)*I(1) + cl(2)*I(2)), r.R, r.cell);
  }

  template <class REAL>
  rotrans<REAL> rotrans2frac(const rotrans<REAL> & r){
    const auto & cl = *r.cell;
    matrix3<REAL> A = {cl(0),cl(1),cl(2)}, B;
    B = A.transpose();
    A = B;
    B = A.inverse();
    return rotrans<REAL>(B*r.T, B*r.R*A);
  }

  template <class REAL>
  rotrans<REAL> rotrans2cart(const rotrans<REAL> & f, periodic_cell<REAL> & cl){
    matrix3<REAL> A = {cl(0),cl(1),cl(2)}, B;
    B = A.transpose();
    A = B;
    B = A.inverse();
    return rotrans<REAL>(A*f.T, A*f.R*B, &cl);
  }
  
  //!\brief Inverse of rotrans R, P=R^(-1). Means R*P==1.
  //  template<class REAL, bool BOUND>
  //  rotrans<REAL,BOUND> invert(const rotrans<REAL,BOUND> & R){
  //    matrix3<REAL> A = (R.R).inverse();
  //    vector3<REAL> t = - A*R.T;
  //    return rotrans<REAL,BOUND>(t, A, R.cell);
  //  }

  //!\brief Power n of rotrans R
  //template<class REAL, bool BOUND>


  // ----------------------------------------------------------------
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

  //  template<class REAL,bool BOUND>
  //  inline rotrans<REAL,BOUND> py_invert_rt(const rotrans<REAL,BOUND> & R) {
  //    return invert(R);}

  //  template<class REAL,bool BOUND>
  //  inline rotrans<REAL,BOUND> py_pow_rt(const rotrans<REAL,BOUND> & R, int n) {
  //    return pow(R,n);}

#endif



  // ----------------------------------------------------------------
  /*
  template<class REAL>
  struct translation
  {
    //typedef typename BOUNDARY::real real;
    static REAL tol_trans;

    static translation<REAL> unity;

    vector3<REAL> T;

    translation()
    { T = 0e0; }

    translation( const translation<REAL> & a)
    {
      T = a.T;
    }

    translation(const vector3<REAL> & t)
    {
      T = t;
    }

    inline translation<REAL> operator*(const translation<REAL> & b) const
    {
      vector3<REAL> res = T+b.T;
      return translation<REAL>(res);
    }

    inline bool operator==(const translation<REAL> & b) const
    {
      return norm(T - b.T) <= tol_trans;
    }

    inline bool operator!=(const translation<REAL> & b) const
    {
      return norm(T - b.T) > tol_trans;
    }

    inline vector3<REAL> operator*(const vector3<REAL> & v) const
    {
      vector3<REAL> res = T+v;
      return res;
    }

    inline operator vector3<REAL>() const
    {
      return T;
    }

  };

  template<class REAL>
  REAL translation<REAL>::tol_trans = 1e-10;

  template<class REAL>
  translation<REAL>  translation<REAL>::unity(vector3<REAL>(0e0));

  template<class REAL>
  translation<REAL> invert(const translation<REAL> & R)
  {
    return translation<REAL>(-R.T);
  }

  template<class REAL>
  translation<REAL> pow(const translation<REAL> & R, int n)
  {
    return translation<REAL>(R.T*n);
  }
  */

}

#endif
