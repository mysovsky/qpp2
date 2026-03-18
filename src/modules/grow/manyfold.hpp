#ifndef _QPP_MANYFOLD
#define _QPP_MANYFOLD

#include <mathf/lace3d.hpp>
#include <consts.hpp>
#include <geom/geom.hpp>
//#include <io/qppdata.hpp>
#include <mathf/specfunc.hpp>
#include <cmath>

namespace qpp{

  //#define v2d lace::vector2d<VALTYPE>
  //#define v3d lace::vector3d<VALTYPE>

  //double const default_geomtol = 1e-10;
  inline int const default_maxiter=100;

  template <class REAL>
  class parametric_surface {

  protected:
    
    STRING_EX _name, _error;

    REAL geomtol;
    int maxiter{default_maxiter};

  public:
    typedef parametric_surface<REAL> SELF;
    // limits of the parametric domain    
    virtual vector2<REAL> parmin() const {return {0,0};}
    virtual vector2<REAL> parmax() const {return {0,0};}

    virtual vector3<REAL> map(const vector2<REAL> &v) const {return {0,0,0};}
    // get the point on the surface corresponding to xi and eta parameters values

    vector3<REAL> parm2cart(const REAL & xi, const REAL & eta) const
    // get the point on the surface corresponding to xi and eta parameters values
    {
      return map(vector2<REAL>(xi,eta));
    }

    virtual vector2<REAL> project(const vector3<REAL> & r) const {return {0,0};}
    // get the point on surface closest to r

    virtual REAL dx_dxi(const REAL & xi, const REAL & eta, int i, int j) const{ return 0e0;}
    // derivatives of cartesian coondinates over parameters

    virtual vector3<REAL> normal(const vector2<REAL> &v) const{return {0,0,1};}
    // get the vector normal to the surface

    virtual qpp::matrix2<REAL> gtensor(const vector2<REAL> &v) const{return matrix2<REAL>();}
    // Covariant metric tensor;

    virtual vector2<REAL> ruler(const vector2<REAL> & from, const vector2<REAL> & direction, 
				const REAL & distance) const{return vector2<REAL>(0,0);}
    // find the point at "distance" from the given "from" point in given "direction"

    virtual vector2<REAL> protract(vector2<REAL> a, vector2<REAL> b,REAL distance, REAL angle) const{return vector2<REAL>(0,0);
    }
    // find the point at "distance" from point1 and at "angle" (surface angle - angle between the tangents) 
    // to the point1-point2 direction

    virtual vector2<REAL> triangul(vector2<REAL> a, vector2<REAL> b, REAL distance, REAL angle) const {return vector2<REAL>(0,0);}
    // the same, but the "angle" is between the direction in 3d space
    
    virtual std::vector<vector2<REAL>> triangul2b(vector2<REAL> a1, REAL b1, vector2<REAL> a2, REAL b2)const{
      return std::vector<vector2<REAL>>();
    }
    // find the points situated at distance b1 from point a1 and at distance b2 from the point a2

    virtual REAL surface_angle(vector2<REAL> p1, vector2<REAL> p2, vector2<REAL> p3) const
    {
      return REAL(0);
    }

    virtual STRING_EX name() const // fixme - const
    {
      return _name;
    }
    /*
    virtual int gettype() // fixme - const
    {
      return qppdata_manyfold;
    }
    */
    virtual STRING_EX error()
    {
      return _error;
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    static void py_export(py::module m, const char * pyname) {
      py::class_<SELF>(m,pyname)
	.def("parmin", &SELF::parmin)
	.def("parmax", &SELF::parmax)
	.def("parm2cart", &SELF::parm2cart)
	.def("map", &SELF::map)
	.def("project", &SELF::project)
	.def("gtensor", &SELF::gtensor)
	.def("normal", &SELF::normal)
	.def("ruler", &SELF::ruler)
	.def("protract", &SELF::protract)
	.def("triangul", &SELF::triangul)
	.def("triangul2b", &SELF::triangul2b)
	.def_readwrite("geomtol", &SELF::geomtol)
	;
    }
    
#endif
    
  };


  // --------------------------------------------------------------------
  // The important class of surfaces with diagonal metric tensor
  template <class REAL>
  class orthoparametric_surface : public parametric_surface<REAL>
  {

    using parametric_surface<REAL>::dx_dxi;

  public:

    using parametric_surface<REAL>::geomtol;
    using parametric_surface<REAL>::maxiter;
    using parametric_surface<REAL>::map;
    using parametric_surface<REAL>::project;

    virtual vector2<REAL> gdiag(const vector2<REAL> &v) const {return {1e0,1e0};}
    // Square root of the diagonal components of the metric tensor

    virtual qpp::matrix2<REAL> gtensor(const vector2<REAL> &v) const override
    // Covariant metric tensor;
    {
      vector2<REAL> g = gdiag(v);
      matrix2<REAL> A;
      A(0,0) = g(0)*g(0);
      A(1,1) = g(1)*g(1);
      A(0,1) = 0e0;
      A(1,0) = 0e0;
      return A;
    }

    virtual vector3<REAL> normal(const vector2<REAL> &v) const override
    // get the vector normal to the surface
    {
      vector3<REAL> n1, n2;
      for (int j=0; j<3; j++)
	{
	  n1(j) = dx_dxi(v[0], v[1], j, 0);
	  n2(j) = dx_dxi(v[0], v[1], j, 1);
	}
      return n1.cross(n2) /(n1.norm()*n2.norm());
    }
    
    virtual vector2<REAL> ruler(const vector2<REAL> & from, const vector2<REAL> & direction, 
					  const REAL & distance) const override
    // find the point at "distance" from the given "from" point in given "direction"
    {      
      vector3<REAL> r0 = map(from);
      vector3<REAL> nperp = normal(from), ntan=0;
      for (int i=0; i<3; i++)
	for (int j=0; j<2; j++)
	  ntan(i) += dx_dxi(from(0),from(1),i,j)*direction(j);
      ntan /= ntan.norm();
      
      vector3<REAL> n3 = nperp.cross(ntan);
      
      vector3<REAL> r1, r2, dr;
      vector2<REAL> p1, p2;

      r1 = r0 + ntan*distance;
      p1 = project(r1);
      std::cout << r0 << " " << r1 << " " << p1 << "\n";
      while (true)
	{
	  r2 = map(p1);
	  dr = r2 - r0;
	  dr -= n3*n3.dot(dr);
	  dr *= distance/dr.norm();
	  r2 = r0 + dr;
	  p2 = project(r2);

	  //std::cout << r1 << " " << r2 << " " << p2 << "\n";
	  if ( (r1-r2).norm() < geomtol ) break;

	  r1 = r2;
	  p1 = p2;
	}

      return p2;
    }

    virtual vector2<REAL> protract(vector2<REAL> a, vector2<REAL> b, REAL distance, REAL angle) const override
    // find the point at "distance" from point1 and at "angle" to the point1-point2 direction
    {
      vector3<REAL> A = map(a), B = map(b);
      vector3<REAL> nperp = normal(b), n, nba,nbc;

      nba = A-B;
      nba -= nperp*nba.dot(nperp);
      nba /= nba.norm();

      nbc =  nba*std::cos(angle) + (nperp.cross(nba))*std::sin(angle);
      n =  - nba*std::sin(angle) + (nperp.cross(nba))*std::cos(angle);
      
      //std::cout << "np= " << nperp << " norm = " << norm(nperp) << "\n";
      //std::cout << "nbc= " << nbc << " norm = " << norm(nbc) << "\n";
      //std::cout << "n= " << n << " norm = " << norm(n) << "\n";

      vector3<REAL> C = B + distance*nbc, C2;
      vector2<REAL> c = project(C), c2;
      
      int ii = maxiter;
      
      while (ii-->0)
	{
	  C2 = map(c);
	  C2 -= n*n.dot(C2-B);
	  C2 = B + distance*(C2-B)/(C2-B).norm();
	  c2 = project(C2);

	  //std::cout << norm(C-C2) << "\n";

	  if ((C-C2).norm() < geomtol*pi/180) break;

	  C=C2;
	  c=c2;
	}
      return c2;
    }
    
    virtual vector2<REAL> triangul(vector2<REAL> a, vector2<REAL> b, REAL distance, REAL angle) const override
    {
      vector3<REAL> A = map(a), B = map(b);
      vector3<REAL> nperp = normal(b), n, n0;

      n0 = A-B;
      n0 /= n0.norm();
      n = nperp.cross( n0);
      n /= n.norm();

      REAL ca = std::cos(angle), sa = std::sin(angle), AB = (A-B).norm();
      vector3<REAL> C = B + distance*((A-B)*ca/AB + n*sa), C2;
      vector2<REAL> c = project(C), c2;
      C2 = map(c);
      
      while (true)
	{
	  /*
	  std::cout << A << B << C << "\n";
	  std::cout << distance << " " << norm(B-C)<< "\n";
	  std::cout << angle*180/pi << " " << std::acos(scal(C-B,A-B)/(norm(B-C)*norm(B-A)))*180/pi << "\n";	  
	  */
	  // distance
	  C2 = B + distance*(C2-B)/(C2-B).norm();
	  C2 = map(project(C2));

	  //angle
	  n = n0 - (C2-B)*(C2-B).dot(n0)/(vector3<REAL>(C2-B)).norm2();
	  n /= n.norm();

	  REAL r = (C2-B).norm(),
	    x = (r*ca - (C2-B).dot(n0))/(n.dot(n0)-n.dot(C2-B)*ca/r);
	  C2 += x*n;

	  c2 = project(C2);
	  C2 = map(c2);
	  
	  if ((C-C2).norm() < geomtol*pi/180) break;

	  C=C2;
	  c=c2;
	}
      return c2;
    }

    virtual std::vector<vector2<REAL>> triangul2b(vector2<REAL> a1, REAL b1, vector2<REAL> a2, REAL b2) const override
    {
      vector3<REAL> A1 = map(a1), A2 = map(a2);      
      vector3<REAL> DA = A2 - A1, N1 = DA/DA.norm();

      REAL A = DA.norm();
      REAL x = (b1*b1 - b2*b2 + A*A)/(2*A);
      REAL y = std::sqrt(b1*b1 - x*x);      
      
      std::vector<vector2<REAL>> res;

      if (b1+b2<A)
	{
	  //std::cout << "triangle rule violated: " << b1 << " " << b2 << " " << A << "\n";
	  return res;
	}

      for (int sgn=-1; sgn<=1; sgn+=2)
	{
	  vector2<REAL> r, r1;	  
	  vector3<REAL> R, R1;

	  vector3<REAL> N2 = N1 % normal((a1+a2)/2) * sgn;

	  int it = maxiter;

	  while(it-- > 0)
	    {
	      R = A1 + N1*x + N2*y;
	      r = project(R);
	      
	      R1 = map(r);
	      N2 = R1 - A1;
	      N2 = N2 - N1*N2.dot(N1);
	      N2 /= N2.norm();

	      //std::cout << norm(R-R1) << "\n";
	      if ( (R-R1).norm() < geomtol ) break;
	      
	      r = r1;
	      R = R1;
	    }

	  res.push_back(r);
	}
      return res;
    }

    virtual REAL surface_angle(vector2<REAL> p1, vector2<REAL> p2, vector2<REAL> p3) const override
    {
      vector3<REAL> n = normal(p2);
      vector3<REAL> r1 = map(p1), r2 = map(p2), r3 = map(p3);
      vector3<REAL> n1 = r1-r2, n2 = r3-r2;

      n1 = n1 - n*n.dot(n1);
      n2 = n2 - n*n.dot(n2);

      n1 /= n1.norm();
      n2 /= n2.norm();

      return std::acos(n1.dot(n2));
    }

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    static void py_export(py::module m, const char * pyname) {
      py::class_<orthoparametric_surface<REAL>, parametric_surface<REAL>>(m,pyname);
    }
    
#endif
    
  };

  // -----------------------------------------------------------------------------------
  
  template <class REAL>
  class parametric_sphere : public orthoparametric_surface<REAL>{
    
    REAL R;

    // using typename qpp_object<charT,traits>::string;

  protected:
    using parametric_surface<REAL>::_name;
    
  public:

    using parametric_surface<REAL>::geomtol;
    using parametric_surface<REAL>::maxiter;

    parametric_sphere(const REAL & _R)
    { 
      R=_R;
      geomtol=geometry<REAL>::tol_geom_default;
      maxiter = default_maxiter;
    }

    // limits of the parametric domain    
    virtual vector2<REAL> parmin() const override
    {
      return vector2<REAL>( 0, 0);
    }

    virtual vector2<REAL> parmax() const override
    {
      return vector2<REAL>( pi, 2*pi);
    }

    virtual vector3<REAL> map(const vector2<REAL> &v) const override
    // get the point on the surface corresponding to xi and eta parameters values
    {
      REAL st = std::sin(v(0));
      return vector3<REAL>(R*st*std::cos(v(1)), R*st*std::sin(v(1)), R*std::cos(v(0)));
    }

    virtual vector2<REAL> project(const vector3<REAL> & r) const override
    // get the point on surface closest to r
    {
      std::cout << r.z()/r.norm()<<" " <<  r.x() << " " << r.y()<<" \n";
      return vector2<REAL>(std::acos(r.z()/r.norm()), atanxy(r.x(),r.y()) );
    }

    virtual REAL dx_dxi(const REAL & theta, const REAL & phi, int i, int j) const override
    // derivatives of cartesian coondinates over parameters
    {
      if (j==0)
	{
	  if (i==0)
	    return R*std::cos(theta)*std::cos(phi);
	  else if (i==1)
	    return R*std::cos(theta)*std::sin(phi);
	  else if (i==2)
	    return -R*std::sin(theta);
	}
      else if (j==1)
	{
	  if (i==0)
	    return -R*std::sin(theta)*std::sin(phi);
	  else if (i==1)
	    return R*std::sin(theta)*std::cos(phi);
	  else if (i==2)
	    return REAL(0);
	}
      //fixme - i,j range checking
      IndexError("i,j out of valid range");
      return REAL(0);

    }
    
    vector2<REAL> gdiag(const vector2<REAL> &v) const override
    {	
      return vector2<REAL>(R, R*std::sin(v.x()));
    }
 
    // derived from qpp_object
    virtual STRING_EX category() //const
    {
      return "parametric sphere";
    }

    virtual void write(OSTREAM &os, int offset=0)
    {
      // fixme - write all data 
      for (int k=0; k<offset; k++) os << " ";
      os << "sphere";
      if (_name != "")
	os << " " << _name;
      os << "(" << R << ")";
      if (offset>0)
	os << ";\n";
      
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    static void py_export(py::module m, const char * pyname) {
      py::class_<parametric_sphere<REAL>, orthoparametric_surface<REAL>>(m,pyname)
	.def(py::init<REAL>());
    }    
#endif
    
  };

  // -----------------------------------------------------------------------------------------------

  template <class REAL>
  class parametric_torus : public orthoparametric_surface<REAL>{

    //using typename qpp_object<charT,traits>::string;
    //    using parametric_surface<REAL,charT,traits>::map;

    REAL a,b;
    vector3<REAL> origin, axis;
    matrix3<REAL> Rot;
    
  protected:
    using parametric_surface<REAL>::_name;

  public:

    using parametric_surface<REAL>::geomtol;
    using parametric_surface<REAL>::maxiter;

    parametric_torus( REAL _a, REAL _b, const STRING_EX &__name = "") : origin(0,0,0), axis(0,0,1)
    {
      a = _a; 
      b = _b;
      Rot = REAL(1);
      _name = __name; 
      geomtol=geometry<REAL>::tol_geom_default;
      maxiter = default_maxiter;
      Rot = matrix3<REAL>::unity();
    }

    parametric_torus( REAL _a, REAL _b, vector3<REAL> _origin, 
		      vector3<REAL> _axis, const STRING_EX &__name = "") 
    {
      a = _a; 
      b = _b;
      origin = _origin;
      axis = _axis;
      
      axis = axis/axis.norm();
      vector3<REAL> n( axis.y(), -axis.x(), 0);
      REAL phi = std::acos(axis.z());

      Rot = qpp::RotMtrx(n,phi);
      _name = __name;
      geomtol=geometry<REAL>::tol_geom_default;
      maxiter = default_maxiter;
    }

    virtual vector2<REAL> parmin() const override
    {
      return vector2<REAL>(REAL(0),REAL(0));
    }
    virtual vector2<REAL> parmax() const override
    {
      return vector2<REAL>(2*pi,2*pi);
    }

    virtual vector3<REAL> map(const vector2<REAL> &v) const override
    // get the point on the surface corresponding to xi and eta parameters values
    {
      vector3<REAL> r;
      REAL theta = v[0], phi = v[1];
      r.z() = b*std::sin(theta);
      r.x() = a - b*std::cos(theta);
      r.y() = r.x()*std::sin(phi);
      r.x() *= std::cos(phi);

      return origin + Rot*r;
    }

    virtual vector2<REAL> project(const vector3<REAL> & r) const override
    // get the point on surface closest to r
    {
      vector3<REAL> r1 = Rot.inverse()*(r-origin);
      std::cout << r << " " << r1 << "\n";
      REAL phi = atanxy(r1.x(),r1.y());
      REAL theta = atanxy(a-std::sqrt(r1.x()*r1.x()+r1.y()*r1.y()), r1.z());
      return vector2<REAL>(theta,phi);
    }

    virtual REAL dx_dxi(const REAL & theta, const REAL & phi, int i, int j) const override
    // derivatives of cartesian coondinates over parameters
    {
      if (j==0)
	{
	  if (i==0)
	    return  b*std::sin(theta)*std::cos(phi);
	  else if (i==1)
	    return  b*std::sin(theta)*std::sin(phi);
	  else if (i==2)
	    return  b*std::cos(theta);
	}
      else if (j==1)
	{
	  if (i==0)
	    return -(a - b*std::cos(theta))*std::sin(phi);
	  else if (i==1)
	    return  (a - b*std::cos(theta))*std::cos(phi);
	  else if (i==2)
	    return REAL(0);
	}
      IndexError("i,j out of valid range");
      return REAL(0);
      //fixme - i,j range checking      
    }

    vector2<REAL> gdiag(const vector2<REAL> &v) const override
    {	
      return vector2<REAL>(b, a - b*std::cos(v.x()));
    }
 
    // derived from qpp_object
    virtual STRING_EX category() const //const
    {
      return "parametric torus";
    }

    virtual void write(OSTREAM &os, int offset=0)
    {
      // fixme - write all data 
      for (int k=0; k<offset; k++) os << " ";
      os << "torus";
      if (_name != "")
	os << " " << _name;
      os << "(" << a << "," << b << ")";
      if (offset>0)
	os << ";\n";
      
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
    static void py_export(py::module m, const char * pyname) {
      py::class_<parametric_torus<REAL>, orthoparametric_surface<REAL>>(m,pyname)
	.def(py::init< REAL, REAL, const STRING_EX &>())
	.def(py::init< REAL, REAL, vector3<REAL>, 
	     vector3<REAL>, const STRING_EX &>());
	;
    }    
#endif    

  };
  
  // ----------------------------------------------------------------------------

  template <class REAL>
  class parametric_plane : public orthoparametric_surface<REAL>{

    // using typename qpp_object<charT,traits>::string;

    vector3<REAL> origin, axis;
    matrix3<REAL> Rot;
    
  protected:
    using parametric_surface<REAL>::_name;

  public:

    using parametric_surface<REAL>::geomtol;
    using parametric_surface<REAL>::maxiter;

    parametric_plane(STRING_EX __name = "") : origin(0,0,0), axis(0,0,1)
    {
      Rot = REAL(1);
      _name = __name;
      geomtol=geometry<REAL>::to_geom_default;
      maxiter = default_maxiter;
    }

    parametric_plane( vector3<REAL> _origin, vector3<REAL> _axis, STRING_EX __name = "") 
    {
      origin = _origin;
      axis = _axis;
      
      axis = axis/norm(axis);
      vector3<REAL> n( axis.y(), -axis.x(), 0);
      REAL phi = std::acos(axis.z());

      Rot = qpp::RotMtrx(n,phi);
      _name = __name;
      geomtol=geometry<REAL>::to_geom_default;
      maxiter = default_maxiter;
    }

    virtual vector2<REAL> parmin() const
    {
      return vector2<REAL>(REAL(0),REAL(0));
    }
    virtual vector2<REAL> parmax() const
    {
      return vector2<REAL>(REAL(0),REAL(0));
    }

    virtual vector3<REAL> map(const vector2<REAL> &v) const
    // get the point on the surface corresponding to xi and eta parameters values
    {
      vector3<REAL> r;
      r.x() = v.x;
      r.y() = v.y;
      r.z() = REAL(0);

      return origin + Rot*r;
    }

    virtual vector2<REAL> project(const vector3<REAL> & r) const override
    // get the point on surface closest to r
    {
      vector3<REAL> r1 = invert(Rot)*(r-origin);
      return vector2<REAL>(r1.x(),r1.y());
    }

    virtual REAL dx_dxi(const REAL & theta, const REAL & phi, int i, int j) const
    // derivatives of cartesian coondinates over parameters
    {
      return Rot(i,j);
      //fixme - i,j range checking      
    }

    vector2<REAL> gdiag(const vector2<REAL> &v) const
    {	
      return vector2<REAL>(REAL(1),REAL(1));
    }
 
    // derived from qpp_object
    virtual STRING_EX category() //const
    {
      return "plane";
    }

    virtual void write(OSTREAM &os, int offset=0)
    {
      // fixme - write all data 
      for (int k=0; k<offset; k++) os << " ";
      os << "plane";
      if (_name != "")
	os << " " << _name;
      if (offset>0)
	os << ";\n";
      
    }

  };

  /*
  template <class REAL, class charT=std::string::value_type , class traits = std::char_traits<charT> >
  class parametric_plane : public parametric_surface<REAL,charT,traits>{

    using typename qpp_object<charT,traits>::string;
    //    using parametric_surface<REAL,charT,traits>::parm2xyz;

    vector3<REAL> origin,axis1,axis2;
    
  protected:
    using parametric_surface<REAL,charT,traits>::_name;

  public:

    parametric_plane( string __name = "") : origin(0,0,0), axis1(1,0,0), axis2(0,1,0)
    {
      _name = __name;
    }

    parametric_plane( vector3<REAL> _origin, vector3<REAL> _axis1,  
		      vector3<REAL> _axis2, string __name = "")
    {
      origin = _origin;
      axis1 = _axis1;
      axis2 = _axis2;
      _name = __name;
    }

    virtual vector3<REAL> parm2xyz(REAL x, REAL y)
    {
      return origin + axis1*x + axis2*y;
    }

    virtual vector2<REAL> xyz2parm(vector3<REAL> r)
    {}
    // get the point on surface closest to r

    virtual vector2<REAL> ruler(vector2<REAL> from, vector2<REAL> direction, 
				    REAL distance)
    {
      return from + distance*direction/norm(direction);
    }

    virtual vector2<REAL> protractor(vector2<REAL> point1, vector2<REAL> point2, 
					 REAL distance, REAL alpha)
    {
      vector2<REAL> delta0 = point1-point2;
      REAL beta = atanxy(delta0.x, delta0.y);
      vector2<REAL> n( std::cos(alpha+beta), std::sin(alpha+beta) );
      return ruler(point2,n,distance);
    }
    
    virtual REAL xi_min()
    {
      return 0;
    }

    virtual REAL xi_max()
    {
      return 0;
    }

    virtual REAL eta_min()
    {
      return 0;
    }

    virtual REAL eta_max()
    {
      return 0;
    }

    virtual lace::matrix3d<REAL> gtensor(REAL theta, REAL phi)
    {

    }


    // derived from qpp_object
    virtual string category()
    {
      return "plane";
    }

    virtual void write(std::basic_ostream<charT,traits> &os, int offset=0)
    {
      // fixme - write all data 
      for (int k=0; k<offset; k++) os << " ";
      os << "plane";
      if (_name != "")
	os << " " << _name;
      os << "(" << ")";
      if (offset>0)
	os << ";\n";
      
    }

  };
  */


};

#endif
