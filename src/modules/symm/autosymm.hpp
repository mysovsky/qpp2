#ifndef QPP_AUTOSYMM_H
#define QPP_AUTOSYMM_H

//#include <geom/shape.hpp>
#include <symm/cell.hpp>
#include <symm/gen_cell.hpp>
#include <symm/transform.hpp>
#include <symm/group_theory.hpp>
#include <symm/point_groups.hpp>
#include <symm/permut.hpp>
//#include <symm/subspace.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <optional>

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pyqpp/py_indexed_property.hpp>
#include <geom/ngbr.hpp>
#include <geom/shape.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif


namespace qpp {
  
  template<class REAL>
  bool has_symmetry (geometry<REAL> & geom, periodic_cell<REAL> & symm,
                     REAL R = geometry<REAL>::tol_geom_default) {

    auto cell = geom.cell;
    geom.cell = std::shared_ptr<periodic_cell<REAL>>(&symm);
    //geom.set_DIM(symm.DIM);

    bonding_table<REAL> b;
    b.default_distance = R;
    neighbours_table<REAL> ngbr(geom,b);

    ngbr.build();

    bool res = true;
    int Ng = 1;
    for (int d = 0; d<symm.DIM; d++)
      Ng *= symm.end()(d) - symm.begin()(d) + 1;

    //std::cout << "Group orders: " << symm.begin() << symm.end() << "\n";
    //std::cout << "Group size = " << Ng << "\n";

    for (int i=0; i<geom.nat(); i++) {
      int found = 0;

      //std::cout << "atom " << i;

      for (int j=0; j<ngbr.n(i); j++) {
	index J = ngbr(i,j);
	//std::cout << " " << J;

	if ( geom.atom(i) == geom.atom(J) &&
	     J != index::D(symm.DIM).atom(J(0))) {
	  found++;
	  //std::cout << "EQ";
	}
	if (found == Ng-1)
	  break;
      }

      //std::cout << "\nnfound= " << found << "\n";

      if ( found != Ng-1){
	res = false;
	//break;
      }
    }

    geom.cell = cell;
    //geom.set_DIM(cell.DIM);
    return res;
  }

  // -------------------------------------------------------------
  /*
    template <class REAL>
    void add_subspace(std::vector<subspace3<REAL> > & subspaces,
    std::vector<std::vector<rotrans<REAL,false> > > &elements,
    const subspace3<REAL> &s,
    const std::vector<rotrans<REAL,false> > &g){
    int i=0;
    while (i<subspaces.size() && subspaces[i] != s) i++;


    if (i<subspaces.size())
    for (const auto & gg : g){
    if ( std::find(elements[i].begin(),elements[i].end(),gg) ==
    elements[i].end() )
    elements[i].push_back(gg);

    std::cout << "adding element to " << i <<
    " sub= (" << s.dim << s.point << s.axis  << ") g= " << gg << "\n";

    }
    else if (s.dim > -1){
    subspaces.push_back(s);
    elements.push_back(g);

    std::cout << "adding subspace (" << s.dim << s.point << s.axis  << ") g= ";
    for (const auto & gg : g)   std::cout << gg;
    std::cout << "\n";
    }
    }
  */
  enum{
    bravais_triclinic = 1,
    bravais_monoclinic = 2,
    bravais_monoclinic_bas = 3,
    bravais_orthorhombic_prim = 4,
    bravais_orthorhombic_base = 5,
    bravais_orthorhombic_body = 6,
    bravais_orthorhombic_face = 7,
    bravais_tetragonal_prim = 8,
    bravais_tetragonal_body = 9,
    bravais_hexagonal = 10,
    bravais_trigonal = 11,
    bravais_cubic_prim = 12,
    bravais_fcc = 13,
    bravais_bcc = 14
  };


  /*! \brief Find all point symmetry operations of 3 lattice vectors
   *  comprising a periodic cell
   @param[out] G the resulting point symmetry group in array form
   @param[in] cell the periodic cell
   @param[in] R tolerance radius
  */
  template<class REAL>
  void bravais_point_group(array_group<matrix3<REAL> > & G,
                           const periodic_cell<REAL> & cell,
                           REAL R = geometry<REAL >::tol_geom_default){
    if (cell.DIM != 3)
      IndexError("bravais_point_group:: works only for 3d-periodic lattice");

    REAL amax = cell(0).norm();
    if (cell(1).norm() > amax) amax = cell(1).norm();
    if (cell(2).norm() > amax) amax = cell(2).norm();

    geometry<REAL> points(0);
    shape_sphere<REAL> S(amax+R);

    vector3<REAL> fmax = S.fmax(cell);
    int fx = int(fmax(0))+1,
      fy = int(fmax(1))+1,
      fz = int(fmax(2))+1;
    for (iterator I({-fx,-fy,-fz}, {fx,fy,fz}); !I.end(); I++){
      vector3<REAL> r = cell.transform({0.,0.,0.},I);
      if (S.within(r))
	points.add(std::string("point"), r);
    }

    vector3<REAL> new_centre;
    find_point_symm(G, points, new_centre, R);
  }

  // ------------------------------------------------------------------------------------------

  template<class REAL>
  void find_translations(std::vector<vector3<REAL> > & transl,
                         std::vector<permutation> & perm,
                         geometry<REAL > & g1,
                         geometry<REAL> & g2,
                         const periodic_cell<REAL> &cell,
                         REAL R = geometry<REAL >
                         ::tol_geom_default){
    transl.clear();
    g1.build_typetable();
    g2.build_typetable();
    if (g1.nat() != g2.nat() || g1.typetable()->n_types() != g2.typetable()->n_types())
      return;
    int nt = g1.typetable()->n_types();
    std::vector<std::vector<int> > t1(nt),t2(nt);
    for (int i=0; i<g1.nat(); i++){
      t1[g1.typetable()->type(i)].push_back(i);
      int t = g1.typetable()->type_of_atom(g2.atom(i));
      if (t==-1) return;
      t2[t].push_back(i);
    }

    // debug
    /*
      std::cout << "find_translations\n";
      for (int t=0; t<t1.size(); t++)
      {
      std::cout << "(" << g1.atom_of_type(t);
      for (int i=0; i<t1[t].size(); i++) std::cout << "," << t1[t][i];
      std::cout << ")";
      }
      std::cout << "\n";
      for (int t=0; t<t2.size(); t++)
      {
      std::cout << "(" << g1.atom_of_type(t);
      for (int i=0; i<t2[t].size(); i++) std::cout << "," << t2[t][i];
      std::cout << ")";
      }
      std::cout << "\n";
    */

    if (! std::equal(t1.begin(), t1.end(), t2.begin()))
      return;

    //reduce to uc
    for (int i=0; i<g1.nat(); i++){
      if (g1.frac)
	for (int j=0; j<3; j++)
	  g1.coord(i)(j) -= floor(g1.coord(i)(j));
      else
	g1.coord(i) = cell.reduce(g1.coord(i));

      if (g2.frac)
	for (int j=0; j<3; j++)
	  g2.coord(i)(j) -= floor(g2.coord(i)(j));
      else
	g2.coord(i) = cell.reduce(g2.coord(i));
    }

    int t=0;
    for (int i=0; i<t1.size(); i++)
      if (t1[t].size()>t1[i].size())
        t=i;

    //std::cout << "t= " << t << "\n";

    for (int i=0; i<t1[t].size(); i++){
      vector3<REAL> v = g2.pos(t2[t][i]) - g1.pos(t1[t][0]), vs=v;
      if (g1.frac) vs = cell.cart2frac(v);

      geometry<REAL> g(g1);
      for (int j=0; j<g.nat(); j++)
	g.coord(j) += vs;


      //std::cout << i << " v= " << v << "\n";


      bool is_transl = true;
      permutation P(g.nat());

      for (int j=0; j<g.nat(); j++){
	bool found = false;
	for (iterator J({-1,-1,-1},{1,1,1}); !J.end(); J++)
	  for (int k=0; k<g.nat(); k++)
	    if ((g2.pos(k) - g.pos(j,J)).norm() < 2*R && g2.atom(k)==g.atom(j)){
	      //std::cout << j << g.coord(j);

	      found = true;
	      if (g.frac){
		g.coord(j)(0) += J(0);
		g.coord(j)(1) += J(1);
		g.coord(j)(2) += J(2);
	      }
	      else
		g.coord(j) += cell(0)*J(0)+cell(1)*J(1)+cell(2)*J(2);

	      P[j] = k;
	      //std::cout << g.coord(j) << "\n";

	      goto FOUND;
	    }
      FOUND:
	if (!found){
	  is_transl = false;
	  break;
	}
      }
      if (is_transl){
	vector3<REAL> dv = vector3<REAL>::Zero();
	for (int j=0; j < g.nat(); j++)
	  dv += g2.pos(j) - g.pos(j);
	dv /= g.nat();
	transl.push_back(v+dv);
	perm.push_back(P);
      }
    }

    //debug
    /*
      std::cout << "alive after all!\n";
      for (int i=0; i<transl.size(); i++)
      {
      std::cout << i << transl[i];
      std::cout << perm[i].to_string() << "\n";
      }
    */

  }

    
  template<class REAL>
  rotrans<REAL> rotrans_sub(const rotrans<REAL> & R1, const rotrans<REAL> & R2)
  {
    auto cell = R1.cell;
    vector3<REAL> df = cell->cart2frac(R1.T - R2.T);
    for (int i=0; i<3; i++)
      df[i] -= round(df[i]);
    return rotrans<REAL>(cell->frac2cart(df), R1.R - R2.R, cell);
  }

    
  template<class REAL>
  void rotrans_diff(REAL & rdiff, REAL & tdiff,
		    const rotrans<REAL> & R1, const rotrans<REAL> & R2)
  {
    rotrans<REAL> dR = rotrans_sub( R1, R2);
    rdiff = dR.R.norm();
    tdiff = dR.T.norm();
  }

  
  template<class REAL>
  void fix4_cryst_group(array_group<rotrans<REAL> > & G, const static_table<int> & M){
    int N = G.size();
    static_table<rotrans<REAL> > F(N,N);

    //std::cout << "fix4: alive 1\n";

    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++){
	F(i,j) = rotrans_sub(G[i]*G[j], G[M(i,j)]);
          //std::cout << "Fij " << i << " " << j << " " << F(i,j) << "\n";
        }

    std::vector<rotrans<REAL> > Ginv(N);

    for (int j=0; j<N; j++){
      matrix3<REAL> S = G[j].R.inverse();
      Ginv[j] = rotrans<REAL>(-1e0*S*G[j].T, S, G[0].cell);
        //std::cout << j << " inv " << Ginv[j] << "\n";
    }

    //std::cout << "fix4: alive 2\n";

    std::vector<rotrans<REAL> > Di(N);
    for (int i=0; i<N; i++){

      //std::cout << "Di " << i << "\n";
      
      rotrans<REAL> D(vector3<REAL>(0),G[0].cell);
      for (int j=0; j<N; j++) {
	rotrans<REAL> DD;
	DD = Ginv[j]*F(j,i);
	D.R += DD.R;
	D.T += DD.T;
	DD = F(i,j)*Ginv[j];
	D.R += DD.R;
	D.T += DD.T;	
        //std::cout << i << " Di " << Di[i] << "\n";
      }
      D.R /= 2*N;
      D.T /= 2*N;
      Di[i] = D;

    }

    //std::cout << "fix4: alive 3\n";

    for (int i=0; i<N; i++){

        //std::cout << "unitarize " << i << "\n";

        G[i].R = G[i].R - Di[i].R;
        G[i].T = G[i].T - Di[i].T;

        //std::cout << "new Gi " << i << " " << G[i] << "\n";

        unitarize(G[i].R);

        //std::cout << "unitarized "  << G[i] << "\n";
      }
    //std::cout << "fix4: alive 4\n";
  }

  
  template<class REAL>
  void reconstruct_cryst_group(array_group<rotrans<REAL> > & G,
			       const static_table<int> & M){
    auto cell = G[0].cell;

    int maxit = 100;
    int N = G.size();
    int it = 0;
    REAL reps = (*cell->tol_rot)*N;
    REAL teps = (*cell->tol_transl)*N;
    REAL rdiff, tdiff, rerr, terr;

    std::cout << "reconstruct_cryst_group\n";

    while (true) {

      rerr = REAL(0);
      terr = REAL(0);
      for (int i=0; i<N; i++)
	for (int j=0; j<N; j++){
	  rotrans_diff(rdiff, tdiff, G[i]*G[j], G[M(i,j)]);
	  rerr += rdiff*rdiff;
	  terr += tdiff*tdiff;
	}
      rerr = std::sqrt(rerr);
      terr = std::sqrt(terr);

      std::cout << "fix_cryst_group: iteration = " << it << " rterror = " << rerr << " " << terr << "\n";

      if (rerr < reps && terr < teps) break;      
      fix4_cryst_group(G,M);
      
      //if (++it > maxit) OverflowError("Too many iterations in fix_point_group");
      if (++it > maxit) {
	std::cout << "Too many iterations in fix_point_group\n";
	break;
      }
    }
  }

  
  /*! \brief Find the translation that turns 3D-periodic geometry g1 into g2
    @param[out] t the translation vectors
    @param[in] g1 initial geometry
    @param[in] g2 final geometry
    @param[in] cell the periodic cell of DIM==3
    @param[in] R the tolerance radius
  */
  template<class REAL>
  void find_translations(std::vector<vector3<REAL> > & transl,
                         geometry<REAL > & g1,
                         geometry<REAL > & g2,
                         const periodic_cell<REAL> &cell,
                         REAL R = geometry<REAL>
                         ::tol_geom_default){
    std::vector<permutation> perm;
    find_translations(transl,perm,g1,g2,cell,R);
  }

  /*! \brief Find the crystalline symmetry group
    @param[out] G - the crystalline group of bound rotranslational
    operations in array form.
    Bound rotranslational operations are used to make the group finite
    @param[in] geom - the geometry of unit cell together with lattice vectors,
    which should be
    stored in geom.cell object. This geometry must be 3D - periodic
    @param[in] R( - the tolerance radius. Symmetry operation is considered
    valid, if the displacement of atom due to
    this operation is less than R
  */
  template<class REAL>
  void find_cryst_symm(array_group<rotrans<REAL> > & G,
		       std::vector<permutation> & P,
                       geometry<REAL> & geom,
                       REAL R = geometry<REAL >
                       ::tol_geom_default){
    //spgw_get_symmetry(G,geom,R);

    array_group<matrix3<REAL> > B;
    bravais_point_group(B,*geom.cell.get(),R);
    G.group.clear();
    G.group.push_back(rotrans<REAL>(matrix3<REAL>::unity(),  geom.cell));

    //std::vector<permutation> P;
    P.clear();
    P.push_back(permutation(geom.nat()));

    for (int i = 0; i < B.size(); i++){
      geometry<REAL> geom1(geom);
      for (int j = 0; j < geom1.nat(); j++)
	if (geom1.frac)
	  geom1.coord(j) = geom1.cell->cart2frac(B[i]*geom1.pos(j));
	else
	  geom1.coord(j) = B[i]*geom1.coord(j);

      std::vector<vector3<REAL> > T;
      std::vector<permutation> P1;

      //debug
      //std::cout << i << " matrix " << B[i] << "\n";
      
      find_translations(T, P1, geom1, geom, *geom.cell.get(), R);
      //std::cout << "find_transl finished\n";
      
      for (int j=0; j < T.size(); j++){
	//std::cout << i << " " << j << "\n";
	rotrans<REAL> S(T[j],B[i], geom.cell);
	//int n = P1[j].order();
	//std::cout << i << " " << T[j] << P1[j].to_string() << n << "\n";

	int ng_was = G.size();
	auto p = P1[j];
	auto idx = std::find(P.begin(),P.end(),p);
	if (idx==P.end()){
	  P.push_back(p);
	  G.group.push_back(S);	
	}

	//debug

	//G.add(S);
	//G.add(S);
	
	std::cout << "new rotranses:\n" << G.size() << " " << P.size() << "\n";
	for (int ii =ng_was; ii<G.size(); ii++)
	  std::cout << ii << G[ii] << "\n";
	
	
      }
    }

    /*
    std::cout << "Before complete group\n";
    for (int j=0; j<P.size(); j++)
      std::cout << P[j].to_string() << "\n";
    */
    std::cout << "Sizes:\n" << G.size() << " " << P.size() << "\n";
    complete_group(G,P);
    
    std::cout << "After complete group\n";
    std::cout << "Sizes:\n" << G.size() << " " << P.size() << "\n";

    group_analyzer<permutation> AP(P);
    static_table<int> multab(G.size(),G.size());
    for (int i=0; i<G.size(); i++)
      for (int j=0; j<G.size(); j++)
	multab(i,j) = AP.multab(j,i);
    
    //reconstruct_cryst_group(G,multab);

    std::cout << "After reconstruct group\n";

  }

  template<class REAL>
  void find_cryst_symm(array_group<rotrans<REAL> > & G,
                       geometry<REAL > & geom,
                       REAL R = geometry<REAL >
                       ::tol_geom_default){
    std::vector<permutation>  P;
    find_cryst_symm(G,P,geom,R);
  }

  /*! \brief Find the crystalline symmetry group
    @param G(OUT) - the crystalline group of bound rotranslational operations in positionary generator form.
    Bound rotranslational operations are used to make the group finite
    @param geom(IN) - the geometry of unit cell together with lattice vectors, which should be
    stored in geom.cell object. This geometry must be 3D - periodic
    @param R(IN) - the tolerance radius. Symmetry operation is considered valid, if the displacement of atom due to
    this operation is less than R
  */
  template<class REAL>
  void find_cryst_symm(genform_group<rotrans<REAL> > & G,
                       geometry<REAL > & geom,
                       REAL R = geometry<REAL>::tol_geom_default){
    array_group<rotrans<REAL> > G1;
    find_cryst_symm(G1,geom,R);
    generator_form(G,G1);
  }

  // ----------------------------------------------------------------------------
  
  template <class REAL>
  std::vector< rotrans<REAL> > rotrans_circle(const rotrans<REAL> & r, REAL eps)
  {
    rotrans<REAL> F(r.T,r.R),
      G = F;
    std::vector< rotrans<REAL> > circle={G};
    while ( (G.R - matrix3<REAL>::Identity()).norm() > eps )
      {
	G = G*F;
	circle.push_back(G);
      }
    return circle;
  }

  template <class REAL>
  std::vector< index > inspect_rotrans(const rotrans<REAL> & r, REAL eps)
  {
    std::vector<index> I = { index{1,0,0}, index{0,1,0}, index{0,0,1}, index{0,0,0} };
    matrix<REAL>  trns(3,4), cntr(3,4);
    for (int i=0; i<4; i++)
      {
	auto cicrle = rotrans_circle(rotrans_shift(r,I[i]),eps);
	vector3<REAL> c(0);
	for (const auto & g : cicrle)
	  c += g.T;
	c = c/cicrle.size();
	cntr.col(i) = c;
	trns.col(i) = cicrle.back().T;
      }
    for (int i=0; i<3; i++){
      trns.col(i) -= trns.col(3);
      cntr.col(i) -= cntr.col(3);
    }

    //std::cout << "trns\n" << trns << "\ncntr\n" << cntr << "\n";
    
    std::vector<int> indep, dep;
    matrix<REAL> ns;
    nullspace(indep, dep, ns, trns, eps);

    /*
    std::cout << "NS:\n";
    for (int i:indep)
      std::cout << i << " ";
    std::cout << "\n" << ns << "\n";
    */
    
    int d = ns.cols();
    if ( indep.size()>0 && indep.back()==3 ) return {};
    std::vector<index> ret;
    index idx({0,0,0});
    for (int j=0; j<3; j++)
      if ( std::abs(ns.col(d-1)(j) - round(ns.col(d-1)(j)) ) < eps )
	idx(j) = - round(ns.col(d-1)(j));
      else
	return {};
    ret.push_back(idx);
    for (int i=0; i<d-1; i++)
      {
	bool succ = true;
	for (int j=0; j<3; j++)
	  if ( std::abs(ns.col(i)(j) - round(ns.col(i)(j)) ) < eps )
	    idx(j) = - round(ns.col(i)(j));
	  else
	    succ = false;
	if (succ)
	  ret.push_back(idx);
      }
    return ret;
  }

  /*
  template <class REAL>
  std::vector<index> rotrans_grid( const rotrans<REAL> & r, const std::vector<index> & idx)
  {
    int d = idx.size();
    rotrans<REAL> f = rotrans2frac(r);
    std::vector<rotrans<REAL> > f_shft;
    for (const index & I : idx){
      rotrans<REAL> fs(vector3<REAL>(f.T(0)+I(0), f.T(1)+I(1), f.T(2)+I(2)), f.R);
      f_shft.push_back(fs);
    }
    vector3<REAL> pt_cntr = invariant_subspace(f).point;
    std::vector< vector3<REAL> > pt_shft;
    for (const auto & fs : f_shft)
      pt_shft.push_back(invariant_subspace(fs).point);
    periodic_cell<REAL> cl(3);
    for (int i=0; i<d; i++)
      cl(i) = pt_shft[i] - pt_cntr;
    if (d==1){
      int i=0;
      for (int j=1; j<3; j++)
	if (std::abs(cl(0)(i))>std::abs(cl(0)(j)))
	  i=j;
      cl(1) = vector3<REAL>(0,0,0);
      cl(1)(i) = 1e0;
    }
    if (d<=2)
      cl(2) = cl(0).cross(cl(1));
    shape_sphere<REAL> S(std::sqrt(3e0)/2,-1e0*pt_cntr);
    vector3<REAL> fmin = S.fmin(cl);
    vector3<REAL> fmax = S.fmax(cl);
    index grid0{0,0,0},grid1{0,0,0};
    for (int i=0; i<d; i++){
      grid0(i) = std::floor(fmin[i]);
      grid1(i) = std::floor(fmax[i]) + 1;
    }

    /*
    std::cout << "cl\n";
    for (int i=0; i<3; i++)
      std::cout << cl(i) << "\n";
    std::cout << " fmin= " << fmin << " fmax= " << fmax << "\n";
    std::cout << "grid0= " << grid0 << "grid1= " << grid1 << "\n";
    
    
    std::vector<index> grd;
    for ( iterator g(grid0,grid1); !g.end(); g++){
      index I{0,0,0};
      for (int i=0; i<d; i++)
	{
	  I = I + idx[i]*g(i);
	  //std::cout << i << idx[i] << g(i) << I << "\n";
	}
      //std::cout << g << I << "\n";
      auto s = invariant_subspace(rotrans2frac(rotrans_shift(r,I)));
      if (s.point.norm() < std::sqrt(3e0)/2)
	grd.push_back(I);
    }
    return grd;
  }
  */
  // ----------------------------------------------------------------------------
  /*
  template <class REAL>
  std::vector<rotrans<REAL> > relevant_symmetries( const rotrans<REAL> & r, REAL eps){
    auto I = inspect_rotrans(r,eps);
    std::vector<rotrans<REAL> > R;
    if (I.size()==0)
      return R;
    auto rr = rotrans_shift(r,I[0]);
    I.erase(I.begin());
    auto grd = rotrans_grid(rr,I);
    for (auto g : grd){
      auto gr = rotrans_shift(rr,g);
      R.push_back(gr);
    }
    return R;
  }
  
  // ----------------------------------------------------------------------------

  template<class REAL, class TRANSF = int>
  struct highsymmsite{

    subspace3<REAL> S;
    std::vector<TRANSF> group;
    //bool maximal;
    std::vector<int> next;

    highsymmsite(const subspace3<REAL> & s, const std::vector<TRANSF> & g):
      S(s), group(g) {}
    
  };

  template<class T>
  bool belongs_to(const T & t, const std::vector<T> G)
  {
    return std::find(G.begin(),G.end(),t) != G.end();
  }
  
  template <class T>
  bool is_subset(const std::vector<T> G1, const std::vector<T> G){
    for (const T & i : G1)
      if ( ! belongs_set(i,G) )
	return false;
    return true;
  }

  template<class T>
  void merge_to_set( std::vector<T> & A, const std::vector<T> & B){
    for (const T & b:B)
      if (!belongs_to(b,A))
	A.push_back(b);
  }
*/
  /*
  template<class TRANSF>
  void complete_group(std::vector<TRANSF> & S){
    std::vector<TRANSF> N = S;
    while (N.size()>0)
      {
	std::vector<TRANSF> NN;
	for (auto *c : { &S, &N } )
	  for (const auto & i : *c)
	    for (const auto & j : N)
	      {
		TRANSF k = i*j;
		if ( !belongs_to(k,S) && !belongs_to(k,N) && !belongs_to(k,NN))
		  NN.push_back(k);
		k = j*i;
		if ( !belongs_to(k,S) && !belongs_to(k,N) && !belongs_to(k,NN))
		  NN.push_back(k);
	      }
	for (const auto &i : N)
	  S.push_back(i);
	N = NN;
      }
  }
  */
  /*
  template<class REAL>
  int fps_find(const std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
	       const subspace3<REAL> & s){
    int i;
    bool found = false;
    for (i=0; i<H.size(); i++)
      if (H[i].S == s){
	found = true;
	break;
      }
    if (found)
      return i;
    else
      return -1;
  }

  template <class REAL>
  void fps_new( std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
		const subspace3<REAL> & s,
		const std::vector< rotrans<REAL> > & ops){
    int i = fps_find(H,s);
    if (i>-1)
      for (const auto & g: ops)
	{
	  H[i].group.push_back(g);
	}
    else
      H.push_back(highsymmsite<REAL, rotrans<REAL> >(s,ops));
  }

  
  template<class REAL>
  void fps_merge( std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
		  int to, const std::vector<int> & from){
    for (int i:from)
      {
	H[i].next.push_back(to);
	//H[to].group.insert(H[to].group.end(), H[i].group.begin(), H[i].group.end());
      }
  }

  template<class REAL>
  void fps_add( std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
		const subspace3<REAL> & s, int i1, int i2){
    //std::vector< rotrans<REAL> > g(H[i1].group);
    std::vector< rotrans<REAL> > g;
    //merge_to_set(g,H[i2].group);
    //g.insert(g.end(),H[i2].group.begin(), H[i2].group.end());
    H.push_back( highsymmsite<REAL, rotrans<REAL> >(s,g) );
    H[i1].next.push_back(H.size()-1);
    H[i2].next.push_back(H.size()-1);      
  }
  */ /*
  template <class REAL>
  index closest_rotrans(const rotrans<REAL> & g, const vector3<REAL> & r)
  {
    vector3<REAL> t = r - g.R*r - g.T;
    t = g.cell->cart2frac(t);
    return index({round(t(0)), round(t(1)), round(t(2))});
  }

  template <class REAL>
  bool subspace_incell(const subspace3<REAL> & S, const periodic_cell<REAL> & cell)
  {
    REAL eps = geometry<REAL >::tol_geom_default;
    if (S.dim==-1)
      return false;
    if (S.dim==3)
      return true;
    if (S.dim==0){
      vector3<REAL> f = cell.cart2frac(S.point);
      return -eps < f(0) && f(0) < 1e0 + eps &&
				   -eps < f(1) && f(1) < 1e0 + eps &&
							 -eps < f(2) && f(2) < 1e0 + eps;
    }
    else{
      vector3<REAL> f = cell.cart2frac(S.point);
      f -= vector3<REAL>(.5,.5,.5);
      return f.norm() < .5*std::sqrt(3e0)+eps;
    }
  }

  template <class REAL>
  void elementary_subspaces(std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
			    const array_group<rotrans<REAL> > & G, REAL eps)
  {
    group_analyzer<rotrans<REAL>, array_group<rotrans<REAL> > > A(G);
    //for (int g: A.abelian_division()){
    for (int g=1; g < G.size(); g++){
      auto rg = relevant_symmetries(G[g],eps);
      for (const auto & r: rg){
	subspace3<REAL> S = invariant_subspace(r);
	if (S.dim == 3 || S.dim == -1)
	  continue;
	fps_new(H,S,{rg});
      }
    }
    

    /*
    auto cell = G[0].cell;
    for (iterator I({0,0,0},{n,n,n}); !I.end(); I++){
      vector3<REAL> f(I(0), I(1), I(2));
      f = f/n;
      vector3<REAL> r = cell -> frac2cart(f);
      for (int i=0; i<G.size(); i++){
	auto g = G[i];
	index C = closest_rotrans(g,r);
	vector3<REAL> Tprime = g.T + C(0)*(*cell)(0) + C(1)*(*cell)(1) + C(2)*(*cell)(2);
	rotrans<REAL> gprime(Tprime, g.R, cell);
	subspace3<REAL> S = invariant_subspace(gprime);
	if (S.dim == 3 || S.dim == -1)
	  continue;
	//if (subspace_incell(S,*cell))
	  fps_new(H, S, {gprime});
	//else
	//  std::cout << "rejected " << r << " " << f << S.dim << S.point << S.axis << "\n";
      }
      }
  }
*/  /*
  template <class REAL>
  void fps_intersection(std::vector<highsymmsite<REAL, rotrans<REAL> > > & H,
			const periodic_cell<REAL> & cell,
			int d1, int d2)
  {
    for (int i=0; i<H.size(); i++)
      if (H[i].S.dim == d1)
	for (int j=0; j < ( d1==d2 ? i: H.size()) ; j++)
	  if (i!=j && H[j].S.dim == d2){
	    
	    //std::cout << d1<<d2<< " " << i << " " << j << "\n";
	    
	    auto S = H[i].S & H[j].S;
	    //if (S.dim == -1  || !subspace_incell(S,cell) )
	    if (S.dim == -1)
	      continue;	    
	    int k = fps_find(H,S);
	    if (k==-1){
	      std::cout << "before add " << S.dim << S.point << S.axis << "\n";
	      fps_add(H,S,i,j);
	    }
	    else{
	      //std::cout << "before merge " << k << " " << S.dim << S.point << S.axis << "\n";
	      fps_merge(H,k,{i,j});
	    }
	  }    
  }
  
  /*! \brief Finds all point subgroups of crystalline symmetry group.
    Can be used to list all high symmetry sites in the lattice.
    @param groups (OUT)    - std::vector containing point subgroups
    @param subspaces (OUT) - std::vector containing the central points of the point groups
    @param G (IN)          - crystalline symmetry group in array form
  
  //  template<class REAL, bool BOUND>
  template<class REAL>
  void find_point_subgroups(std::vector<array_group<rotrans<REAL> > > & groups,
			    std::vector<subspace3<REAL> > & subspaces,
			    const array_group<rotrans<REAL> > & G,
			    REAL eps)
  {
    auto cell = G[0].cell;
    std::vector<highsymmsite<REAL, rotrans<REAL> > > H;
    elementary_subspaces(H,G,eps);
    fps_intersection(H, *cell, 2, 2);
    fps_intersection(H, *cell, 2, 1);
    fps_intersection(H, *cell, 1, 1);

    /*
    double t1=0e0, t2=0e0, t3=0e0, t4=0e0;
      
    for (int i=0; i<H.size(); i++)
      if (H[i].S.dim == 2)
	for (int j=0; j<i; j++)
	  if (H[j].S.dim == 2){

	    std::clock_t time1 = std::clock();
	    
	    auto S = H[i].S & H[j].S;
	    if (S.dim == -1)
	      continue;

	    std::clock_t time2 = std::clock();
	    t1 += 1e0*(time2-time1)/CLOCKS_PER_SEC*1000;
	    
	    //std::cout << i << " " << j << "\n";
	    int k = fps_find(H,S);

	    std::clock_t time21 = std::clock();
	    t4 += 1000e0*(time21-time2)/CLOCKS_PER_SEC;
	    
	    if (k==-1) {
	      fps_add(H,S,i,j);
	      
	      std::clock_t time3 = std::clock();
	      t2 += 1000e0*(time3-time21)/CLOCKS_PER_SEC;

	    }
	    else{
	      fps_merge(H,k,{i,j});

	      std::clock_t time3 = std::clock();
	      t3 += 1000e0*(time3-time21)/CLOCKS_PER_SEC;

	    }
	  }
    std::cout << " intersection: " << t1 << " fps_add: " << t2
	      << " fps_merge: " << t3 << " fps_find: " << t4 << "\n";
    */
    // Assign subsets
  /*
    for (int i=0; i<H.size(); i++)
      for (int j=0; j<H.size(); j++)
	if (i!=j && H[i].S.within(H[j].S) )
	  fps_merge(H,j,{i});
    for (int i=0; i<H.size(); i++)
      //if (H[i].next.size()==0)
	{
	  subspaces.push_back(H[i].S);
	  /*
	  array_group<matrix3<REAL> > pg;
	  for (auto const & g : H[i].group)
	  pg.add(g.R);
	  array_group<rotrans<REAL> > gg;
	  for (auto & g : H[i].group)
	    gg.add(g);
	  groups.push_back(gg);
	}
  }
  */
  /*
  //double the UC
  periodic_cell<REAL> cell8(*G[0].cell);
  for (int i=0; i<3; i++)
  cell8(i) *= 2;
    
  array_group<rotrans<REAL> > G8("",rotrans<REAL>(matrix3<REAL>::unity, & cell8));

  //array_group<rotrans<REAL> > G8(G);

  int N = G.size();
    
  for (iterator j({0,0,0},{1,1,1}); !j.end(); j++)
  for (int i=0; i < N; i++)
  {
  std::cout << i << j << std::endl;
  G8.add(rotrans<REAL>(G[0].cell->transform(G[i].T,j),G[i].R, & cell8));
  }
      
  std::cout << "size= " << G8.size() << "\n";
  //group_analyzer<rotrans<REAL>, array_group<rotrans<REAL> > > A(G8);
  //group_analyzer<rotrans<REAL>, array_group<rotrans<REAL> > > B(G);
  //auto C = double_group(double_group(double_group(B)));

  /*
  for (int i=0; i<G8.size(); i++)
  for (int j=0; j<G8.size(); j++)
  if (A.multab(i,j)!=C.multab(i,j))
  std::cout << i << " " << j << " " << A.multab(i,j) << " " << C.multab(i,j) << std::endl;
    

  // Form all possible invariant subspaces of abelian subgroups
  // Combine (multiply) subgroups with coinciding subspaces
  std::vector<highsymmsite<REAL, rotrans<REAL> > > H;
  for (int i=0; i < G8.size(); i++)
  {
  auto s = invariant_subspace(G8[i]);
  if (s.dim == -1)
  continue;
  bool found = false;
  for (int j=0; j<H.size(); j++)
  if ( H[j].S == s )
  {
  if (!belongs_to(G8[i],H[j].group))
  H[j].group.push_back(G8[i]);
  found = true;
  break;
  }
  if (!found)
  H.push_back(highsymmsite<REAL, rotrans<REAL> >(s,{G8[i]}));
  }
  //for (int i=0; i<H.size(); i++)
  //  complete_group(H[i].group);

  int inew = 0;

  while (inew < H.size())
  {
  int inewest = H.size();
  for (int ig1 = 0; ig1 < inewest; ig1++)
  for (int ig2 = inew; ig2 < inewest; ig2++)
  if (ig1!=ig2)
  {
  subspace3<REAL> s = H[ig1].S & H[ig2].S;

  bool found = false;
	      
  if (s == H[ig1].S)
  {
  fps_merge(H,ig1,{ig2});
  found = true;
  }
  if (s == H[ig2].S)
  {
  fps_merge(H,ig2,{ig1});
  found = true;
  }

  for (int i=0; i<H.size(); i++)
  if (i!=ig1 && i!=ig2 && s == H[i].S)
  {
  merge_to_set(H[i].group,H[ig1].group);
  merge_to_set(H[i].group,H[ig2].group);
  }
	      
	      
  }
  inew = inewest;
  }

  for (int i=0; i<H.size(); i++)
  //if (H[i].next.size()==0)
  {
  /*
  std::cout << "{";
  for (const auto & j:H[i].group)
  std::cout << j << ","; 
  std::cout << "} "
  std::cout << i << " dim= " << H[i].S.dim << " " << G[0].cell->reduce(H[i].S.point)
  << " " << H[i].S.axis << H[i].next.size();
  for (int j:H[i].next)
  std::cout << " " << j ;
  std::cout << "\n";
  }
  //G = G8;
  }
  */
  /*
    template<class REAL>
    void find_point_subgroups(std::vector<array_group<matrix3<REAL> > > & groups,
    std::vector<vector3<REAL> > &cntrs,
    std::vector<int> & dims,
    const array_group<rotrans<REAL> > & G){}
  */

  /*
    template<class REAL>
    void find_point_subgroups1(std::vector<array_group<matrix3<REAL> > > & groups,
    std::vector<subspace3<REAL> > & subspaces,
    const array_group<rotrans<REAL> > & G){
    
    //std::vector<subspace3<REAL> > subspaces;
    std::vector<std::vector<rotrans<REAL> > > elements;

    //std::cout << "find_point_subs:\n";

    for (const auto & g : G.group){
    auto s = invariant_subspace(g);
    //std::cout << "g= " << g << " s= " << "(" << s.dim << "," << s.pt << "," << s.n << ")\n";
    add_subspace(subspaces,elements, s,{g});
    }

    //debug
    std::setprecision(4);
    std::cout << std::fixed;
    for (int i=0; i<subspaces.size(); i++)
    std::cout << i << " d= " << subspaces[i].dim << " pt= " << subspaces[i].point << " n= " << subspaces[i].axis
    << " ng= " << elements[i].size() << "\n";


    for (int i=0; i<subspaces.size(); i++)
    if ( subspaces[i].dim == 0)
    for (int j=0; j<subspaces.size(); j++)
    if (subspaces[j].within(subspaces[i].point))
    for (const auto & gg : elements[j])
    if ( std::find(elements[i].begin(),elements[i].end(),gg)
    == elements[i].end())
    elements[i].push_back(gg);

    //debug

    std::cout << "\n\n";
    for (int i=0; i<subspaces.size(); i++)
    std::cout << i << "d= " << subspaces[i].dim << " pt= "
    << subspaces[i].point << " n= " << subspaces[i].axis
    << " ng= " << elements[i].size() << "\n";

    /*
    for (int i=0; i<subspaces.size(); i++)
    {
    std::cout << "d= " << subspaces[i].dim << " pt= " << subspaces[i].pt << " n= " << subspaces[i].n;
    for (int j=0; j<elements[i].size(); j++)
    std::cout << elements[i][j];
    std::cout << "\n";
    }
    

    int nnew = 0, n = subspaces.size();
    bool contin = true;
    while (contin){
    int nnewnew = subspaces.size();
    for (int i=nnew; i<nnewnew; i++)
    if (subspaces[i].dim > 0)
    for (int j=0; j<n; j++)
    if (subspaces[j].dim > 0){
    add_subspace(subspaces,elements,
    subspaces[i] & subspaces[j], elements[i]);
    add_subspace(subspaces,elements,
    subspaces[i] & subspaces[j], elements[j]);
    }
    contin = nnewnew < subspaces.size();
    nnew = nnewnew;
    }

    /*
    std::vector<int> idx;
    for (int i=0; i<subspaces.size(); i++) idx.push_back(i);

    std::sort(idx.begin(), idx.end(),
    [&subspaces](int i, int j) -> bool
    { return subspaces[i].dim < subspaces[j].dim; }
    );
    reorder(subspaces,idx);
    reorder(elements, idx);
    
    for (int i=0; i<subspaces.size(); i++){
    groups.push_back(array_group<matrix3<REAL> >());
    for (int j=0; j<elements[i].size(); j++)
    groups[i].generate(elements[i][j].R);
    }

    /*
    for (int d=0; d<3; d++)
    for (int i=0; i<subspaces.size(); i++)
    if (subspaces[i].dim==d){
    cntrs.push_back(subspaces[i].point);
    dims.push_back(d);
    subs.push_back(array_group<matrix3<REAL> >());
    int n=subs.size()-1;
    for (int j=0; j<elements[i].size(); j++)
    groups[n].generate(elements[i][j].R);
    }
    
    }
  */
  /*
    template<class REAL>
    void find_point_subgroups(std::vector<array_group<matrix3<REAL> > > & groups,
    std::vector<subspace3<REAL> > & subspaces,
    const array_group<rotrans<REAL> > & G)
    {
    array_group<rotrans<REAL> > G1;
    G1.group.clear();
    for (const auto & x : G.group)
    G1.group.push_back(rotrans<REAL>(x.T,x.R));

    find_point_subgroups(groups,subspaces,G1);
    }
  */

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

  template<class REAL>
  void py_fix4_cryst_group(array_group<rotrans<REAL> > & G, const std::vector<std::vector<int> > & M){
    int N = G.size();
    static_table<int> MM(N,N);
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
	MM(i,j) = M[i][j];
    fix4_cryst_group(G,MM);
  }

  template<class REAL>
  void py_reconstruct_cryst_group(array_group<rotrans<REAL> > & G,
				  std::vector<permutation> & P){
    group_analyzer<permutation> AP(P);
    static_table<int> multab(G.size(),G.size());
    for (int i=0; i<G.size(); i++)
      for (int j=0; j<G.size(); j++)
	multab(i,j) = AP.multab(j,i);
    
    reconstruct_cryst_group(G,multab);

  }

  
  template<class REAL>
  void py_find_cryst_symm1(array_group<rotrans<REAL> > & G,
			   std::vector<permutation> &P,
                           geometry<REAL > & geom,
                           REAL R = geometry<REAL>
                           ::tol_geom_default)
  {  find_cryst_symm(G,P,geom,R); }

  template<class REAL>
  void py_find_cryst_symm2(genform_group<rotrans<REAL> > & G,
                           geometry<REAL > & geom,
                           REAL R = geometry<REAL>
                           ::tol_geom_default)
  { find_cryst_symm(G,geom,R); }

  template<class REAL>
  py::object py_find_cryst_symm3( geometry<REAL > & geom,
                           REAL R = geometry<REAL>
                           ::tol_geom_default)
  {
    array_group<rotrans<REAL> >  * G = new array_group<rotrans<REAL> >;
    std::vector<permutation> * P = new std::vector<permutation>;
    find_cryst_symm(*G,*P,geom,R);
    return py::make_tuple(*G,*P);
  }
  /*
  template<class REAL>
  void py_find_point_subgroups(py::list & groups, py::list &subspaces,
			       const array_group<rotrans<REAL> > & G,
			       REAL eps){
    //std::vector<array_group<matrix3<REAL> > >  vgroups;
    std::vector<array_group<rotrans<REAL> > >  vgroups;
    std::vector<subspace3<REAL> > vsubspaces;
    find_point_subgroups(vgroups,vsubspaces,G,eps);
    for (int i = 0; i < vgroups.size(); i++){
      groups.append(vgroups[i]);
      subspaces.append(vsubspaces[i]);
    }
  }
  */
  /*
    template<class REAL, bool BOUND>
    void py_find_point_subgroups2(py::list & subs,
    py::list &cntrs,
    py::list & dims,
    const array_group<rotrans<REAL,BOUND> >
    & G){
    std::vector<array_group<matrix3<REAL> > >  vsubs;
    std::vector<vector3<REAL> > vcntrs;
    std::vector<int> vdims;
    find_point_subgroups(vsubs,vcntrs,vdims,G);
    for (int i=0; i<vsubs.size(); i++){
    subs.append(vsubs[i]);
    cntrs.append(vcntrs[i]);
    dims.append(vdims[i]);
    }
    }
  */

  /*
    template<class REAL>
    void py_find_point_subgroups2(bp::list & subs, bp::list &cntrs,
    const array_group<rotrans<REAL> > & G)
    {
    std::vector<array_group<matrix3d<REAL> > >  vsubs;
    std::vector<vector3<REAL> > vcntrs;
    find_point_subgroups(vsubs,vcntrs,G);
    for (int i=0; i<vsubs.size(); i++)
    {
    subs.append(vsubs[i]);
    cntrs.append(vcntrs[i]);
    }
    }
  */

  template<class REAL>
  py::object py_find_translations( geometry<REAL > & g1,
				   geometry<REAL > & g2,
				   const periodic_cell<REAL> &cell,
				   REAL R,
				   bool return_permut = false){
    std::vector<vector3<REAL> > t;
    std::vector<permutation> perm;
    find_translations(t, perm, g1, g2, cell, R);
    py::list transl;
    py::object res;
    for (const auto & tt : t)
      transl.append(tt);
    if (return_permut)
      res = py::make_tuple(transl,perm);
    else
      res = transl;
    return res;
  }


#endif

}

#endif
