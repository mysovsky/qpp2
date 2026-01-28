#include <iostream>
#include <functional>
//#include <data/altref.hpp>
#include <geom/geom.hpp>
#include <mathf/lace3d.hpp>



int main(){
  qpp::geometry<double> g;
  g.add("Co",0,0,0);
  g.add("O",1,0,0);
  g.add("O",0,1,0);
  g.add("O",0,0,1);
  g.build_typetable();
  for (int t=0; t< g.typetable()->n_atom_types(); t++)
    std::cout << t << " " << g.typetable()->atomic_type(t)<< "\n";
}
