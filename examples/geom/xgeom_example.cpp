//#define PY_EXPORT
#undef PY_EXPORT
#undef QPPCAD_PY_EXPORT
#include <geom/xgeom.hpp>
#include <string>
#include <vector>
#include <iostream>

int main(){
  std::vector<std::string> fn({"mass","optx","opty","optz", "cmnt"});
  std::vector<std::string> ft({"real","bool","bool","bool", "str"});
  qpp::periodic_cell<double> cl(3);
  //std::cout << fn << "\n";
 qpp::xgeometry<double> g(fn,ft, cl);

 std::vector<std::string> ffn;
 std::vector<int> fft;
 g.get_format(ffn,fft);
 std::cout<< g.findex("mass") << " " << g._xfields[0].index() << " "<< (qpp::type_array&g.xftype("real"))<< " " << (11&8) <<"\n";
 g.add("O",.0, .0, .0, 16.0, false, false, false, "Ocntr");
 g.add("H",1, .0, .0, 1.0, true, true, true, "Hcut");
 g.add("H",.0, 1., .0, 1.0, true, true, true, "Hcut");
 g.insert(1,"X",.1,.1,.1, .52, false, true, false, "Whoknowswhat");
 for (int f=0; f<g.nfields(); f++)
   std::cout << g._xfields[f].index() << " fld  ";
 std::cout << "\n";
 //g.fillxfields(1,.042,true,false,true);
 //g.xfield<double>(0,1)=.042; 
 // g.xfield<bool>(1,1)=false;  
 //g.xfield<bool>(2,1)=true;
 //g.xfield<bool>(3,1)=false;
 /*
 for (int i=0; i<g.nat(); i++)
   std::cout << g.atom(i) << g.pos(i) << g.xfield<double>(0,i) << " "
	     <<  g.xfield<bool>(1,i)<< " "
	     <<  g.xfield<bool>(2,i) << " "
	     <<  g.xfield<bool>(3,i) << " "
	     <<  g.xfield<std::string>(4,i) << "\n";
 */
 fn={"mass","optf", "cmnt"};
 ft = {"real","b", "l(i)"};

 qpp::xgeometry<double>::xtrafieldtype a = g.newcolumn(g.xftype("l(i)"));
 qpp::xgeometry<double> g1(fn,ft, cl);
 std::vector<int> v = {5,4,3,2,1};
 std::cout << "before adding\n";
 g1.add("Se",.1,.1,.1,3.14,true,v);
 //std::cout << "added 1\n";
 v={1,2,3};
 g1.add("S",.11,.11,.31,-.8,true,v);
 v={2,1};
 g1.add("H",.01,.01,-31,.4,true,v);
 //std::cout << "added 2\n";
 //std::cout << "added 3\n";
 std::cout << "before printout\n";
 for (int i =0; i< g1.nat(); i++){   
   std::cout << g1.atom(i) << " " <<g1.pos(i) << g1.xfield<double>(0,i) <<" " <<  g1.xfield<bool>(1,i);
   std::cout << "[";
   for (int k: g1.xfield<std::vector<int>>(2,i))
     std::cout << k << " ";
   std::cout << "]\n";       
 }
 g1.insert(1,"O",.11,.11,.31,-1.2,false,v);
 for (int i =0; i< g1.nat(); i++){   
   std::cout << g1.atom(i) << " " <<g1.pos(i) << g1.xfield<double>(0,i) <<" " <<  g1.xfield<bool>(1,i);
   std::cout << "[";
   for (int k: g1.xfield<std::vector<int>>(2,i))
     std::cout << k << " ";
   std::cout << "]\n";       
 }

 std::vector<qpp::xgeometry<double>::fieldtypes> f;
 g1.get_fields(0,f);
 /*
 for (int i =0; i< g1.nat(); i++){
   auto l = g1.py_getitem(i);
   py::print(l);
 }
 */
}
