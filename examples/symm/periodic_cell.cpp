#include<iostream>
#include<symm/cell.hpp>

int main(){
  qpp::vector3<double> a(0,2.73,2.73);
  qpp::vector3<double> b(2.73,0,2.73);
  qpp::vector3<double> c(2.73,2.73,0);

  qpp::periodic_cell<double> cl(a,b,c);
  qpp::index_range R({-2,-2,-2},{2,2,2,});
  auto i = R.__iter__();
  for (; !i.end(); i++)
    std::cout << cl.transform(qpp::vector3<double>(0),i)<< "\n";
}

