#include <iostream>
#include <mathf/lace3d.hpp>

using namespace qpp;

// This example shows how to manually build D4 and D4h point groups

std::vector<double> understand(const std::vector<double> & vals){
  return vals;
}

int main()
{
  double a=42, b=3.14;
  vector3<double> n={a,0,0}, n1={b,a,0};

  vector3<double> n2 = n%n1;
  n /= n.norm();
  n1 /= n1.norm();
  n2 /= n2.norm();

  vector2<double> r(0,1);
  std::vector<double> v=understand({0., 1., 1., 0.});
  matrix2<double> M = matrix2<double>::create({0.,1.,1.,0.});
  //std::cout << fmt::format("{}{}{}\n", n,n1,n2);
  std::cout << n << n1 << n2 <<  "\n";
  std::cout << n.x() << n.y()<< n.z() << "\n";
  std::cout << r.x() << " " <<  r.y()<< "\n"<< M << "\n";
}
