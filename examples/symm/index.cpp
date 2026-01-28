//#define PY_EXPORT
#include <iostream>
#include <data/errors.hpp>
#include <symm/index.hpp>
//#include <symm/group_theory.hpp>
#include <mathf/lace3d.hpp>

//using namespace qpp;

// This example shows how to manually build D4 and D4h point groups

int main()
{
  qpp::index a({0,0,0});
  qpp::index b({2,2,2});
  qpp::index_range R({0,0,0,0},{5,3,2,1});
  for (
       qpp::iterator i({0,0,0,0},{5,3,2,1}); !i.end(); i++)
    std::cout << i<< "\n";
}
