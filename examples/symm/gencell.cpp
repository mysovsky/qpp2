#include <iostream>
#include <symm/groups.hpp>
#include <symm/gen_cell.hpp>
#include <mathf/lace3d.hpp>
#include <symm/shoenflis.hpp>
#include <memory>

int main()
{
  qpp::array_group<qpp::matrix3<float> > Garr = qpp::shnfl<float>::Ih();
  qpp::genform_group<qpp::matrix3<float> > G;
  generator_form(G,Garr);
 // Print the generators
  std::cout << "Number of Ih group generators = " << G.DIM << std::endl << "Generators are:" << std::endl;
  for (int i=0; i<G.DIM; i++)
    std::cout << G.generators[i] << std::endl << "the order is " << G.end()(i)+1 << std::endl;  

  std::shared_ptr<qpp::periodic_cell<float>> a;
  qpp::gen_cell<float,  qpp::matrix3<float>> b(G);

  a=std::shared_ptr<qpp::periodic_cell<float>>(new qpp::gen_cell<float,  qpp::matrix3<float>>(G), [](qpp::periodic_cell<float> *p){ delete p;});

  qpp::vector3<float> r0(1,.1,.1);

  std::cout << a->symmetrize(r0,.5)<< "\n";
}
