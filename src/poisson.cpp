#include <ts/system/System.h>
#include "poisson-fragment.h"
#include "mesh.h"
#include <tuple>
#include <vector>
#include <string>
#include <set>
#include <fstream>

using ts::system::System;
using ts::type::ID;


System* createSystem() {
  FragmentTools* ct = new FragmentTools;
  ReduceDataTools* rt = new ReduceDataTools;
  return new System(ct, rt);
}

std::map<int, double> balancer(uint64_t o, std::map<int, uint64_t> ns) {
  uint64_t co = o;
  std::map<int, double> r;

  for(auto& n : ns) {
    if(n.second < co) {
      uint64_t m = (co - n.second) / 2;
      co -= m;
      r[n.first] = m / o;
    }
  }

  return r;
}

std::map<int, double> lazybalancer(uint64_t, std::map<int, uint64_t>) {
  return std::map<int, double>();
}

int main()
{
  System* s = createSystem();
  Fragment* fragment = new Fragment(ID(s->id(), 0, 0));
  fragment->setbnd(1, 1, 1);

  Mesh* mesh = new Mesh(10, 10, 10);
  if(s->id() == 0) {
      mesh->setCorners(true, true, true, true, true, false);
      fragment->addNeighbour(ID(1, 0, 0), 1);
  }
  else if (s->id() == 1) {
      mesh->setCorners(true, true, true, true, false, true);
      fragment->addNeighbour(ID(0, 0, 0), 0);
  }
  mesh->sethxyz(1, 1, 1);
  mesh->setw(1);
  mesh->setcoef(1);

  fragment->setMesh(mesh);

  s->addFragment(fragment);

  s->setBalancer(lazybalancer);
  s->run();

  delete s;
  return 0;
}
