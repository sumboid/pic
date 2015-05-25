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
  return 0;
}
