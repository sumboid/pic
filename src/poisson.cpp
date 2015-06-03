#include <ts-ng/system/System.h>
#include "poisson-fragment.h"
#include "mesh.h"
#include <tuple>
#include <vector>
#include <string>
#include <set>
#include <fstream>

using ts::system::System;
using ts::type::ID;

namespace {
    const int nx = 32;
    const int ny = 32;
    const int nz = 32;

    const double sx = 4;
    const double sy = 4;
    const double sz = 4;

    const double hx = sx / (nx - 2);
    const double hy = sy / (ny - 2);
    const double hz = sz / (nz - 2);

    const double cx = sx / 2;
    const double cy = sy / 2;
    const double cz = sz / 2;

    const double PM = 1.0;

    const int SPLITY = 8;

}


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

std::tuple<int, int> interval(int current, int all, int size) {
   int end = size % all;
   int partSize = size / all;

   if(current >= all - end) {
       size_t some = current - (all - end);
       return std::tuple<int, int>(partSize * current + some, partSize * current + some + partSize + 1);
   }
   else {
       return std::tuple<int, int>(partSize * current, partSize * current + partSize);
   }
}

void fillParticles(Mesh* mesh, int bx, int ex, int by, int ey) {
    std::ifstream file("boom.dat");
    double x, y, z;

    for(int i = 0; i < 10000; ++i) {
        file >> x >> y >> z;
        double nx = x / hx;
        int ix = nx;

        double ny = y / hy;
        int iy = ny;

        if(ix >= bx && ix < ex && iy >= by && iy < ey) {
            mesh->addParticle(new Particle(x - (bx * hx) + hx, y - (by * hy) + hy, z));
        }
    }

    file.close();
}

Fragment* createFragment(int cx, int ax, int cy, int ay) {
    int bx, ex, by, ey;
    std::tie(bx, ex) = interval(cx, ax, nx);
    std::tie(by, ey) = interval(cy, ay, ny);

    int sfx = ex - bx;
    int sfy = ey - by;

    Mesh* mesh = new Mesh(sfx, sfy, nz);
    mesh->setID(cx, cy, 0);
    mesh->setSizes(sx, sy, sz, nx, ny, nz);

    fillParticles(mesh, bx, ex, by, ey);


    Fragment* f = new Fragment(ts::type::ID(cx, cy, 0));
    if(bx > 0 && ex < nx && by > 0 && ey < ny) {
        mesh->setSides(true, true, false, false, false, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }
    else if(bx == 0 && ex < nx && by > 0 && ey < ny) {
        mesh->setSides(true, true, false, false, true, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }
    else if(bx > 0 && ex == nx && by > 0 && ey < ny) {
        mesh->setSides(true, true, false, false, false, true);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }
    else if(bx > 0 && ex < nx && by == 0 && ey < ny) {
        mesh->setSides(true, true, true, false, false, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
    }
    else if(bx == 0 && ex < nx && by == 0 && ey < ny) {
        mesh->setSides(true, true, true, false, true, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
    }
    else if(bx > 0 && ex == nx && by == 0 && ey < ny) {
        mesh->setSides(true, true, true, false, false, true);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy + 1, 0), cx);
    }
    else if(bx > 0 && ex < nx && by > 0 && ey == ny) {
        mesh->setSides(true, true, false, true, false, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }
    else if(bx == 0 && ex < nx && by > 0 && ey == ny) {
        mesh->setSides(true, true, false, true, true, false);
        f->addNeighbour(ts::type::ID(cx + 1, cy, 0), cx + 1);

        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }
    else if(bx > 0 && ex == nx && by > 0 && ey == ny) {
        mesh->setSides(true, true, false, true, false, true);
        f->addNeighbour(ts::type::ID(cx - 1, cy, 0), cx - 1);

        f->addNeighbour(ts::type::ID(cx, cy - 1, 0), cx);
    }


    f->setMesh(mesh);
    std::cout << cx << ": Particles number " << f->weight() << std::endl;
    std::cout << cx << ": [" << bx << ", " << ex << ")" << std::endl;
    return f;
}

int main()
{
  System* s = createSystem();
  auto id = s->id();
  auto size = s->size();
  s->setBalancer(lazybalancer);

  for(uint64_t i = 0; i < SPLICY; ++i) {
        Fragment* f = createFragment(id, size, );
        s->addFragment(f);
  }

  s->run();

  return 0;
}
