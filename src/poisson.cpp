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

namespace {
    const int nx = 20;
    const int ny = 20;
    const int nz = 20;

    const double sx = 4;
    const double sy = 4;
    const double sz = 4;

    const double hx = sx / (nx);
    const double hy = sy / (ny);
    const double hz = sz / (nz);

    const double cx = sx / 2;
    const double cy = sy / 2;
    const double cz = sz / 2;

    const double PM = 1.0;

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

void fillFiBoom(Mesh* mesh, int b, int e) {
    double x2,y2,z2,r;
    int i,k,l;

    if(b == 0) {
        x2 = (hx / 2 - cx) * (hx / 2 - cx);
        for (k=0;k<ny;k++)
        {
            y2 = (k * hy + hy / 2 - cy) * (k * hy + hy / 2 - cy);
            for (l=0;l<nz;l++)
            {
                z2 = (l * hz + hz / 2 - cz) * (l * hz + hz / 2 - cz);
                r=sqrt(x2 + y2 + z2);
                mesh->setFi(0, k, l, -PM / r);
            }
        }
    }

    if(e == nx) {
        x2 = (sx - hx / 2 - cx) * (sx - hx / 2 - cx);
        for (k=0;k<ny;k++)
        {
            y2 = (k * hy + hy / 2 - cy) * (k * hy + hy / 2 - cy);
            for (l=0;l<nz;l++)
            {
                z2 = (l * hz + hz / 2 - cz) * (l * hz + hz / 2 - cz);
                r=sqrt(x2 + y2 + z2);
                mesh->setFi(nx - 1, k, l, -PM / r);
            }
        }
    }

    y2 = (hy / 2 - cy) * (hy / 2 - cy);
    for (i=b;i<e;i++)
    {
        x2 = (i * hx + hx / 2 - cx) * (i * hx + hx / 2 - cx);
        for (l=0;l<nz;l++)
        {
            z2 = (l * hz + hz / 2 - cz) * (l * hz + hz / 2 - cz);
            r=sqrt(x2 + y2 + z2);
            mesh->setFi(i, 0, l, -PM / r);
        }
    }

    y2 = (sy - hy / 2 - cy) * (sy - hy / 2 - cy);
    for (i=b;i<e;i++)
    {
        x2 = (i * hx + hx / 2 - cx) * (i * hx + hx / 2 - cx);
        for (l=0;l<nz;l++)
        {
            z2 = (l * hz + hz / 2 - cz) * (l * hz + hz / 2 - cz);
            r=sqrt(x2 + y2 + z2);
            mesh->setFi(i, ny - 1, l, -PM / r);
        }
    }

    z2 = (hz / 2 - cz) * (hz / 2 - cz);
    for (i=b;i<e;i++)
    {
        x2 = (i * hx + hx / 2 - cx) * (i * hx + hx / 2 - cx);
        for (k=0;k<ny;k++)
        {
            y2 = (k * hy + hy / 2 - cy) * (k * hy + hy / 2 - cy);
            r=sqrt(x2 + y2 + z2);
            mesh->setFi(i, k, 0, -PM / r);
        }
    }

    z2 = (sz - hz / 2 - cz) * (sz - hz / 2 - cz);
    for (i=b;i<e;i++)
    {
        x2 = (i * hx + hx / 2 - cx) * (i * hx + hx / 2 - cx);
        for (k=0;k<ny;k++)
        {
            y2 = (k * hy + hy / 2 - cy) * (k * hy + hy / 2 - cy);
            r=sqrt(x2 + y2 + z2);
            mesh->setFi(i, k, nz - 1, -PM / r);
        }
    }
}

void fillFi(Mesh* mesh, int b, int e) {
    double x,y,z,x2,y2,z2,x2py2,x2pz2,r;
    int i,k,l;

    if(b == 0) {
        x=-0.5*hx-cx;
        x2=x*x;
        for (k=0;k<ny;k++)
        { y=(k-0.5)*hy-cy;
            x2py2=x2+y*y;
            for (l=0;l<nz;l++)
            { z=(l-0.5)*hz-cz;
                r=sqrt(x2py2+z*z);
                mesh->setFi(0, k, l, -PM / r);
            }
        }
    }

    if(e == nx) {
        x=sx+0.5*hx-cx;
        x2=x*x;
        for (k=0;k<ny;k++)
        { y=(k-0.5)*hy-cy;
            x2py2=x2+y*y;
            for (l=0;l<nz;l++)
            { z=(l-0.5)*hz-cz;
                r=sqrt(x2py2+z*z);
                mesh->setFi(nx - 1, k, l, -PM / r);
            }
        }
    }

    y=-0.5*hy-cy;
    y2=y*y;
    for (i=b;i<e;i++)
    { x=(i-0.5)*hx-cx;
        x2py2=x*x+y2;
        for (l=0;l<nz;l++)
        { z=(l-0.5)*hz-cz;
            r=sqrt(x2py2+z*z);
            mesh->setFi(i, 0, l, -PM / r);
        }
    }

    y=sy+0.5*hy-cy;
    y2=y*y;
    for (i=b;i<e;i++)
    { x=(i-0.5)*hx-cx;
        x2py2=x*x+y2;
        for (l=0;l<nz;l++)
        { z=(l-0.5)*hz-cz;
            r=sqrt(x2py2+z*z);
            mesh->setFi(i, ny - 1, l, -PM / r);
        }
    }

    z=-0.5*hz-cz;
    z2=z*z;
    for (i=b;i<e;i++)
    { x=(i-0.5)*hx-cx;
        x2pz2=x*x+z2;
        for (k=0;k<ny;k++)
        { y=(k-0.5)*hy-cy;
            r=sqrt(x2pz2+y*y);
            mesh->setFi(i, k, 0, -PM / r);
        }
    }

    z=sz+0.5*hz-cz;
    z2=z*z;
    for (i=b;i<e;i++)
    { x=(i-0.5)*hx-cx;
        x2pz2=x*x+z2;
        for (k=0;k<ny;k++)
        { y=(k-0.5)*hy-cy;
            r=sqrt(x2pz2+y*y);
            mesh->setFi(i, k, nz - 1, -PM / r);
        }
    }
}

void fillParticles(Mesh* mesh, int b, int e) {
    std::ifstream file("circle.dat");
    double x, y, z;

    for(int i = 0; i < 10000; ++i) {
        file >> x >> y >> z;
        int xi = x / hx;
        if(xi >= b && xi < e) {
            mesh->addParticle(new Particle(x - (b * hx) + (b == 0 ? 0 : hx), y, z));
        }
    }

    file.close();
}

Fragment* createFragment(int cur, int all, int size) {
    int b, e;
    std::tie(b, e) = interval(cur, all, size);
    int s = e - b;

    Mesh* mesh = new Mesh(s, ny, nz);
    fillFi(mesh, b, e);
    fillParticles(mesh, b, e);
    mesh->setID(cur, 0, 0);
    Fragment* f = new Fragment(ts::type::ID(cur, 0, 0));
    if(b > 0 && e < size) {
        mesh->setSides(true, true, true, true, false, false);
        f->addNeighbour(ts::type::ID(cur + 1, 0, 0), cur + 1);
        f->addNeighbour(ts::type::ID(cur - 1, 0, 0), cur - 1);
    }
    else if(b > 0) {
        mesh->setSides(true, true, true, true, false, true);
        f->addNeighbour(ts::type::ID(cur - 1, 0, 0), cur - 1);
    }
    else if(e < size) {
        mesh->setSides(true, true, true, true, true, false);
        f->addNeighbour(ts::type::ID(cur + 1, 0, 0), cur + 1);
    }

    f->setMesh(mesh);
    return f;
}

int main()
{
  System* s = createSystem();
  auto id = s->id();
  auto size = s->size();

  Fragment* f = createFragment(id, size, nx);

  s->setBalancer(lazybalancer);
  s->addFragment(f);
  s->run();

  return 0;
}
