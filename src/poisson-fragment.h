#pragma once

#include <ts/types/Fragment.h>
#include <ts/types/FragmentTools.h>
#include <ts/types/ID.h>
#include <ts/types/ReduceData.h>
#include <ts/types/ReduceDataTools.h>
#include <ts/util/Uberlogger.h>
#include <ts/util/Arc.h>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <math.h>

#include "mesh.h"

using ts::type::ID;

// Reduction
class ReduceData: public ts::type::ReduceData {
private:
    double max;
public:
  ReduceData(): max(0) {}
  ReduceData(double _max): max(_max) {}
  double get() { 
      return max; 
  }
  virtual ~ReduceData() {}
  virtual ts::type::ReduceData* copy() { return new ReduceData(max); }
};

class ReduceDataTools: public ts::type::ReduceDataTools {
public:
  ~ReduceDataTools() {}
  ts::Arc* serialize(ts::type::ReduceData* rd) {
    ts::Arc* arc = new ts::Arc;
    ts::Arc& a = *arc;

    ReduceData* nrd = reinterpret_cast<ReduceData*>(rd);
    a << nrd->get();
    return arc;
  }

  ts::type::ReduceData* deserialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    double max;
    a >> max;
    ReduceData* rd = new ReduceData(max);
    return rd;
  }

  ts::type::ReduceData* reduce(ts::type::ReduceData* rd1,
                               ts::type::ReduceData* rd2) {
    ReduceData* nrd1 = reinterpret_cast<ReduceData*>(rd1);
    ReduceData* nrd2 = reinterpret_cast<ReduceData*>(rd2);

    return nrd1->get() > nrd2->get() ? rd1->copy() : rd2->copy();
  }
};


class Fragment: public ts::type::Fragment {
friend class FragmentTools;
private:
  Mesh* mesh;
  unsigned int bx, by, bz;

  double max;

  //Boundary stuff
  MeshBoundary* xy1;
  MeshBoundary* xy2;
  MeshBoundary* xz1;
  MeshBoundary* xz2;
  MeshBoundary* yz1;
  MeshBoundary* yz2;
  bool Fi;
  bool Ro;

public:
  Fragment(ts::type::ID id): ts::type::Fragment(id), max(0) {
    Fi = false;
    Ro = false;
  }

  ~Fragment() {
    if(isBoundary()) {
        delete xy1;
        delete xy2;
        delete xz1;
        delete xz2;
        delete yz1;
        delete yz2;
    } else {
        delete mesh;
    }
  }

  void setMesh(Mesh* _mesh) {
    mesh = _mesh;
  }

  void setbnd(int _bx, int _by, int _bz) {
    bx = _bx;
    by = _by;
    bz = _bz;
  }

  void runStep(std::vector<ts::type::Fragment*> fs) override {
      std::cout << "State: (" << iteration()  << ", " << progress()  << ")" << std::endl;
      if(iteration() == 0 && progress() == 0) {
        max = mesh->process();
      }

      else {
          auto x = id().c[0];
          auto y = id().c[1];
          auto z = id().c[2];

          struct _ {
            bool side;
            Fragment* fragment;
            _ () : side(false) {}
            void set(Fragment* f) {
                fragment = f;
                side = true;
            }
          } nbh[6]; // xy1, xy2, yz1, yz2, xz1, xz2

          for(auto nb : fs) {
              auto nid = nb->id();
              auto nx = nid.c[0];
              auto ny = nid.c[1];
              auto nz = nid.c[2];

              Fragment* rnb = reinterpret_cast<Fragment*>(nb);
              std::cout << "(" << x << ", " << y << ", " << z << ") -> (" << nx << ", " << ny << ", " << nz << ")" << std::endl;
              if      (nx < x) nbh[2].set(rnb);
              else if (nx > x) nbh[3].set(rnb);
              else if (ny < y) nbh[4].set(rnb);
              else if (ny > y) nbh[5].set(rnb);
              else if (nz < z) nbh[0].set(rnb);
              else if (nz > z) nbh[1].set(rnb);
          }

          for(int i = 0; i < 6; ++i)
              if(nbh[i].side != false) {
                  if      (i == 0) mesh->setBoundaryXY1(nbh[i].fragment->xy1);
                  else if (i == 1) mesh->setBoundaryXY2(nbh[i].fragment->xy2);
                  else if (i == 2) mesh->setBoundaryYZ1(nbh[i].fragment->yz1);
                  else if (i == 3) mesh->setBoundaryYZ2(nbh[i].fragment->yz2);
                  else if (i == 4) mesh->setBoundaryXZ1(nbh[i].fragment->xz1);
                  else if (i == 5) mesh->setBoundaryXZ2(nbh[i].fragment->xz2);
              }

          max = mesh->process();
      }

      nextIteration();
      saveState();
      setUpdate();
      setReduce();
      setNeighbours(iteration(), progress());
  }

  ReduceData* reduce() override {
    return new ReduceData(max);
  }

  ReduceData* reduce(ts::type::ReduceData* rd) override {
    ReduceData* rrd = reinterpret_cast<ReduceData*>(rd);
    double amax = rrd->get();
    return amax > max ? new ReduceData(amax) : new ReduceData(max);
  }

  void reduceStep(ts::type::ReduceData* rd) override {
    max = reinterpret_cast<ReduceData*>(rd)->get();
  }

  Fragment* getBoundary() override {
    MeshBoundary** boundary = mesh->getBoundary();
    Fragment* fragment = new Fragment(id());
    fragment->xy1 = boundary[0];
    fragment->xy2 = boundary[1];
    fragment->xz1 = boundary[2];
    fragment->xz2 = boundary[3];
    fragment->yz1 = boundary[4];
    fragment->yz2 = boundary[5];

    fragment->setBoundary();
    fragment->Fi = true;
    fragment->Ro = false;
    // TODO: some boundaries may be empty
    return fragment;
  }

  void serializeBoundary(ts::Arc* arc) {
    ts::Arc& a = *arc;
    a << id().c[0];
    a << id().c[1];
    a << id().c[2];

    a << Fi;
    a << Ro;

    if(Fi == true) {
      xy1->serializeFi(arc);
      xy2->serializeFi(arc);
      yz1->serializeFi(arc);
      yz2->serializeFi(arc);
      xz1->serializeFi(arc);
      xz2->serializeFi(arc);
    } else if(Ro == true) {
      xy1->serializeRo(arc);
      xy2->serializeRo(arc);
      yz1->serializeRo(arc);
      yz2->serializeRo(arc);
      xz1->serializeRo(arc);
      xz2->serializeRo(arc);
    }
  }

  static Fragment* deserializeBoundary(ts::Arc* arc) {
    ts::Arc& a = *arc;
    uint64_t x, y, z;
    a >> x;
    a >> y;
    a >> z;

    Fragment* f = new Fragment(ts::type::ID(x, y, z));
    f->setBoundary();
    a >> f->Fi;
    a >> f->Ro;

    if(f->Fi == true) {
      f->xy1 = MeshBoundary::deserializeFi(arc);
      f->xy2 = MeshBoundary::deserializeFi(arc);
      f->yz1 = MeshBoundary::deserializeFi(arc);
      f->yz2 = MeshBoundary::deserializeFi(arc);
      f->xz1 = MeshBoundary::deserializeFi(arc);
      f->xz2 = MeshBoundary::deserializeFi(arc);
    } else if(f->Ro == true) {
      f->xy1 = MeshBoundary::deserializeRo(arc);
      f->xy2 = MeshBoundary::deserializeRo(arc);
      f->yz1 = MeshBoundary::deserializeRo(arc);
      f->yz2 = MeshBoundary::deserializeRo(arc);
      f->xz1 = MeshBoundary::deserializeRo(arc);
      f->xz2 = MeshBoundary::deserializeRo(arc);
    }

    return f;
  }

  void serialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    a << id().c[0] << id().c[1] << id().c[2];
    a << bx << by << bz;
    mesh->serialize(arc);
  }

  static Fragment* deserialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    int x, y, z, bx, by, bz;

    a >> x;
    a >> y;
    a >> z;
    a >> bx;
    a >> by;
    a >> bz;

    Fragment* f = new Fragment(ts::type::ID(x, y, z));
    f->setbnd(bx, by, bz);
    f->setMesh(new Mesh(arc));

    return f;
  }

  Fragment* copy() override {
    if(!isBoundary()) {
        Fragment* f = new Fragment(id());
        f->setMesh(mesh->copy());
        f->bx = bx;
        f->by = by;
        f->bz = bz;
        return f;
    } else {
        Fragment* f = new Fragment(id());
        f->xy1 = xy1->copy();
        f->xy2 = xy2->copy();
        f->xz1 = xz1->copy();
        f->xz2 = xz2->copy();
        f->yz1 = yz1->copy();
        f->yz2 = yz2->copy();
        
        f->Fi = Fi;
        f->Ro = Ro;
        return f;
    }

  }

  uint64_t weight() {
    return 1;
  }
};

class FragmentTools: public ts::type::FragmentTools {
friend class Fragment;
private:
public:
  FragmentTools() {
  }

  ~FragmentTools() {}

  void bserialize(ts::type::Fragment* fragment, ts::Arc* arc) {
    Fragment* f = (Fragment*) fragment;
    f->serializeBoundary(arc);
  }

  ts::type::Fragment* bdeserialize(ts::Arc* arc) {
    Fragment* result = Fragment::deserializeBoundary(arc);
    return result;
  }

  void fserialize(ts::type::Fragment* fragment, ts::Arc* arc) {
    Fragment* f = (Fragment*) fragment;
    f->serialize(arc);
  }

  ts::type::Fragment* fdeserialize(ts::Arc* arc) {
    Fragment* result = Fragment::deserialize(arc);
    return result;
  }

  ts::type::Fragment* createGap(const ID& id) override {
    return new Fragment(id);
  }
};


