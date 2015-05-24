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

  double max;

  //Boundary stuff
  FiBoundary* Fi;
  RoBoundary* Ro;
  ParticleBoundary* particles;

  bool fib;
  bool rob;
  bool pab;

  enum State {
      MOVE = 0,
      DENSITY = 1,
      POTENTIAL = 2,
      FORCE = 3
  };

  State state;

  void next() {
      if(state == MOVE) state = DENSITY;
      else if(state == DENSITY) state = POTENTIAL;
      else if(state == POTENTIAL) state = FORCE;
      else if(state == FORCE) state = MOVE;
  }


public:
  Fragment(ts::type::ID id): ts::type::Fragment(id), max(0) {
    fib = false;
    rob = false;
    pab = false;

    state = POTENTIAL;
  }

  ~Fragment() {
    if(isBoundary()) {
        if(fib) delete Fi;
        if(rob) delete Ro;
        if(pab) delete particles;
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

      if(state == POTENTIAL) {
          if(process() != 0) {
              std::vector<FiBoundary*> bs;
              for(auto f: fs) { //it's ok no neighbours here
                  Fragment* rf = reinterpret_cast<Fragment*> f;
                  bs.push_back(rf->Fi);
              }

              mesh->setBoundaryFi(bs);
          } else {
              std::vector<RoBoundary*> bs;
              for(auto f: fs) { //it's ok no neighbours here
                  Fragment* rf = reinterpret_cast<Fragment*> f;
                  bs.push_back(rf->Ro);
              }
              mesh->setBoundaryRo(bs);
          }


          max = mesh->processPotential();

          fib = true;
          Fi = mesh->getFiBoundary();

          saveState();
          setUpdate();
          setReduce();
          setNeighbours(iteration(), progress());
          return;
      }
      else if(state == FORCE) {
          std::vector<FiBoundary*> bs;
          for(auto f: fs) {
              Fragment* rf = reinterpret_cast<Fragment*> f;
              bs.push_back(rf->Fi);
          }

          mesh->setBoundaryFi(bs);

          mesh->processForces();
          next();
          return;
      }
      else if(state == MOVE) {
          particles = mesh->moveParticles();
          pab = true;

          saveState();
          setUpdate();
          setNeighbours(iteration(), progress());
          next();
          return;
      }
      else if(state == DENSITY) {
          std::vector<ParticlesBoundary*> bs;
          for(auto f: fs) {
              Fragment* rf = reinterpret_cast<Fragment*> f;
              bs.push_back(rf->particles);
          }

          mesh->setBoundaryParticle(bs);

          mesh->processDensity();

          rob = true;
          Ro = mesh->getRoBoundary();
          saveState();
          setUpdate();
          setNeighbours(iteration(), progress());
          next();
          nextIteration();
          return;
      }
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
    if(max < 0.001) next(); // ???
  }

  Fragment* getBoundary() override {
    Fragment* fragment = new Fragment(id());
    fragment->setBoundary();

    if(fib == true) {
        fragment->fib = true;
        fragment->Fi = Fi->copy();

        fib = false;
        delete Fi;
    }
    if(rob == true) {
        fragment->rob = true;
        fragment->Ro = Ro->copy();

        rob = false;
        delete Ro;
    }
    if(pab == true) {
        fragment->pab = true;
        fragment->particles = particles->copy();

        pab = false;
        delete particles;
    }

    return fragment;
  }

  void serializeBoundary(ts::Arc* arc) {
    ts::Arc& a = *arc;
    a << id().c[0];
    a << id().c[1];
    a << id().c[2];

    a << fib;
    a << rob;
    a << pab;

    if(fib == true) {
        Fi->serialize(arc);
    }
    else if(rob == true) {
        Ro->serialize(arc);
    }
    else if(pab == true) {
        particles->serialize(arc);
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
    a >> f->fib;
    a >> f->rob;
    a >> f->pab;

    if(fib == true) {
        f->Fi = FiBoundary::deserialize(arc);
    }
    else if(rob == true) {
        f->Ro = RoBoundary::deserialize(arc);
    }
    else if(pab == true) {
        f->particles = ParticleBoundary::deserialize(arc);
    }

    return f;
  }

  void serialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    a << id().c[0] << id().c[1] << id().c[2];
    a << (int)state;
    mesh->serialize(arc);
  }

  static Fragment* deserialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    int x, y, z, s;

    a >> x;
    a >> y;
    a >> z;
    a >> s;

    Fragment* f = new Fragment(ts::type::ID(x, y, z));
    f->state = (State) s;
    f->setMesh(new Mesh(arc));

    return f;
  }

  Fragment* copy() override {
    if(!isBoundary()) {
        Fragment* f = new Fragment(id());
        f->setMesh(mesh->copy());
        f->state = state;
        return f;
    } else {
        Fragment* f = new Fragment(id());
        f->fib = fib;
        f->rob = rob;
        f->pab = pab;

        if(fib == true) {
            f->Fi = Fi->copy();
        }
        else if(rob == true) {
            f->Ro = Ro->copy();
        }
        else if(pab == true) {
            f->particles = particles->copy();
        }

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


