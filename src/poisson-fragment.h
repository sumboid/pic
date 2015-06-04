#pragma once

#include <ts-ng/types/Fragment.h>
#include <ts-ng/types/FragmentTools.h>
#include <ts-ng/types/ID.h>
#include <ts-ng/types/ReduceData.h>
#include <ts-ng/types/ReduceDataTools.h>
#include <ts-ng/util/Uberlogger.h>
#include <ts-ng/util/Arc.h>
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

  double* RoMesh;
  uint64_t msize;

  bool fib;
  bool rob;
  bool pab;

  bool rom;

  enum State {
      MOVE = 0,
      LU_MOVE = 1,
      GU_MOVE = 2,

      DENSITY = 3,
      MU_DENSITY = 4,
      GU_DENSITY = 5,
      SU_DENSITY = 6,

      POTENTIAL = 7,
      SU_POTENTIAL = 8,

      FORCE = 9
  };

  int potential_c;

  State state;

  void next() {
      int cur = static_cast<int>(state);
      if(cur == 9) state = static_cast<State>(0);
      else state = static_cast<State>(state + 1);
  }


public:
  Fragment(ts::type::ID id): ts::type::Fragment(id), max(100) {
    fib = false;
    rob = false;
    pab = false;
    rom = false;

    state = DENSITY;
    potential_c = 0;
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
    msize = mesh->rsize();
  }


  void runStep(std::vector<ts::type::Fragment*> fs) override {
      std::cout << id().c[0] << "-" << id().c[1] << "-" << id().c[2] << "-" << id().c[3] << ": State: (" << iteration()  << ", " << progress()  << ", " << fs.size() <<  ")" << std::endl;
      std::ofstream loadfile(std::to_string(nodeID()) + "-" + std::to_string(iteration()) + "-" +
                             std::to_string(progress()) + "-" + std::to_string(id().c[0]) + "-" +
                             std::to_string(id().c[1]) + "-" + std::to_string(id().c[2]) + "-" + std::to_string(id().c[3]) + ".load");

      loadfile << weight();
      if(state == DENSITY) {
          mesh->processDensity();
          if(isVirtual()) {
              RoMesh =  mesh->getRoM();
              rom = true;
              saveState();
              setUpdate();
          }
          if(isMaster()) {

              setVNeighbours(iteration(), progress());
          }
          next();
      } else if(state == MU_DENSITY) {
          if(isMaster()) {
              for(auto f: fs) {
                  ULOG(error) << "I'm master" << UEND;
                  Fragment* rf = reinterpret_cast<Fragment*>(f);
                  mesh->updateRoM(rf->RoMesh, false);
              }
          }
          if(!isVirtual()) {
              Ro = mesh->getRoBoundary();
              rob = true;
              saveState();
              setUpdate();
              setNeighbours(iteration(), progress());
          }
          next();
      } else if(state == GU_DENSITY) {
          if(!isVirtual()) {

              std::vector<RoBoundary*> bs;
              for(auto f : fs) { //it's ok no neighbours here
                  Fragment* rf = reinterpret_cast<Fragment*>(f);
                  bs.push_back(rf->Ro);
              }
              mesh->setBoundaryRo(bs);
              mesh->printRo(iteration());

          }
          if(isMaster()) {
              RoMesh = mesh->getRoM();
              rom = true;
              saveState();
              setVUpdate();
          }
          if(isVirtual()) {
              setNeighbours(iteration(), progress());
          }
          next();
      } else if(state == SU_DENSITY) {
          if(isVirtual()) {
              for(auto f: fs) {
                  Fragment* rf = reinterpret_cast<Fragment*>(f);
                  mesh->updateRoM(rf->RoMesh, true);
              }
          }
          next();
      } else if(state == POTENTIAL) {
          next();
      } else if(state == SU_POTENTIAL) {
          next();
      } else if(state == FORCE) {
          mesh->uberprocessForces();
          next();
      } else if(state == MOVE) {
          particles = mesh->moveParticles(iteration());
          if(isVirtual()) {
              pab = true;
              saveState();
              setUpdate();
          }
          if(isMaster()) {
              setVNeighbours(iteration(), progress());
          }
          next();
      } else if(state == LU_MOVE) {
          if(isMaster()) {
              for(auto f: fs) {
                  Fragment* rf = reinterpret_cast<Fragment*>(f);
                  particles->merge(rf->particles);
              }
          }
          if(!isVirtual()) {
              pab = true;
              saveState();
              setUpdate();
              setNeighbours(iteration(), progress());
          }
          next();
      } else if(state == GU_MOVE) {
          if(!isVirtual()) {
             std::vector<ParticleBoundary*> bs;
             for(auto f: fs) {
                 Fragment* rf = reinterpret_cast<Fragment*>(f);
                 bs.push_back(rf->particles);
             }
             mesh->setBoundaryParticle(bs);
          }
          next();
          nextIteration();
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
  }

  Fragment* getBoundary() override {
    Fragment* fragment = new Fragment(id());
    fragment->setBoundary();
    fragment->msize = msize;

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
    if(rom == true) {
        fragment->rom = true;
        fragment->RoMesh = new double[msize];
        memcpy(fragment->RoMesh, RoMesh, msize * sizeof(double));

        rom = false;
        delete[] RoMesh;
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
    a << rom;
    a << msize;

    if(fib == true) {
        Fi->serialize(arc);
    }
    else if(rob == true) {
        Ro->serialize(arc);
    }
    else if(pab == true) {
        particles->serialize(arc);
    }
    else if(rom == true) {
        for(uint64_t i = 0; i < msize; ++i)
            a << RoMesh[i];
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
    a >> f->rom;
    a >> f->msize;

    if(f->fib == true) {
        f->Fi = FiBoundary::deserialize(arc);
    }
    else if(f->rob == true) {
        f->Ro = RoBoundary::deserialize(arc);
    }
    else if(f->pab == true) {
        f->particles = ParticleBoundary::deserialize(arc);
    }
    else if(f->rom == true) {
        f->RoMesh = new double[f->msize];

        for(uint64_t i = 0; i < f->msize; ++i)
            a >> f->RoMesh[i];
    }

    return f;
  }

  void serialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    a << (int)state;
    a << fib;
    a << rob;
    a << pab;
    a << rom;
    a << msize;

    mesh->serializeWithParticles(arc);

    if(fib == true) {
        Fi->serialize(arc);
    }
    else if(rob == true) {
        Ro->serialize(arc);
    }
    else if(pab == true) {
        particles->serialize(arc);
    }
    else if(rom == true) {
        for(uint64_t i = 0; i < msize; ++i)
            a << RoMesh[i];
    }
  }

  static Fragment* deserialize(ts::Arc* arc) {
    ts::Arc& a = *arc;
    int s;
    a >> s;

    Fragment* f = new Fragment(ts::type::ID(0,0,0));
    f->state = (State) s;
    a >> f->fib;
    a >> f->rob;
    a >> f->pab;
    a >> f->rom;
    a >> f->msize;

    f->setMesh(new Mesh(arc));

    if(f->fib == true) {
        f->Fi = FiBoundary::deserialize(arc);
    }
    else if(f->rob == true) {
        f->Ro = RoBoundary::deserialize(arc);
    }
    else if(f->pab == true) {
        f->particles = ParticleBoundary::deserialize(arc);
    }
    else if(f->rom == true) {
        f->RoMesh = new double[f->msize];

        for(uint64_t i = 0; i < f->msize; ++i)
            a >> f->RoMesh[i];
    }
    return f;
  }

  Fragment* copy() override {
    if(!isBoundary()) {
        Fragment* f = new Fragment(id());
        f->setMesh(mesh->copyWithParticles());
        f->state = state;
        f->fib = fib;
        f->rob = rob;
        f->pab = pab;
        f->rom = rom;
        f->msize = msize;

        if(fib == true) {
            f->Fi = Fi->copy();
        }
        else if(rob == true) {
            f->Ro = Ro->copy();
        }
        else if(pab == true) {
            f->particles = particles->copy();
        }
        else if(rom == true) {
            f->RoMesh = new double[msize];
            memcpy(f->RoMesh, RoMesh, msize* sizeof(double));
        }
        return f;
    } else {
        Fragment* f = new Fragment(id());
        f->fib = fib;
        f->rob = rob;
        f->pab = pab;
        f->rom = rom;
        f->msize = msize;

        if(fib == true) {
            f->Fi = Fi->copy();
        }
        else if(rob == true) {
            f->Ro = Ro->copy();
        }
        else if(pab == true) {
            f->particles = particles->copy();
        }
        else if(rom == true) {
            ULOG(error) << "CREATE COPY OF BOUNDARY" << UEND;
            f->RoMesh = new double[msize];
            memcpy(f->RoMesh, RoMesh, msize* sizeof(double));
        }

        return f;
    }

  }

  uint64_t weight() {
    return mesh->particlesNumber();
  }

  Fragment* split() override {
      Fragment* f = new Fragment(id());

      f->fib = fib;
      f->rob = rob;
      f->pab = pab;
      f->rom = rom;
      f->msize = msize;

      f->setMesh(mesh->split());

      return f;
  }

  bool canSplit() override {
      return (10000. / (weight() + 1)) < 2 && !isVirtual();
  }

  void merge(ts::type::Fragment*) override {}

  bool canMove() override {
      return state == DENSITY || progress() == 0;
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


