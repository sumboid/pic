#ifndef PARTICLE
#define PARTICLE

#include <iostream>

class Particle {
private:
    double _x, _y, _z;
    double _vx, _vy, _vz;

    Particle() {}
public:
    Particle(double __x, double __y, double __z)
    {
        _x = __x;
        _y = __y;
        _z = __z;
        _vx = 0;
        _vy = 0;
        _vz = 0;
    }

    ~Particle() {
    }

    double x() { return _x; }
    double y() { return _y; }
    double z() { return _z; }
    double vx() { return _vx; }
    double vy() { return _vy; }
    double vz() { return _vz; }

    void x(double __x) { _x = __x; }
    void y(double __y) { _y = __y; }
    void z(double __z) { _z = __z; }
    void vx(double __vx) { _vx = __vx; }
    void vy(double __vy) { _vy = __vy; }
    void vz(double __vz) { _vz = __vz; }

    Particle* copy() {
        Particle* r = new Particle(_x, _y, _z);
        r->vx(_vx);
        r->vy(_vy);
        r->vz(_vz);
        return r;
    }

    void serialize(ts::Arc* arc) {
      ts::Arc& a = *arc;
      a << _x << _y << _z << _vx << _vy << _vz;
    }

    static Particle* deserialize(ts::Arc* arc) {
      ts::Arc& a = *arc;
      double x, y, z, vx, vy, vz;
      a >> x;
      a >> y;
      a >> z;
      a >> vx;
      a >> vy;
      a >> vz;
      Particle* r = new Particle(x, y, z);
      r->vx(vx);
      r->vy(vy);
      r->vz(vz);
      return r;
    }
};

#endif // PARTICLE

