#ifndef PARTICLE
#define PARTICLE

class Particle {
private:
    double _x, _y, _z;
    double _vx, _vy, _vz;
    double _m;

public:
    Particle(double __x, double __y, double __z, double __m): _x(__x), _y(__y), _z(__z), _vx(0), _vy(0), _vz(0), _m(__m) {}
    Particle(double __x, double __y, double __z): _x(__x), _y(__y), _z(__z), _vx(0), _vy(0), _vz(0), _m(1) {}
    Particle(): _x(0), _y(0), _z(0), _vx(0), _vy(0), _vz(0), _m(1) {}

    inline double m() { return _m; }
    inline double x() { return _x; }
    inline double y() { return _y; }
    inline double z() { return _z; }
    inline double vx() { return _vx; }
    inline double vy() { return _vy; }
    inline double vz() { return _vz; }

    inline void m(double __m) { _m = __m; }
    inline void x(double __x) { _x = __x; }
    inline void y(double __y) { _y = __y; }
    inline void z(double __z) { _z = __z; }
    inline void vx(double __vx) { _vx = __vx; }
    inline void vy(double __vy) { _vy = __vy; }
    inline void vz(double __vz) { _vz = __vz; }

    inline Particle* copy() {
        Particle* r = new Particle(_x, _y, _z, _m);
        r->vx(_vx);
        r->vy(_vy);
        r->vz(_vz);
        return r;
    }

    void serialize(ts::Arc* arc) {
      ts::Arc& a = *arc;
      a << _x << _y << _z << _vx << _vy << _vz << _m;
    }

    static Particle* deserialize(ts::Arc* arc) {
      ts::Arc& a = *arc;
      double x, y, z, vx, vy, vz, m;
      a >> x;
      a >> y;
      a >> z;
      a >> vx;
      a >> vy;
      a >> vz;
      a >> m;
      Particle* r = new Particle(x, y, z, m);
      r->vx(vx);
      r->vy(vy);
      r->vz(vz);
      return r;
    }
};

#endif // PARTICLE

