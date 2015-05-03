#pragma once
#include <math.h>
#include <iostream>
#include <cstring>
#include <ts/util/Arc.h>

struct MeshBoundary {
    double* Fi;
    double* Ro;

    int _size;


    ~MeshBoundary() {
        delete[] Fi;
        delete[] Ro;
    }

    MeshBoundary(int size) {
        Fi = new double[size];
        Ro = new double[size];
        _size = size;
    }

    MeshBoundary* copy() {
        MeshBoundary* result = new MeshBoundary(_size);
        memcpy(result->Fi, Fi, _size * sizeof(double));
        memcpy(result->Ro, Ro, _size * sizeof(double));
        return result;
    }

    void serializeFi(ts::Arc* arc) {
        ts::Arc& a = *arc;
        std::cout << "Serialize " << _size << " elements" << std::endl;
        a << _size;
        for(int i = 0; i < _size; ++i) a << Fi[i];
    }

    void serializeRo(ts::Arc* arc) {
        ts::Arc& a = *arc;
        a << _size;
        for(int i = 0; i < _size; ++i) a << Ro[i];
    }

    static MeshBoundary* deserializeFi(ts::Arc* arc) {
        ts::Arc& a = *arc;
        int size;
        a >> size;
        std::cout << "Deserialize " << size << " elements" << std::endl;
        MeshBoundary* r = new MeshBoundary(size);
        for(int i = 0; i < size; ++i) a >> r->Fi[i];
        return r;
    }


    static MeshBoundary* deserializeRo(ts::Arc* arc) {
        ts::Arc& a = *arc;
        int size;
        a >> size;
        MeshBoundary* r = new MeshBoundary(size);
        for(int i = 0; i < size; ++i) a >> r->Ro[i];
        return r;
    }
};

struct Force {
  double _1;
  double _2;
};

class Mesh {
    private:
        double* Fi;
        double* Ro;
        Force* Fx;
        Force* Fy;
        Force* Fz;

        int size[3]; // X, Y, Z

        double hx, hy, hz;
        double w;
        double coef;

        bool corners[6]; //xy1 xy2 xz1 xz2 yz1 yz2
        inline int element(int x, int y, int z) {
            return size[0]*size[1]*z + size[0]*y + x;
        }




    public:
        Mesh(int x, int y, int z) {
            size[0] = x + 2;
            size[1] = y + 2;
            size[2] = z + 2;

            int mesh_size = (x + 2) * (y + 2) * (z + 2);
            Fi = new double[mesh_size];
            Ro = new double[mesh_size];

            memset(Fi, 1, mesh_size * sizeof(double));
            memset(Ro, 1, mesh_size * sizeof(double));
            for(int i = 0; i < 6; ++i) corners[i] = false;
        }

        Mesh(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a >> size[0];
            a >> size[1];
            a >> size[2];

            Fi = new double[size[0] * size[1] * size[2]];
            Ro = new double[size[0] * size[1] * size[2]];

            a >> corners[0];
            a >> corners[1];
            a >> corners[2];
            a >> corners[3];
            a >> corners[4];
            a >> corners[5];

            a >> hx;
            a >> hy;
            a >> hz;

            a >> w;
            a >> coef;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a >> Fi[i];
                a >> Ro[i];
            }
        }

        ~Mesh() {
            delete[] Fi;
            delete[] Ro;
        }

        inline Mesh* copy() {
           Mesh* c = new Mesh(size[0] - 2, size[1] - 2, size[2] - 2);
           c->setCorners(corners[0],
                         corners[1],
                         corners[2],
                         corners[3],
                         corners[4],
                         corners[5]);
           memcpy(c->Fi, Fi, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Ro, Ro, size[0] * size[1] * size[2] * sizeof(double));
           return c;
        }

        inline void sethxyz(double _hx, double _hy, double _hz) {
            hx = _hx;
            hy = _hy;
            hz = _hz;
        }

        inline void setw(double _w) {
            w = _w;
        }

        inline void setcoef(double _coef) {
            coef = _coef;
        }

        inline void setCorners(bool xy1, bool xy2, 
                               bool xz1, bool xz2,
                               bool yz1, bool yz2) {
            corners[0] = xy1;
            corners[1] = xy2;
            corners[2] = xz1;
            corners[3] = xz2;
            corners[4] = yz1;
            corners[5] = yz2;
        }

        inline void setFi(int x, int y, int z, double _Fi) {
            Fi[element(x + 1, y + 1, z + 1)] = _Fi;
        }

        inline void setRo(int x, int y, int z, double _Ro) {
            Ro[element(x + 1, y + 1, z + 1)] = _Ro;
        }

        inline double getFi(int x, int y, int z) {
            return Fi[element(x + 1, y + 1, z + 1)];
        }

        inline double getRo(int x, int y, int z) {
            return Ro[element(x + 1, y + 1, z + 1)];
        }

        int getSizeX() {
            return size[0] - 2;
        }

        int getSizeY() {
            return size[1] - 2;
        }

        int getSizeZ() {
            return size[2] - 2;
        }

        void setBoundaryXY1(MeshBoundary* boundary) {
            for(int i = 1; i < size[0] - 1; i++)
                for(int j = 1; j < size[1] - 1; j++) {
                    Fi[element(i, j, 0)] = boundary->Fi[(j - 1) * size[0] + (i - 1)];
                    Ro[element(i, j, 0)] = boundary->Ro[(j - 1) * size[0] + (i - 1)];
                }
        }

        void setBoundaryXY2(MeshBoundary* boundary) {
            for(int i = 1; i < size[0] - 1; i++)
                for(int j = 1; j < size[1] - 1; j++) {
                    Fi[element(i, j, size[2] - 1)] = boundary->Fi[(j - 1) * size[0] + (i - 1)];
                    Ro[element(i, j, size[2] - 1)] = boundary->Ro[(j - 1) * size[0] + (i - 1)];
                }
        }

        void setBoundaryYZ1(MeshBoundary* boundary) {
            for(int i = 1; i < size[1] - 1; i++)
                for(int j = 1; j < size[2] - 1; j++) {
                    Fi[element(0, i, j)] = boundary->Fi[(j - 1) * size[1] + (i - 1)];
                    Ro[element(0, i, j)] = boundary->Ro[(j - 1) * size[1] + (i - 1)];
                }
        }

        void setBoundaryYZ2(MeshBoundary* boundary) {
            for(int i = 1; i < size[1] - 1; i++)
                for(int j = 1; j < size[2] - 1; j++) {
                    Fi[element(size[0] - 1, i, j)] = boundary->Fi[(j - 1) * size[1] + (i - 1)];
                    Ro[element(size[0] - 1, i, j)] = boundary->Ro[(j - 1) * size[1] + (i - 1)];
                }
        }

        void setBoundaryXZ1(MeshBoundary* boundary) {
            for(int i = 1; i < size[0] - 1; i++)
                for(int j = 1; j < size[2] - 1; j++) {
                    Fi[element(i, 0, j)] = boundary->Fi[(j - 1) * size[0] + (i - 1)];
                    Ro[element(i, 0, j)] = boundary->Ro[(j - 1) * size[0] + (i - 1)];
                }
        }

        void setBoundaryXZ2(MeshBoundary* boundary) {
            for(int i = 1; i < size[0] - 1; i++)
                for(int j = 1; j < size[2] - 1; j++) {
                    Fi[element(i, size[1] - 1, j)] = boundary->Fi[(j - 1) * size[0] + (i - 1)];
                    Ro[element(i, size[1] - 1, j)] = boundary->Ro[(j - 1) * size[0] + (i - 1)];
                }
        }

        double process() {
            double max = 0;
            int ib = (corners[4] ? 2 : 1);
            int ie = (corners[5] ? size[0] - 2 : size[0] - 1);
            int jb = (corners[2] ? 2 : 1);
            int je = (corners[3] ? size[1] - 2 : size[1] - 1);
            int kb = (corners[0] ? 2 : 1);
            int ke = (corners[1] ? size[2] - 2 : size[2] - 1);

            Mesh* c = copy();

            for(int i = ib; i < ie; i++)
                for(int j = jb; j < je; j++)
                    for(int k = kb; k < ke; k++) {
                        double Fi0 = c->Fi[element(i, j, k)];
                        double Fi1 = c->Fi[element(i - 1, j, k)];
                        double Fi2 = c->Fi[element(i + 1, j, k)];
                        double Fi3 = c->Fi[element(i, j - 1, k)];
                        double Fi4 = c->Fi[element(i, j + 1, k)];
                        double Fi5 = c->Fi[element(i, j, k - 1)];
                        double Fi6 = c->Fi[element(i, j, k + 1)];

                        setFi(i, j, k, coef * ((Fi1 + Fi2) * hx +
                                               (Fi3 + Fi4) * hy +
                                               (Fi5 + Fi6) * hz -
                                               4 * M_PI * getRo(i, j, k)) -
                                               (1 - w) * Fi0);

                        double nmax = fabs(Fi0 - getFi(i, j, k));
                        if(max < nmax)
                            max = nmax;
                    }
            delete c;
            return max;
        }

        MeshBoundary** getBoundary() {
            MeshBoundary** result = new MeshBoundary*[6];
            std::cout << "Heyho " << size[0] << " " << size[1] << std::endl;
            MeshBoundary* xy1 = new MeshBoundary(size[0] * size[1]);
            MeshBoundary* xy2 = new MeshBoundary(size[0] * size[1]);
            MeshBoundary* xz1 = new MeshBoundary(size[0] * size[2]);
            MeshBoundary* xz2 = new MeshBoundary(size[0] * size[2]);
            MeshBoundary* yz1 = new MeshBoundary(size[1] * size[2]);
            MeshBoundary* yz2 = new MeshBoundary(size[1] * size[2]);

            for(int i = 1; i < size[0] - 1; ++i) // x
                for(int j = 1; j < size[1] - 1; ++j)  { // y
                    xy1->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, j, 0)];
                    xy1->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, j, 0)];

                   // std::cout << "ELEMENT: " << element(i, j, size[2] - 1)  << " = (" << i << ", " << j << ", " << size[0] << ")" << std::endl;
                    xy2->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, j, size[2] - 1)];
                    xy2->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, j, size[2] - 1)];
                }

            for(int i = 1; i < size[0] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    xz1->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, 0, j)];
                    xz1->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, 0, j)];

                    xz2->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, size[1] - 1, j)];
                    xz2->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, size[1] - 1, j)];
                }

            for(int i = 1; i < size[1] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    yz1->Fi[(j - 1) * size[1] + (i - 1)] = Fi[element(0, i, j)];
                    yz1->Ro[(j - 1) * size[1] + (i - 1)] = Ro[element(0, i, j)];

                    yz2->Fi[(j - 1) * size[1] + (i - 1)] = Fi[element(size[0] - 1, i, j)];
                    yz2->Ro[(j - 1) * size[1] + (i - 1)] = Ro[element(size[0] - 1, i, j)];
                }

            result[0] = xy1;
            result[1] = xy2;
            result[2] = xz1;
            result[3] = xz2;
            result[4] = yz1;
            result[5] = yz2;

            return result;
        }

        void serialize(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a << size[0] << size[1] << size[2];
            a << corners[0] << corners[1] << corners[2] <<
                 corners[3] << corners[4] << corners[5];

            a << hx << hy << hz << w << coef;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a << Fi[i] << Ro[i];
            }
        }
};

