#pragma once
#include <math.h>
#include <iostream>
#include <cstring>
#include <ts/util/Arc.h>
#include <vector>

#include "particle.h"
#include "boundary.h"


class Mesh {
    private:
        double* Fi;
        double* Ro;

        double* Fx;
        double* Fy;
        double* Fz;

        int size[3]; // X, Y, Z

        double hx, hy, hz;
        double w;
        double coef;
        double am;

        bool side[6]; //xy1 xy2 xz1 xz2 yz1 yz2
        bool corners[12]; //xz1xy1 xz1yz1 xz1xy2 xz1yz2
                          //xz2xy1 xz2yz1 xz2xy2 xz2yz2
                          //xy1yz1 xy1yz2 xy2yz1 xy2yz2

        int64_t id[3];

        std::vector<Particle*> ps;

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

            int forces_size = (x + 1) * (y + 1) * (z + 1);
            Fx = new double[forces_size];
            Fy = new double[forces_size];
            Fz = new double[forces_size];

            memset(Fi, 1, mesh_size * sizeof(double));
            memset(Ro, 1, mesh_size * sizeof(double));
            for(int i = 0; i < 6; ++i) side[i] = false;
            for(int i = 0; i < 12; ++i) corners[i] = false;
        }

        Mesh(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a >> size[0];
            a >> size[1];
            a >> size[2];

            Fi = new double[size[0] * size[1] * size[2]];
            Ro = new double[size[0] * size[1] * size[2]];

            a >> side[0];
            a >> side[1];
            a >> side[2];
            a >> side[3];
            a >> side[4];
            a >> side[5];

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
           c->setSides(side[0],
                         side[1],
                         side[2],
                         side[3],
                         side[4],
                         side[5]);
           memcpy(c->Fi, Fi, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Ro, Ro, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fx, Fx, (size[0] - 1) * (size[1] - 1) * (size[2] - 1) * sizeof(double));
           memcpy(c->Fy, Fy, (size[0] - 1) * (size[1] - 1) * (size[2] - 1) * sizeof(double));
           memcpy(c->Fz, Fz, (size[0] - 1) * (size[1] - 1) * (size[2] - 1) * sizeof(double));
           return c;
        }

        inline void addParticle(Particle* p) {
            ps.push_back(p);
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

        inline void setSides(bool xy1, bool xy2,
                             bool xz1, bool xz2,
                             bool yz1, bool yz2) {
            side[0] = xy1;
            side[1] = xy2;
            side[2] = xz1;
            side[3] = xz2;
            side[4] = yz1;
            side[5] = yz2;

            if(xz1 == true) {
                if(xy1 == true) corners[0] = true;
                if(yz1 == true) corners[1] = true;
                if(xy2 == true) corners[2] = true;
                if(yz2 == true) corners[3] = true;
            }
            if(xz2 == true) {
                if(xy1 == true) corners[4] = true;
                if(yz1 == true) corners[5] = true;
                if(xy2 == true) corners[6] = true;
                if(yz2 == true) corners[7] = true;
            }
            if(xy1 == true) {
                if(yz1 == true) corners[8] = true;
                if(yz2 == true) corners[9] = true;
            }
            if(xy2 == true) {
                if(yz1 == true) corners[10] = true;
                if(yz2 == true) corners[11] = true;
            }
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

        inline int getSizeX() {
            return size[0] - 2;
        }

        inline int getSizeY() {
            return size[1] - 2;
        }

        inline int getSizeZ() {
            return size[2] - 2;
        }

        inline setBoundaryFi(std::vector<FiBoundary*>& bs) {
            for(auto b : bs) {
                int64_t bid[3];
                bid[0] = b->id[0];
                bid[1] = b->id[1];
                bid[2] = b->id[2];

                if(bid[0] == id[0] &&
                   bid[1] == id[1] &&
                   bid[2]  > id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[1] - 1; j++) {
                            Fi[element(i, j, size[2] - 1)] = b->XY1[(j - 1) * size[0] + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[1] - 1; j++) {
                            Fi[element(i, j, 0)] = b->XY2[(j - 1) * size[0] + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(i, size[1] - 1, j)] = b->XZ1[(j - 1) * size[0] + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(i, 0, j)] = b->XZ2[(j - 1) * size[0] + (i - 1)];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[1] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(size[0] - 1, i, j)] = b->XY1[(j - 1) * size[1] + (i - 1)];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[1] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(0, i, j)] = b->XY2[(j - 1) * size[1] + (i - 1)];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[2] - 1; ++i) {
                        Fi[element(0, 0, i)] = b->XZ2YZ2[i - 1];
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[2] - 1; ++i) {
                        Fi[element(0, size[1] - 1, i)] = b->XZ1YZ2[i - 1];
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[2] - 1; ++i) {
                        Fi[element(size[0] - 1, 0, i)] = b->XZ2YZ1[i - 1];
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[2] - 1; ++i) {
                        Fi[element(size[0] - 1, size[1] - 1, i)] = b->XZ1YZ1[i - 1];
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 1; i < size[0] - 1; ++i) {
                        Fi[element(i, size[1] - 1, size[2] - 1)] = b->XZ1XY1[i - 1];
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 1; i < size[0] - 1; ++i) {
                        Fi[element(i, 0, size[2] - 1)] = b->XZ2XY1[i - 1];
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[0] - 1; ++i) {
                        Fi[element(i, size[1] - 1, 0)] = b->XZ1XY2[i - 1];
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[0] - 1; ++i) {
                        Fi[element(i, 0, 0)] = b->XZ2XY2[i - 1];
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 1; i < size[1] - 1; ++i) {
                        Fi[element(size[0] - 1, i, size[2] - 1)] = b->XY1YZ1[i - 1];
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[1] - 1; ++i) {
                        Fi[element(size[0] - 1, i, 0)] = b->XY1YZ2[i - 1];
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 1; i < size[1] - 1; ++i) {
                        Fi[element(0, i, size[2] - 1)] = b->XY2YZ1[i - 1];
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[1] - 1; ++i) {
                        Fi[element(0, i, 0)] = b->XY2YZ2[i - 1];
                    }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                     Fi[element(0, 0, 0)] = b->XZ2XY2YZ2;
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                     Fi[element(size[0] - 1, 0, 0)] = b->XZ2XY2YZ1;
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                     Fi[element(0, size[1] - 1, 0)] = b->XZ1XY2YZ2;
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                     Fi[element(0, 0, size[2] - 1)] = b->XZ2XY1YZ2;
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                     Fi[element(size[0] - 1, size[1] - 1, 0)] = b->XZ1XY2YZ1;
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                     Fi[element(size[0] - 1, 0, size[2] - 1)] = b->XZ2XY1YZ1;
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                     Fi[element(0, size[1] - 1, size[2] - 1)] = b->XZ1XY1YZ2;
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                     Fi[element(size[0] - 1, size[1] - 1, size[2] - 1)] = b->XZ1XY1YZ1;
                }
            }
        }

        inline setBoundaryRo(std::vector<RoBoundary*>& bs) {
            for(auto b : bs) {
                int64_t bid[3];
                bid[0] = b->id[0];
                bid[1] = b->id[1];
                bid[2] = b->id[2];

                if(bid[0] == id[0] &&
                   bid[1] == id[1] &&
                   bid[2]  > id[2]) {
                    for(int i = 0; i < size[0]; i++)
                        for(int j = 0; j < size[1]; j++) {
                            Ro[element(i, j, size[2] - 1)] += b->XY1[j * size[0] + i];
                            Ro[element(i, j, size[2] - 2)] += b->XY1[(size[0] * size[1]) + j * size[0] + i];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i < size[0]; i++)
                        for(int j = 0; j < size[1]; j++) {
                            Ro[element(i, j, 0)] += b->XY2[j * size[0]];
                            Ro[element(i, j, 1)] += b->XY2[(size[0] * size[1]) + j * size[0] + i];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int j = 0; j < size[2]; j++)
                        for(int i = 0; i < size[0]; i++) {
                            Ro[element(i, size[1] - 1, j)] += b->XZ1[j * size[0] + i];
                            Ro[element(i, size[1] - 2, j)] += b->XZ1[(size[0] * size[2]) + j * size[0] + i];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 0; i < size[0]; i++)
                        for(int j = 0; j < size[2]; j++) {
                            Ro[element(i, 0, j)] += b->XZ2[j * size[0] + i];
                            Ro[element(i, 1, j)] += b->XZ2[(size[0] * size[2]) + j * size[0] + i];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 0; i < size[1]; i++)
                        for(int j = 0; j < size[2]; j++) {
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[i * size[2] + j];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + i * size[2] + j];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 0; i < size[1]; i++)
                        for(int j = 0; j < size[2]; j++) {
                            Ro[element(0, i, j)] += b->YZ2[j * size[1] + i];
                            Ro[element(1, i, j)] += b->YZ2[(size[1] * size[2]) + j * size[1] + i];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {

                    for(int i = 0; i <= 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(0, i, j)] += b->YZ2[i * size[2] + j];
                            Ro[element(1, i, j)] += b->YZ2[(size[1] * size[2]) + i * size[2] + j];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(0, i, j)] += b->YZ2[(size[1] - i + 2) * size[2] + j];
                            Ro[element(1, i, j)] += b->YZ2[(size[1] * size[2]) + (size[1] - i + 2) * size[2] + j];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 0; i <= 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[((size[1] - 1 - i) * size[2] + j];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + (size[1] - 1 - i)  * size[2] + j];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[(size[1] - i + 2) * size[2] + j];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + (size[1] - i + 2) * size[2] + j];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[0]; ++j) {
                            Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - i + 2) * size[0] + j];
                            Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - i + 2) * size[0] + j];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                        for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                            for(int j = 0; j < size[0]; ++j) {
                                Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - 1 - i) * size[0] + j];
                                Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - 1 - i) * size[0] + j];
                            }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[0]; ++j) {
                            Ro[element(j, i, 0)] += b->XY2[(size[1] - i + 2) * size[0] + j];
                            Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - i + 2) * size[0] + j];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = 0; j < size[0]; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 1 - i) * size[0] + j];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 1 - i) * size[0] + j];
                       }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 0; i < size[1]; ++i)
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[i * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + i * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i < size[1]; ++i)
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[i * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + i * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 0; i < size[1]; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[i * size[0] + (size[0] - 1 - j)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + i * size[0] + (size[0] - 1 - j)];
                       }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i < size[1]; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[i * size[0] + (size[0] - 1 - j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + i * size[0] + (size[0] - 1 - j)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                       }
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                    for(int i = 0; i <= 1; ++i) //y
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) { //x
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 2 + i) * size[0] + (size[0] - 2 - j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (size[0] - 2 - j)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i) //y
                       for(int j = 0; j <= 1; ++j) { //x
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (size[0] - 2 - j)];
                       }
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i) //y
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) { //x
                           Ro[element(j, i, 0)] += b->XY2[(j - size[1] + 2) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (j - size[1] + 2) * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i) //y
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) { //x
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(j - size[1] + 2) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (j - size[1] + 2) * size[0] + (j - size[0] + 2)];
                       }
                }
            }
        }

        inline double process() {
            double max = 0;
            int ib = (side[4] ? 2 : 1);
            int ie = (side[5] ? size[0] - 2 : size[0] - 1);
            int jb = (side[2] ? 2 : 1);
            int je = (side[3] ? size[1] - 2 : size[1] - 1);
            int kb = (side[0] ? 2 : 1);
            int ke = (side[1] ? size[2] - 2 : size[2] - 1);

            Mesh* c = copy();
            double hx2 = 1. / (hx * hx);
            double hy2 = 1. / (hy * hy);
            double hz2 = 1. / (hz * hz);

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

                        setFi(i, j, k, coef * ((Fi1 + Fi2) * hx2 +
                                               (Fi3 + Fi4) * hy2 +
                                               (Fi5 + Fi6) * hz2 -
                                               4 * M_PI * getRo(i, j, k)) -
                                               (1 - w) * Fi0);

                        double nmax = fabs(Fi0 - getFi(i, j, k));
                        if(max < nmax)
                            max = nmax;
                    }
            delete c;
            return max;
        }

        void processForces() {
            int ib = (side[4] ? 2 : 1);
            int ie = (side[5] ? size[0] - 2 : size[0] - 1);
            int jb = (side[2] ? 2 : 1);
            int je = (side[3] ? size[1] - 2 : size[1] - 1);
            int kb = (side[0] ? 2 : 1);
            int ke = (side[1] ? size[2] - 2 : size[2] - 1);


            for(int i = ib; i < ie; i++)
                for(int j = jb; j < je; j++)
                    for(int k = kb; k < ke; k++) {
                        Fx[element(i, j, k)] = (Fi[element(i - 1, j, k)] - Fi[element(i, j, k)])/hx;
                        Fy[element(i, j, k)] = (Fi[element(i, j - 1, k)] - Fi[element(i, j, k)])/hy;
                        Fz[element(i, j, k)] = (Fi[element(i, j, k - 1)] - Fi[element(i, j, k)])/hz;
                    }
        }

        void processDensity() {
            double sx, sy, sz;
            int i, j, k;
            for(int _ = 0; _ < size[0]*size[1]*size[2]; ++_) Ro[_] = 0;

            for(auto p : ps) { 
                sx = p->x()/hx - 0.5; i = sx; sx = sx - i; i++;
                sy = p->y()/hy - 0.5; j = sy; sy = sy - j; j++;
                sz = p->z()/hz - 0.5; k = sz; sz = sz - k; k++;

                Ro[element(i, j, k)] += (1 - sx) * (1 - sy) * (1 - sz);
                Ro[element(i, j, k + 1)] += (1 - sx) * (1 - sy) * sz;
                Ro[element(i, j + 1, k)] += (1 - sx) * sy * (1 - sz);
                Ro[element(i, j + 1, k + 1)] += (1 - sx) * sy * sz;
                Ro[element(i + 1, j, k)] += sx * (1 - sy) * (1 - sz);
                Ro[element(i + 1, j, k + 1)] += sx * (1 - sy) * sz;
                Ro[element(i + 1, j + 1, k)] += sx * sy * (1 - sz);
                Ro[element(i + 1, j + 1, k + 1)] += sx * sy * sz;
            }

            double s = am / (hx * hy * hz);
            for(int _ = 0; _ < size[0]*size[1]*size[2]; ++_) Ro[_] *= s;
        }

        std::vector<Particle*> moveParticles() {

        }

        FiBoundary** getBoundary() {
            FiBoundary** result = new FiBoundary*[6];
            std::cout << "Heyho " << size[0] << " " << size[1] << std::endl;
            FiBoundary* xy1 = new FiBoundary(size[0] * size[1]);
            FiBoundary* xy2 = new FiBoundary(size[0] * size[1]);
            FiBoundary* xz1 = new FiBoundary(size[0] * size[2]);
            FiBoundary* xz2 = new FiBoundary(size[0] * size[2]);
            FiBoundary* yz1 = new FiBoundary(size[1] * size[2]);
            FiBoundary* yz2 = new FiBoundary(size[1] * size[2]);

            for(int i = 1; i < size[0] - 1; ++i)
                for(int j = 1; j < size[1] - 1; ++j)  {
                    xy1->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, j, 1)];
                    xy1->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, j, 1)];

                    xy2->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, j, size[2] - 2)];
                    xy2->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, j, size[2] - 2)];
                }

            for(int i = 1; i < size[0] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    xz1->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, 1, j)];
                    xz1->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, 1, j)];

                    xz2->Fi[(j - 1) * size[0] + (i - 1)] = Fi[element(i, size[1] - 2, j)];
                    xz2->Ro[(j - 1) * size[0] + (i - 1)] = Ro[element(i, size[1] - 2, j)];
                }

            for(int i = 1; i < size[1] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    yz1->Fi[(j - 1) * size[1] + (i - 1)] = Fi[element(1, i, j)];
                    yz1->Ro[(j - 1) * size[1] + (i - 1)] = Ro[element(1, i, j)];

                    yz2->Fi[(j - 1) * size[1] + (i - 1)] = Fi[element(size[0] - 2, i, j)];
                    yz2->Ro[(j - 1) * size[1] + (i - 1)] = Ro[element(size[0] - 2, i, j)];
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
            a << side[0] << side[1] << side[2] <<
                 side[3] << side[4] << side[5];

            a << hx << hy << hz << w << coef;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a << Fi[i] << Ro[i];
            }
        }
};

