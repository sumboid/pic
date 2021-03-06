#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <ts-ng/util/Arc.h>
#include <vector>
#include <algorithm>

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
        double tau;

        bool side[6]; //xy1 xy2 xz1 xz2 yz1 yz2
        bool corners[12]; //xz1xy1 xz1yz1 xz1xy2 xz1yz2
                          //xz2xy1 xz2yz1 xz2xy2 xz2yz2
                          //xy1yz1 xy1yz2 xy2yz1 xy2yz2

        int64_t id[3];

        std::vector<Particle*> ps;

        int element(int x, int y, int z) {
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

            int forces_size = (x + 2) * (y + 2) * (z + 2);
            Fx = new double[forces_size];
            Fy = new double[forces_size];
            Fz = new double[forces_size];

            memset(Fi, 0, mesh_size * sizeof(double));
            memset(Ro, 1, mesh_size * sizeof(double));
            for(int i = 0; i < 6; ++i) side[i] = false;
            for(int i = 0; i < 12; ++i) corners[i] = false;

            tau = 0.01;
            am = 1.0 / 10000;
            w = 1.2;
        }

        void setSizes(double x, double y, double z, int ix, int iy, int iz) {
            hx = x / (ix - 2);
            hy = y / (iy - 2);
            hz = z / (iz - 2);

            double fakehx = 4 / (64 - 2);
            double fakehy = 4 / (64 - 2);
            double fakehz = 4 / (64 - 2);
            double hx2 = 1. / (fakehx * fakehx);
            double hy2 = 1. / (fakehy * fakehy);
            double hz2 = 1. / (fakehz * fakehz);
            coef = 0.5 * w / (hx2 + hy2 + hz2);
        }

        Mesh* split() {
            Mesh* s = copy();
            std::vector<Particle*> remove;

            uint64_t size = ps.size() / 2;

            for(uint64_t i = 0; i < size; ++i) {
                s->ps.push_back(ps[i]->copy());
                delete ps[i];
            }

            ps.erase(ps.begin(), ps.begin() + size);

            return s;
        }

        Mesh(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a >> id[0];
            a >> id[1];
            a >> id[2];
            a >> size[0];
            a >> size[1];
            a >> size[2];

            Fi = new double[size[0] * size[1] * size[2]];
            Ro = new double[size[0] * size[1] * size[2]];
            Fx = new double[size[0] * size[1] * size[2]];
            Fy = new double[size[0] * size[1] * size[2]];
            Fz = new double[size[0] * size[1] * size[2]];

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
            a >> am;
            a >> tau;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a >> Fi[i];
                a >> Ro[i];
                a >> Fx[i];
                a >> Fy[i];
                a >> Fz[i];
            }

            int particles = 0;
            a >> particles;
            if(particles == 1) {
                ULOG(error) << "HAVE PARTICLES YAY" << UEND;
                uint64_t number;
                a >> number;
                for(uint64_t i = 0; i < number; ++i) {
                    ps.push_back(Particle::deserialize(arc));
                }
            }
        }

        void serialize(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a << id[0] << id[1] << id[2];
            a << size[0] << size[1] << size[2];
            a << side[0] << side[1] << side[2] <<
                 side[3] << side[4] << side[5];

            a << hx << hy << hz << w << coef << am << tau;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a << Fi[i] << Ro[i] << Fx[i] << Fy[i] << Fz[i];
            }

            int particles = 0;
            a << particles;
        }

        void serializeWithParticles(ts::Arc* arc) {
            ts::Arc& a = *arc;
            a << id[0] << id[1] << id[2];
            a << size[0] << size[1] << size[2];
            a << side[0] << side[1] << side[2] <<
                 side[3] << side[4] << side[5];

            a << hx << hy << hz << w << coef << am << tau;

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                a << Fi[i] << Ro[i] << Fx[i] << Fy[i] << Fz[i];
            }

            int particles = 1;
            a << particles;
            uint64_t number = ps.size();
            a << number;
            for(auto &i: ps) {
                i->serialize(arc);
            }
        }

        ~Mesh() {
            delete[] Fi;
            delete[] Ro;
            delete[] Fx;
            delete[] Fy;
            delete[] Fz;
            for(auto p: ps) {
                delete p;
            }
        }

        void printRo(int iteration) {
            std::ofstream file(std::to_string(id[0]) + "-" + std::to_string(id[1]) + "-" + std::to_string(iteration) + ".density");
            for(int i = 1; i < size[0] - 1; ++i) {
                for(int j = 1; j < size[1] - 1; ++j) {
                    double sum = 0;
                    for(int k = 1; k < size[2] - 1; ++k) {
                        sum += Ro[element(i, j, k)];
                    }
                    file << sum << " ";
                }
                file << std::endl;
            }
            file.close();
        }

        void printPhi(int iteration) {
            std::ofstream file(std::to_string(id[0]) + "-" + std::to_string(iteration) + "-phi");
            for(int i = 0; i < size[0]; ++i) {
                for(int j = 0; j < size[1]; ++j) {
                    double sum = 0;
                    for(int k = 0; k < size[2]; ++k) {
                        sum += Fi[element(i, j, k)];
                    }
                    file << sum << " ";
                }
                file << std::endl;
            }
            file.close();
        }

        Mesh* copy() {
           Mesh* c = new Mesh(size[0] - 2, size[1] - 2, size[2] - 2);
           c->id[0] = id[0];
           c->id[1] = id[1];
           c->id[2] = id[2];

           c->hx = hx;
           c->hy = hy;
           c->hz = hz;

           c->w = w;
           c->coef = coef;
           c->am = am;
           c->tau = tau;

           c->setSides(side[0],
                       side[1],
                       side[2],
                       side[3],
                       side[4],
                       side[5]);

           memcpy(c->Fi, Fi, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Ro, Ro, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fx, Fx, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fy, Fy, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fz, Fz, size[0] * size[1] * size[2] * sizeof(double));

           return c;
        }

        Mesh* copyWithParticles() {
           Mesh* c = new Mesh(size[0] - 2, size[1] - 2, size[2] - 2);
           c->id[0] = id[0];
           c->id[1] = id[1];
           c->id[2] = id[2];

           c->hx = hx;
           c->hy = hy;
           c->hz = hz;

           c->w = w;
           c->coef = coef;
           c->am = am;
           c->tau = tau;

           c->setSides(side[0],
                       side[1],
                       side[2],
                       side[3],
                       side[4],
                       side[5]);

           memcpy(c->Fi, Fi, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Ro, Ro, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fx, Fx, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fy, Fy, size[0] * size[1] * size[2] * sizeof(double));
           memcpy(c->Fz, Fz, size[0] * size[1] * size[2] * sizeof(double));

           for(auto &i : ps) {
               c->ps.push_back(i->copy());
           }

           return c;
        }



        void addParticle(Particle* p) {
            ps.push_back(p);
        }

        void setID(uint64_t x, uint64_t y, uint64_t z) {
            id[0] = x;
            id[1] = y;
            id[2] = z;
        }

        void sethxyz(double _hx, double _hy, double _hz) {
            hx = _hx;
            hy = _hy;
            hz = _hz;
        }

        void setw(double _w) {
            w = _w;
        }

        void setcoef(double _coef) {
            coef = _coef;
        }

        void setSides(bool xy1, bool xy2,
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

        void setFi(int x, int y, int z, double _Fi) {
            Fi[element(x + 1, y + 1, z + 1)] = _Fi;
        }

        void setRo(int x, int y, int z, double _Ro) {
            Ro[element(x + 1, y + 1, z + 1)] = _Ro;
        }

        double getFi(int x, int y, int z) {
            return Fi[element(x + 1, y + 1, z + 1)];
        }

        double getRo(int x, int y, int z) {
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

        void setBoundaryFi(std::vector<FiBoundary*>& bs) {
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
                            Fi[element(i, j, size[2] - 1)] = b->XY1[(j - 1) * (size[0] - 2) + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[1] - 1; j++) {
                            Fi[element(i, j, 0)] = b->XY2[(j - 1) * (size[0] - 2) + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(i, size[1] - 1, j)] = b->XZ1[(j - 1) * (size[0] - 2) + (i - 1)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[0] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(i, 0, j)] = b->XZ2[(j - 1) * (size[0] - 2) + (i - 1)];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[1] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(size[0] - 1, i, j)] = b->YZ1[(j - 1) * (size[1] - 2) + (i - 1)];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 1; i < size[1] - 1; i++)
                        for(int j = 1; j < size[2] - 1; j++) {
                            Fi[element(0, i, j)] = b->YZ2[(j - 1) * (size[1] - 2) + (i - 1)];
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

        void setBoundaryParticle(std::vector<ParticleBoundary*>& bs) {
            for(auto b : bs) {
                int64_t bid[3];
                bid[0] = b->id[0];
                bid[1] = b->id[1];
                bid[2] = b->id[2];

                if(bid[0] == id[0] &&
                   bid[1] == id[1] &&
                   bid[2]  > id[2]) {
                    for(auto p : b->XY1) {
                        auto particle = p->copy();
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(auto p : b->XY2) {
                        auto particle = p->copy();
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ1) {
                        auto particle = p->copy();
                        particle->y((size[1] - 2) * hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ2) {
                        auto particle = p->copy();
                        particle->y(hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ2YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y(hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ1YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ2YZ1) {
                        std::cout << "GOT PARTICLE FROM XZ2YZ1: " << p->y() << std::endl;
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y(hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(auto p : b->XZ1YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  > id[2]) {
                    for(auto p : b->XZ1XY1) {
                        auto particle = p->copy();
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  > id[2]) {
                    for(auto p : b->XZ2XY1) {
                        auto particle = p->copy();
                        particle->y(hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  < id[2]) {
                    for(auto p : b->XZ1XY2) {
                        auto particle = p->copy();
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  < id[2]) {
                    for(auto p : b->XZ2XY2) {
                        auto particle = p->copy();
                        particle->y(hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(auto p : b->XY1YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  > id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(auto p : b->XY1YZ2) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  > id[2]) {
                    for(auto p : b->XY2YZ1) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(auto p : b->XY2YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                    for(auto p : b->XZ2XY2YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y(hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] < id[2]) {
                    for(auto p : b->XZ2XY2YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y(hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(auto p : b->XZ1XY2YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                    for(auto p : b->XZ2XY1YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y(hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(auto p : b->XZ1XY2YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z(hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] > id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                    for(auto p : b->XZ2XY1YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y(hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                    for(auto p : b->XZ1XY1YZ2) {
                        auto particle = p->copy();
                        particle->x(hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] > id[2]) {
                    for(auto p : b->XZ1XY1YZ1) {
                        auto particle = p->copy();
                        particle->x((size[0] - 2) * hx + p->x());
                        particle->y((size[1] - 2) * hy + p->y());
                        particle->z((size[2] - 2) * hz + p->z());
                        ps.push_back(particle);
                    }
                }
            }
        }


        void setBoundaryRo(std::vector<RoBoundary*>& bs) {
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
                            Ro[element(i, j, 0)] += b->XY2[j * size[0] + i];
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
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[j * size[1] + i];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + j * size[1] + i];
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
                            Ro[element(0, i, j)] += b->YZ2[j * size[1] + (size[1] - 2 + i)];
                            Ro[element(1, i, j)] += b->YZ2[(size[1] * size[2]) + j * size[1] + (size[1] - 2 + i)];
                        }
                }

                else if(bid[0]  < id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(0, i, j)] += b->YZ2[j * size[1] + (i - size[1] + 2)];
                            Ro[element(1, i, j)] += b->YZ2[(size[1] * size[2]) + j * size[1] + (i - size[1] + 2)];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  < id[1] &&
                        bid[2] == id[2]) {
                    for(int i = 0; i <= 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[j * size[1] + (size[1] - 2 + i)];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + j * size[1] + (size[1] - 2 + i)];
                        }
                }

                else if(bid[0]  > id[0] &&
                        bid[1]  > id[1] &&
                        bid[2] == id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[2]; ++j) {
                            Ro[element(size[0] - 1, i, j)] += b->YZ1[j * size[1] + (i - size[1] + 2)];
                            Ro[element(size[0] - 2, i, j)] += b->YZ1[(size[1] * size[2]) + j * size[1] + (i - size[1] + 2)];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[0]; ++j) {
                            Ro[element(j, i, size[2] - 1)] += b->XY1[(i - size[1] + 2) * size[0] + j];
                            Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + j];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                        for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                            for(int j = 0; j < size[0]; ++j) {
                                Ro[element(j, i, size[2] - 1)] += b->XY1[(i - size[1] + 2) * size[0] + j];
                                Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + j];
                            }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  > id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i)
                        for(int j = 0; j < size[0]; ++j) {
                            Ro[element(j, i, 0)] += b->XY2[(i - size[1] + 2) * size[0] + j];
                            Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + j];
                        }
                }

                else if(bid[0] == id[0] &&
                        bid[1]  < id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = 0; j < size[0]; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 2 + i) * size[0] + j];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + j];
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
                           Ro[element(j, i, size[2] - 1)] += b->XY1[i * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + i * size[0] + (size[0] - 2 + j)];
                       }
                }

                else if(bid[0]  < id[0] &&
                        bid[1] == id[1] &&
                        bid[2]  < id[2]) {
                    for(int i = 0; i < size[1]; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, 0)] += b->XY2[i * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + i * size[0] + (size[0] - 2 + j)];
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
                           Ro[element(j, i, 0)] += b->XY2[(size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (j - size[0] + 2)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i) //y
                       for(int j = 0; j <= 1; ++j) { //x
                           Ro[element(j, i, 0)] += b->XY2[(i - size[1] + 2) * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + (size[0] - 2 + j)];
                       }
                }

                else if(bid[0] < id[0] &&
                        bid[1] < id[1] &&
                        bid[2] > id[2]) {
                    for(int i = 0; i <= 1; ++i)
                       for(int j = 0; j <= 1; ++j) {
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (size[1] - 2 + i) * size[0] + (size[0] - 2 + j)];
                       }
                }

                else if(bid[0] > id[0] &&
                        bid[1] > id[1] &&
                        bid[2] < id[2]) {
                    for(int i = size[1] - 2; i <= size[1] - 1; ++i) //y
                       for(int j = size[0] - 2; j <= size[0] - 1; ++j) { //x
                           Ro[element(j, i, 0)] += b->XY2[(i - size[1] + 2) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, 1)] += b->XY2[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + (j - size[0] + 2)];
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
                           Ro[element(j, i, size[2] - 1)] += b->XY1[(i - size[1] + 2) * size[0] + (j - size[0] + 2)];
                           Ro[element(j, i, size[2] - 2)] += b->XY1[(size[1] * size[0]) + (i - size[1] + 2) * size[0] + (j - size[0] + 2)];
                       }
                }
            }
        }

        double processPotential() {
            double max = 0;
            int ib = (side[4] ? 2 : 1);
            int ie = (side[5] ? size[0] - 2 : size[0] - 1);
            int jb = (side[2] ? 2 : 1);
            int je = (side[3] ? size[1] - 2 : size[1] - 1);
            int kb = (side[0] ? 2 : 1);
            int ke = (side[1] ? size[2] - 2 : size[2] - 1);

            double hx2 = 1. / (hx * hx);
            double hy2 = 1. / (hy * hy);
            double hz2 = 1. / (hz * hz);

            for(int i = ib; i < ie; i++)
                for(int j = jb; j < je; j++)
                    for(int k = kb; k < ke; k++) {
                        double Fi0 = Fi[element(i, j, k)];
                        double Fi1 = Fi[element(i - 1, j, k)];
                        double Fi2 = Fi[element(i + 1, j, k)];
                        double Fi3 = Fi[element(i, j - 1, k)];
                        double Fi4 = Fi[element(i, j + 1, k)];
                        double Fi5 = Fi[element(i, j, k - 1)];
                        double Fi6 = Fi[element(i, j, k + 1)];

                        Fi[element(i, j, k)] =  coef * ((Fi1 + Fi2) * hx2 +
                                               (Fi3 + Fi4) * hy2 +
                                               (Fi5 + Fi6) * hz2 -
                                               4 * M_PI * getRo(i, j, k)) -
                                               (1 - w) * Fi0;

                        double nmax = fabs(Fi0 - Fi[element(i, j, k)]);
                        if(max < nmax)
                            max = nmax;
                    }
            return max;
        }


        void processForces() {
            int ib = (side[4] ? 1 : 0);
            int ie = (side[5] ? size[0] - 1 : size[0]);
            int jb = (side[2] ? 1 : 0);
            int je = (side[3] ? size[1] - 1 : size[1]);
            int kb = (side[0] ? 1 : 0);
            int ke = (side[1] ? size[2] - 1 : size[2]);


            for(int i = ib; i < ie; i++)
                for(int j = jb; j < je; j++)
                    for(int k = kb; k < ke; k++) {
                        Fx[element(i, j, k)] = (Fi[element(i - 1, j, k)] - Fi[element(i, j, k)])/hx;
                        Fy[element(i, j, k)] = (Fi[element(i, j - 1, k)] - Fi[element(i, j, k)])/hy;
                        Fz[element(i, j, k)] = (Fi[element(i, j, k - 1)] - Fi[element(i, j, k)])/hz;
                    }
        }

        void uberprocessForces() {
            int ib = (side[4] ? 1 : 0);
            int ie = (side[5] ? size[0] - 1 : size[0]);
            int jb = (side[2] ? 1 : 0);
            int je = (side[3] ? size[1] - 1 : size[1]);
            int kb = (side[0] ? 1 : 0);
            int ke = (side[1] ? size[2] - 1 : size[2]);


            for(int i = ib; i < ie; i++)
                for(int j = jb; j < je; j++)
                    for(int k = kb; k < ke; k++) {
                        Fx[element(i, j, k)] = (Ro[element(i - 1, j, k)] - Ro[element(i, j, k)])/hx;
                        Fy[element(i, j, k)] = (Ro[element(i, j - 1, k)] - Ro[element(i, j, k)])/hy;
                        Fz[element(i, j, k)] = (Ro[element(i, j, k - 1)] - Ro[element(i, j, k)])/hz;
                    }
        }

        void processDensity() {
            double sx, sy, sz;
            int i, j, k;
            for(int _ = 0; _ < size[0]*size[1]*size[2]; ++_) Ro[_] = 0;

            std::ofstream out(std::to_string(id[0]) + "-density");
            for(auto &p : ps) {
                sx = p->x()/hx - 0.5; i = sx; sx = sx - i;
                sy = p->y()/hy - 0.5; j = sy; sy = sy - j;
                sz = p->z()/hz - 0.5; k = sz; sz = sz - k;

                out << "(" << i << ", " << j << ", " << k << ")" << std::endl;

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

        void fillVelocities() {
            for(auto &p : ps) {
                double xa,ya,za,xb,yb,zb;
                int ia,ka,la,ib,kb,lb;
                double fx, fy, fz;
                double x,y,z,u,v,w;

                x=p->x();
                y=p->y();
                z=p->z();
                u=p->vx();
                v=p->vy();
                w=p->vz();


                xb=x/hx;
                yb=y/hy;
                zb=z/hz;
                xa=xb-0.5;
                ya=yb-0.5;
                za=zb-0.5;
                ib=xb; xb=xb-ib;
                kb=yb; yb=yb-kb;
                lb=zb; zb=zb-lb;
                ia=xa; xa=xa-ia;
                ka=ya; ya=ya-ka;
                la=za; za=za-la;


                fx=(1-xb)*((1-ya)*((1-za)*Fx[element(ib, ka, la)]+za*Fx[element(ib, ka, la+1)])+
                               ya *((1-za)*Fx[element(ib, ka+1, la)]+za*Fx[element(ib, ka+1, la+1)]))+
                       xb *((1-ya)*((1-za)*Fx[element(ib+1, ka, la)]+za*Fx[element(ib+1, ka, la+1)])+
                               ya *((1-za)*Fx[element(ib+1, ka+1, la)]+za*Fx[element(ib+1, ka+1, la+1)]));

                fy=(1-xa)*((1-yb)*((1-za)*Fy[element(ia, kb, la)]+za*Fy[element(ia, kb, la+1)])+
                               yb *((1-za)*Fy[element(ia, kb+1, la)]+za*Fy[element(ia, kb+1, la+1)]))+
                       xa *((1-yb)*((1-za)*Fy[element(ia+1, kb, la)]+za*Fy[element(ia+1, kb, la+1)])+
                               yb *((1-za)*Fy[element(ia+1, kb+1, la)]+za*Fy[element(ia+1, kb+1, la+1)]));

                fz=(1-xa)*((1-ya)*((1-zb)*Fz[element(ia, ka, lb)]+zb*Fz[element(ia, ka, lb+1)])+
                               ya *((1-zb)*Fz[element(ia, ka+1, lb)]+zb*Fz[element(ia, ka+1, lb+1)]))+
                       xa *((1-ya)*((1-zb)*Fz[element(ia+1, ka, lb)]+zb*Fz[element(ia+1, ka, lb+1)])+
                               ya *((1-zb)*Fz[element(ia+1, ka+1, lb)]+zb*Fz[element(ia+1, ka+1, lb+1)]));

                u += tau*fx;
                v += tau*fy;
                w += tau*fz;

                p->vx(-u * 100);
                p->vy(-v * 100);
                p->vz(-w * 100);
            }
        }

        ParticleBoundary* moveParticles(int iteration) {
            ParticleBoundary* boundary = new ParticleBoundary();
            std::vector<Particle*> remove;
            std::ofstream out(std::to_string(id[0]) + std::to_string(iteration) + "-particles");

            for(std::vector<Particle*>::iterator p = ps.begin(); p != ps.end(); ++p) {
                double xa,ya,za,xb,yb,zb;
                int ia,ka,la,ib,kb,lb;
                double fx, fy, fz;
                double x,y,z,u,v,w;

                double du,dv,dw;
                const double hxt = hx/tau;
                const double hyt = hy/tau;
                const double hzt = hz/tau;

                x=(*p)->x();
                y=(*p)->y();
                z=(*p)->z();
                u=(*p)->vx();
                v=(*p)->vy();
                w=(*p)->vz();


                xb=x/hx;
                yb=y/hy;
                zb=z/hz;
                xa=xb-0.5;
                ya=yb-0.5;
                za=zb-0.5;
                ib=xb; xb=xb-ib;
                kb=yb; yb=yb-kb;
                lb=zb; zb=zb-lb;
                ia=xa; xa=xa-ia;
                ka=ya; ya=ya-ka;
                la=za; za=za-la;


                fx=(1-xb)*((1-ya)*((1-za)*Fx[element(ib, ka, la)]+za*Fx[element(ib, ka, la+1)])+
                               ya *((1-za)*Fx[element(ib, ka+1, la)]+za*Fx[element(ib, ka+1, la+1)]))+
                       xb *((1-ya)*((1-za)*Fx[element(ib+1, ka, la)]+za*Fx[element(ib+1, ka, la+1)])+
                               ya *((1-za)*Fx[element(ib+1, ka+1, la)]+za*Fx[element(ib+1, ka+1, la+1)]));

                fy=(1-xa)*((1-yb)*((1-za)*Fy[element(ia, kb, la)]+za*Fy[element(ia, kb, la+1)])+
                               yb *((1-za)*Fy[element(ia, kb+1, la)]+za*Fy[element(ia, kb+1, la+1)]))+
                       xa *((1-yb)*((1-za)*Fy[element(ia+1, kb, la)]+za*Fy[element(ia+1, kb, la+1)])+
                               yb *((1-za)*Fy[element(ia+1, kb+1, la)]+za*Fy[element(ia+1, kb+1, la+1)]));

                fz=(1-xa)*((1-ya)*((1-zb)*Fz[element(ia, ka, lb)]+zb*Fz[element(ia, ka, lb+1)])+
                               ya *((1-zb)*Fz[element(ia, ka+1, lb)]+zb*Fz[element(ia, ka+1, lb+1)]))+
                       xa *((1-ya)*((1-zb)*Fz[element(ia+1, ka, lb)]+zb*Fz[element(ia+1, ka, lb+1)])+
                               ya *((1-zb)*Fz[element(ia+1, ka+1, lb)]+zb*Fz[element(ia+1, ka+1, lb+1)]));

                du=tau*fx; if (fabs(du)<=hxt) u+=du; else { u+=(1-2*std::signbit(du))*hxt;  }
                dv=tau*fy; if (fabs(dv)<=hyt) v+=dv; else { v+=(1-2*std::signbit(dv))*hyt;  }
                dw=tau*fz; if (fabs(dw)<=hzt) w+=dw; else { w+=(1-2*std::signbit(dw))*hzt;  }

                x+=tau*u;
                y+=tau*v;
                z+=tau*w;

                int ix = x / hx;
                int iy = y / hy;
                int iz = z / hz;

                enum Side {
                    XY1 = 0,
                    XY2 = 1,
                    XZ1 = 2,
                    XZ2 = 3,
                    YZ1 = 4,
                    YZ2 = 5
                };

                bool s[6]; //xy1 xy2 xz1 xz2 yz1 yz2
                for(int i = 0; i < 6; ++i) s[i] = false;

                if (ix <= 0) {
                    if(side[4]) {
                        if(x < 1) { u=-u; x=-x; }
                    } else {
                       s[YZ1] = true;
                    }
                }
                if (iy <= 0) {
                    if(side[2]) {
                       if(y < 1) { v=-v; y=-y; }
                    } else {
                       s[XZ1] = true;
                    }
                }
                if (iz <= 0) {
                    if(side[0]) {
                       if(z < 1) { z=-z; w=-w; }
                    } else {
                       s[XY1] = true;
                    }
                }

                if (ix >= size[0]) {
                    if(side[5]) {
                       if(x > (size[0] - 1) * hx) { x = 2 * size[0] * hx - x; u=-u; }
                    } else {
                       x = x - (size[0] - 1) * hx;
                       s[YZ2] = true;
                    }
                }
                if (iy >= size[1]) {
                    if(side[3]) {
                       if(y > (size[1] - 1) * hy) { y = 2 * size[1] * hy - y; v=-v; }
                    } else {
                       y = y - (size[1] - 1) * hy;
                       s[XZ2] = true;
                    }
                }
                if (iz >= size[2]) {
                    if(side[1]) {
                        if(z > (size[2] - 1) * hz) { z = 2 * size[2] * hz - z; w=-w; }
                    } else {
                        z = z - (size[2] - 1) * hz;
                        s[XY2] = true;
                    }
                }

                out << "(" << (*p)->x() << ", " << (*p)->y() << ", " << (*p)->z() << ") -> " << "(" << x << ", " << y << ", " << z << ")" << std::endl;

                (*p)->x(x);
                (*p)->y(y);
                (*p)->z(z);
                (*p)->vx(u);
                (*p)->vy(v);
                (*p)->vz(w);

                if(s[XZ1] && s[XY1] && s[YZ1]) {
                    boundary->XZ1XY1YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[XY2] && s[YZ1]) {
                    boundary->XZ1XY2YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY1] && s[YZ1]) {
                    boundary->XZ2XY1YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY2] && s[YZ1]) {
                    boundary->XZ1XY1YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[XY1] && s[YZ2]) {
                    boundary->XZ1XY1YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[XY2] && s[YZ2]) {
                    boundary->XZ1XY2YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY1] && s[YZ2]) {
                    boundary->XZ2XY1YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY2] && s[YZ2]) {
                    boundary->XZ2XY2YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[XY1]) {
                    boundary->XZ1XY1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[XY2]) {
                    boundary->XZ1XY2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[YZ1]) {
                    boundary->XZ1YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1] && s[YZ2]) {
                    boundary->XZ1YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY1]) {
                    boundary->XZ2XY1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[XY2]) {
                    boundary->XZ2XY2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[YZ1]) {
                    std::cout << "PUT PARTICLE TO XZ2YZ1: " << std::endl;
                    boundary->XZ2YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2] && s[YZ2]) {
                    boundary->XZ2YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY1] && s[YZ1]) {
                    boundary->XY1YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY1] && s[YZ2]) {
                    boundary->XY1YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY2] && s[YZ1]) {
                    boundary->XY2YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY2] && s[YZ2]) {
                    boundary->XY2YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY1]) {
                    boundary->XY1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XY2]) {
                    boundary->XY2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ1]) {
                    boundary->XZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[XZ2]) {
                    boundary->XZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[YZ1]) {
                    boundary->YZ1.push_back((*p)->copy());
                    remove.push_back((*p));
                }
                else if(s[YZ2]) {
                    boundary->YZ2.push_back((*p)->copy());
                    remove.push_back((*p));
                }
            }

            for(auto r : remove) {
                auto it = std::find(ps.begin(), ps.end(), r);
                if(it != ps.end()) {
                    delete *it;
                    ps.erase(it);
                }
            }

            boundary->id[0] = id[0];
            boundary->id[1] = id[1];
            boundary->id[2] = id[2];

            out.close();
            return boundary;

        }

        uint64_t particlesNumber() {
            uint64_t counter = 0;
            for(auto it = ps.begin(); it != ps.end(); ++it) counter++;
            return counter;
        }

        FiBoundary* getFiBoundary() {
            FiBoundary* result = new FiBoundary(size[0] - 2, size[1] - 2, size[2] - 2);

            for(int i = 1; i < size[0] - 1; ++i)
                for(int j = 1; j < size[1] - 1; ++j)  {
                    result->XY1[(j - 1) * (size[0] - 2) + (i - 1)] = Fi[element(i, j, 1)];
                    result->XY2[(j - 1) * (size[0] - 2) + (i - 1)] = Fi[element(i, j, size[2] - 2)];
                }

            for(int i = 1; i < size[0] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    result->XZ1[(j - 1) * (size[0] - 2) + (i - 1)] = Fi[element(i, 1, j)];
                    result->XZ2[(j - 1) * (size[0] - 2) + (i - 1)] = Fi[element(i, size[1] - 2, j)];
                }

            for(int i = 1; i < size[1] - 1; ++i)
                for(int j = 1; j < size[2] - 1; ++j) {
                    result->YZ1[(j - 1) * (size[1] - 2) + (i - 1)] = Fi[element(1, i, j)];
                    result->YZ2[(j - 1) * (size[1] - 2) + (i - 1)] = Fi[element(size[0] - 2, i, j)];
                }

            for(int i = 1; i < size[0] - 1; ++i) {
                result->XZ1XY1[i - 1] = Fi[element(i, 1, 1)];
                result->XZ1XY2[i - 1] = Fi[element(i, 1, size[2] - 2)];
                result->XZ2XY1[i - 1] = Fi[element(i, size[1] - 2, 1)];
                result->XZ2XY2[i - 1] = Fi[element(i, size[1] - 2, size[2] - 2)];
            }

            for(int i = 1; i < size[1] - 1; ++i) {
                result->XY1YZ1[i - 1] = Fi[element(1, i, 1)];
                result->XY1YZ2[i - 1] = Fi[element(size[0] - 2, i, 1)];
                result->XY2YZ1[i - 1] = Fi[element(1, i, size[2] - 2)];
                result->XY2YZ2[i - 1] = Fi[element(size[0] - 2, i, size[2] - 2)];
            }

            for(int i = 1; i < size[2] - 1; ++i) {
                result->XZ1YZ1[i - 1] = Fi[element(1, 1, i)];
                result->XZ1YZ2[i - 1] = Fi[element(size[0] - 2, 1, i)];
                result->XZ2YZ1[i - 1] = Fi[element(1, size[1] - 2, i)];
                result->XZ2YZ2[i - 1] = Fi[element(size[0] - 2, size[1] - 2, i)];
            }

            result->XZ1XY1YZ1 = Fi[element(1, 1, 1)];
            result->XZ1XY2YZ1 = Fi[element(1, 1, size[2] - 2)];
            result->XZ2XY1YZ1 = Fi[element(1, size[1] - 2, 1)];
            result->XZ2XY2YZ1 = Fi[element(1, size[1] - 2, size[2] - 2)];
            result->XZ1XY1YZ2 = Fi[element(size[0] - 2, 1, 1)];
            result->XZ1XY2YZ2 = Fi[element(size[0] - 2, 1, size[2] - 2)];
            result->XZ2XY1YZ2 = Fi[element(size[0] - 2, size[1] - 2, 1)];
            result->XZ2XY2YZ2 = Fi[element(size[0] - 2, size[1] - 2, size[2] - 2)];

            result->id[0] = id[0];
            result->id[1] = id[1];
            result->id[2] = id[2];

            return result;
        }

        RoBoundary* getRoBoundary() {
            RoBoundary* result = new RoBoundary(size[0], size[1], size[2]);

            for(int i = 0; i < size[0]; ++i)
                for(int j = 0; j < size[1]; ++j)  {
                    result->XY1[j * size[0] + i] = Ro[element(i, j, 1)];
                    result->XY1[(size[0] * size[1]) + j * size[0] + i] = Ro[element(i, j, 0)];
                    result->XY2[j * size[0] + i] = Ro[element(i, j, size[2] - 2)];
                    result->XY2[(size[0] * size[1]) + j * size[0] + i] = Ro[element(i, j, size[2] - 1)];
                }

            for(int i = 0; i < size[0]; ++i)
                for(int j = 0; j < size[2]; ++j) {
                    result->XZ1[j * size[0] + i] = Ro[element(i, 1, j)];
                    result->XZ1[(size[0] * size[2]) + j * size[0] + i] = Ro[element(i, 0, j)];
                    result->XZ2[j * size[0] + i] = Ro[element(i, size[1] - 2, j)];
                    result->XZ2[(size[0] * size[2]) + j * size[0] + i] = Ro[element(i, size[1] - 1, j)];
                }

            for(int i = 0; i < size[1]; ++i)
                for(int j = 0; j < size[2]; ++j) {
                    result->YZ1[j * size[1] + i] = Ro[element(1, i, j)];
                    result->YZ1[(size[1] * size[2]) + j * size[1] + i] = Ro[element(0, i, j)];
                    result->YZ2[j * size[1] + i] = Ro[element(size[0] - 2, i, j)];
                    result->YZ2[(size[1] * size[2]) + j * size[1] + i] = Ro[element(size[0] - 1, i, j)];
                }

            result->id[0] = id[0];
            result->id[1] = id[1];
            result->id[2] = id[2];
            return result;
        }

        void updateRoM(double* aRo, bool clear) {
            if(clear) {
                for(int i = 0; i < size[0] * size[1] * size[2]; ++i){
                    Ro[i] = aRo[i];
                }
            } else {
                for(int i = 0; i < size[0] * size[1] * size[2]; ++i){
                    Ro[i] += aRo[i];
                }
            }
        }

        void updatePhiM(double* aPhi) {
            for(int i = 0; i < size[0] * size[1] * size[2]; ++i){
                Fi[i] = aPhi[i];
            }
        }

        double* getRoM() {
           // ULOG(error) << "SIZE OF ROM = " << rsize() << UEND;

            double* r = new double[size[0] * size[1] * size[2]];

            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                r[i] = Ro[i];
            }
            return r;
        }

        double* getFiM() {
            double* r = new double[size[0] * size[1] * size[2]];
            for(int i = 0; i < size[0] * size[1] * size[2]; ++i) {
                r[i] = Fi[i];
            }
            return r;
        }

        uint64_t rsize() {
            return size[0] * size[1] * size[2];
        }
};

