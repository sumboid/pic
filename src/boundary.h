#ifndef BOUNDARY
#define BOUNDARY

struct FiBoundary {
    double* XY1;
    double* XY2;
    double* XZ1;
    double* XZ2;
    double* YZ1;
    double* YZ2;

    double* XZ1XY1;
    double* XZ1XY2;
    double* XZ1YZ1;
    double* XZ1YZ2;
    double* XZ2XY1;
    double* XZ2XY2;
    double* XZ2YZ1;
    double* XZ2YZ2;
    double* XY1YZ1;
    double* XY1YZ2;
    double* XY2YZ1;
    double* XY2YZ2;

    double XZ1XY1YZ1;
    double XZ1XY2YZ1;
    double XZ2XZ1YZ1;
    double XZ2XY2YZ1;
    double XZ1XY1YZ2;
    double XZ1XY2YZ2;
    double XZ2XZ1YZ2;
    double XZ2XY2YZ2;

    int _size[3];
    int64_t id[3];

    ~FiBoundary() {
        delete[] XY1;
        delete[] XY2;
        delete[] XZ1;
        delete[] XZ2;
        delete[] YZ1;
        delete[] YZ2;

        delete[] XZ1XY1;
        delete[] XZ1XY2;
        delete[] XZ1YZ1;
        delete[] XZ1YZ2;
        delete[] XZ2XY1;
        delete[] XZ2XY2;
        delete[] XZ2YZ1;
        delete[] XZ2YZ2;
        delete[] XY1YZ1;
        delete[] XY1YZ2;
        delete[] XY2YZ1;
        delete[] XY2YZ2;
    }

    FiBoundary(int sizeX, int sizeY, int sizeZ) {
        XY1 = new double[sizeX * sizeY];
        XY2 = new double[sizeX * sizeY];
        XZ1 = new double[sizeX * sizeZ];
        XZ2 = new double[sizeX * sizeZ];
        YZ1 = new double[sizeY * sizeZ];
        YZ2 = new double[sizeY * sizeZ];

        XZ1XY1 = new double[sizeX];
        XZ1XY2 = new double[sizeX];
        XZ1YZ1 = new double[sizeZ];
        XZ1YZ2 = new double[sizeZ];
        XZ2XY1 = new double[sizeX];
        XZ2XY2 = new double[sizeX];
        XZ2YZ1 = new double[sizeZ];
        XZ2YZ2 = new double[sizeZ];
        XY1YZ1 = new double[sizeY];
        XY1YZ2 = new double[sizeY];
        XY2YZ1 = new double[sizeY];
        XY2YZ2 = new double[sizeY];

        size[0] = sizeX;
        size[1] = sizeY;
        size[2] = sizeZ;
    }

    FiBoundary* copy() {
        FiBoundary* result = new FiBoundary(size[0], size[1], size[2]);
        memcpy(result->XY1, XY1, size[0] * size[1] * sizeof(double));
        memcpy(result->XY2, XY2, size[0] * size[1] * sizeof(double));
        memcpy(result->XZ1, XZ1, size[0] * size[2] * sizeof(double));
        memcpy(result->XZ2, XZ2, size[0] * size[2] * sizeof(double));
        memcpy(result->YZ1, YZ1, size[1] * size[2] * sizeof(double));
        memcpy(result->YZ2, YZ2, size[1] * size[2] * sizeof(double));

        memcpy(result->XZ1XY1, XZ1XY1, size[0] * sizeof(double));
        memcpy(result->XZ1XY2, XZ1XY2, size[0] * sizeof(double));
        memcpy(result->XZ1YZ1, XZ1YZ1, size[2] * sizeof(double));
        memcpy(result->XZ1YZ2, XZ1YZ2, size[2] * sizeof(double));
        memcpy(result->XZ2XY1, XZ2XY1, size[0] * sizeof(double));
        memcpy(result->XZ2XY2, XZ2XY2, size[0] * sizeof(double));
        memcpy(result->XZ2YZ1, XZ2YZ1, size[2] * sizeof(double));
        memcpy(result->XZ2YZ2, XZ2YZ2, size[2] * sizeof(double));
        memcpy(result->XY1YZ1, XY1YZ1, size[1] * sizeof(double));
        memcpy(result->XY1YZ2, XY1YZ2, size[1] * sizeof(double));
        memcpy(result->XY2YZ1, XY2YZ1, size[1] * sizeof(double));
        memcpy(result->XY2YZ2, XY2YZ2, size[1] * sizeof(double));

        result->XZ1XY1YZ1 = XZ1XY1YZ1;
        result->XZ1XY2YZ1 = XZ1XY2YZ1;
        result->XZ2XZ1YZ1 = XZ2XZ1YZ1;
        result->XZ2XY2YZ1 = XZ2XY2YZ1;
        result->XZ1XY1YZ2 = XZ1XY1YZ2;
        result->XZ1XY2YZ2 = XZ1XY2YZ2;
        result->XZ2XZ1YZ2 = XZ2XZ1YZ2;
        result->XZ2XY2YZ2 = XZ2XY2YZ2;

        return result;
    }

    void serialize(ts::Arc* arc) {
        ts::Arc& a = *arc;
        a << size[0] << size[1] << size[2];
        a << id[0] << id[1] << id[2];

        for(int i = 0; i < size[0] * size[1]; ++i) {
            a << XY1[i] << XY2[i];
        }

        for(int i = 0; i < size[0] * size[2]; ++i) {
            a << XZ1[i] << XZ2[i];
        }

        for(int i = 0; i < size[1] * size[2]; ++i) {
            a << YZ1[i] << YZ2[i];
        }

        for(int i = 0; i < size[0]; ++i) {
            a << XZ1XY1[i] << XZ1XY2[i] << XZ2XY1[i] << XZ2XY2[i];
        }

        for(int i = 0; i < size[1]; ++i) {
            a << XY1YZ1[i] << XY1YZ2[i] << XY2YZ1[i] << XY2YZ2[i];
        }

        for(int i = 0; i < size[2]; ++i) {
            a << XZ1YZ1[i] << XZ1YZ2[i] << XZ2YZ1[i] << XZ2YZ2[i];
        }

        a << XZ1XY1YZ1;
        a << XZ1XY2YZ1;
        a << XZ2XZ1YZ1;
        a << XZ2XY2YZ1;
        a << XZ1XY1YZ2;
        a << XZ1XY2YZ2;
        a << XZ2XZ1YZ2;
        a << XZ2XY2YZ2;
    }

    static FiBoundary* deserialize(ts::Arc* arc) {
        ts::Arc& a = *arc;

        int size[3];
        a >> size[0];
        a >> size[1];
        a >> size[2];

        FiBoundary* r = new FiBoundary(size[0], size[1], size[2]);
        a >> r->id[0];
        a >> r->id[1];
        a >> r->id[2];

        for(int i = 0; i < size[0] * size[1]; ++i) {
            a >> r->XY1[i];
            a >> r->XY2[i];
        }

        for(int i = 0; i < size[0] * size[2]; ++i) {
            a >> r->XZ1[i];
            a >> r->XZ2[i];
        }

        for(int i = 0; i < size[1] * size[2]; ++i) {
            a >> r->YZ1[i];
            a >> r->YZ2[i];
        }

        for(int i = 0; i < size[0]; ++i) {
            a >> r->XZ1XY1[i];
            a >> r->XZ1XY2[i];
            a >> r->XZ2XY1[i];
            a >> r->XZ2XY2[i];
        }

        for(int i = 0; i < size[1]; ++i) {
            a >> r->XY1YZ1[i];
            a >> r->XY1YZ2[i];
            a >> r->XY2YZ1[i];
            a >> r->XY2YZ2[i];
        }

        for(int i = 0; i < size[2]; ++i) {
            a >> r->XZ1YZ1[i];
            a >> r->XZ1YZ2[i];
            a >> r->XZ2YZ1[i];
            a >> r->XZ2YZ2[i];
        }

        a >> r->XZ1XY1YZ1;
        a >> r->XZ1XY2YZ1;
        a >> r->XZ2XZ1YZ1;
        a >> r->XZ2XY2YZ1;
        a >> r->XZ1XY1YZ2;
        a >> r->XZ1XY2YZ2;
        a >> r->XZ2XZ1YZ2;
        a >> r->XZ2XY2YZ2;

        return r;
    }
};

struct RoBoundary {
    double* XY1;
    double* XY2;
    double* XZ1;
    double* XZ2;
    double* YZ1;
    double* YZ2;

    int _size[3];
    int64_t id[3];

    ~RoBoundary() {
        delete[] XY1;
        delete[] XY2;
        delete[] XZ1;
        delete[] XZ2;
        delete[] YZ1;
        delete[] YZ2;
    }

    RoBoundary(int sizeX, int sizeY, int sizeZ) {
        XY1 = new double[2 * sizeX * sizeY];
        XY2 = new double[2 * sizeX * sizeY];
        XZ1 = new double[2 * sizeX * sizeZ];
        XZ2 = new double[2 * sizeX * sizeZ];
        YZ1 = new double[2 * sizeY * sizeZ];
        YZ2 = new double[2 * sizeY * sizeZ];

        size[0] = sizeX;
        size[1] = sizeY;
        size[2] = sizeZ;
    }

    RoBoundary* copy() {
        RoBoundary* result = new RoBoundary(size[0], size[1], size[2]);
        memcpy(result->XY1, XY1, 2 * size[0] * size[1] * sizeof(double));
        memcpy(result->XY2, XY2, 2 * size[0] * size[1] * sizeof(double));
        memcpy(result->XZ1, XZ1, 2 * size[0] * size[2] * sizeof(double));
        memcpy(result->XZ2, XZ2, 2 * size[0] * size[2] * sizeof(double));
        memcpy(result->YZ1, YZ1, 2 * size[1] * size[2] * sizeof(double));
        memcpy(result->YZ2, YZ2, 2 * size[1] * size[2] * sizeof(double));

        return result;
    }

    void serialize(ts::Arc* arc) {
        ts::Arc& a = *arc;
        a << size[0] << size[1] << size[2];
        a << id[0] << id[1] << id[2];

        for(int i = 0; i < 2 * size[0] * size[1]; ++i) {
            a << XY1[i] << XY2[i];
        }

        for(int i = 0; i < 2 * size[0] * size[2]; ++i) {
            a << XZ1[i] << XZ2[i];
        }

        for(int i = 0; i < 2 * size[1] * size[2]; ++i) {
            a << YZ1[i] << YZ2[i];
        }
    }

    static RoBoundary* deserialize(ts::Arc* arc) {
        ts::Arc& a = *arc;

        int size[3];
        a >> size[0];
        a >> size[1];
        a >> size[2];

        RoBoundary* r = new RoBoundary(size[0], size[1], size[2]);
        a >> r->id[0];
        a >> r->id[1];
        a >> r->id[2];

        for(int i = 0; i < 2 * size[0] * size[1]; ++i) {
            a >> r->XY1[i];
            a >> r->XY2[i];
        }

        for(int i = 0; i < 2 * size[0] * size[2]; ++i) {
            a >> r->XZ1[i];
            a >> r->XZ2[i];
        }

        for(int i = 0; i < 2 * size[1] * size[2]; ++i) {
            a >> r->YZ1[i];
            a >> r->YZ2[i];
        }

        return r;
    }
};

#endif // BOUNDARY

