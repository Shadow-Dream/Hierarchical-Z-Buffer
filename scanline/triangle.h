#ifndef TRIANGLE_H
#define TRIANGLE_H

class Triangle
{
public:
    Triangle();
    Triangle(double *v0, double *v1, double *v2, int yStart, int yEnd);

    bool operator<(const Triangle &other) const { return yStart < other.yStart; }

    double *v[3];
    int yStart, yEnd;
};

class ScanTri{
public:
    ScanTri();
    ScanTri(Triangle& tri);
    void increment(int y);
    void swap();
    void init();

    int yMid;
    double leftYRes, rightYRes;
    double leftX, leftXD, leftZ, leftZD;
    double rightX, rightXD, rightZ, rightZD;

    bool isRightLong;
    bool needsInit;

    double lastXD, lastZD, lastYRes;

    double dZ;
    // double *v0, *v1, *v2;
    // double interRatio, interZ, interX;
};

#endif
