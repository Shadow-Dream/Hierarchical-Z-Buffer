#ifndef TRIANGLE_H
#define TRIANGLE_H

class Triangle
{
public:
    Triangle();
    Triangle(double* v0, double* v1, double* v2,
        int yStart, int yEnd, int xStart, int xEnd);

    void init();
    void increment(int y);
    void reachMid();
    void reachStart();

    int yStart, yMid, yEnd;
    int xStart, xEnd;

    int& minX() noexcept { return xStart; }
    int& maxX() noexcept { return xEnd; }
    int& minY() noexcept { return yStart; }
    int& maxY() noexcept { return yEnd; }

    double leftYRes, rightYRes;
    double leftX, leftXD, leftZ, leftZD;
    double rightX, rightXD, rightZ, rightZD;

    bool isRightLong;
    bool needsInit;

    double lastXD, lastZD, lastYRes;
    double dZ = 0, minZ, maxZ;
};

#endif
