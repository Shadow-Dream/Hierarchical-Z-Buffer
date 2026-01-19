#include "triangle.h"
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <limits>

Triangle::Triangle()
    : v{nullptr, nullptr, nullptr}
{
}

Triangle::Triangle(double *v0, double *v1, double *v2, int yStart, int yEnd)
    : v{v0, v1, v2}, yStart(yStart), yEnd(yEnd)
{
}

ScanTri::ScanTri()
{
}

ScanTri::ScanTri(Triangle &tri)
// : v0{tri.v[0]}, v1{tri.v[1]}, v2{tri.v[2]}
{
    double *v0 = tri.v[0];
    double *v1 = tri.v[1];
    double *v2 = tri.v[2];

    int yStart = tri.yStart;
    int yEnd = tri.yEnd;
    yMid = (int)v1[1];

    double xDelta01 = v1[0] - v0[0];
    double xDelta02 = v2[0] - v0[0];
    double xDelta12 = v2[0] - v1[0];

    double yDelta01 = v1[1] - v0[1];
    double yDelta02 = v2[1] - v0[1];
    double yDelta12 = v2[1] - v1[1];

    double zDelta01 = v1[2] - v0[2];
    double zDelta02 = v2[2] - v0[2];
    double zDelta12 = v2[2] - v1[2];

    if (yMid - yStart == 0)
    {
        if (yEnd - yStart == 0)
        {
            needsInit = false;
            isRightLong = xDelta12 > 0;
            if (isRightLong)
            {
                leftX = v0[0];
                leftXD = xDelta01;
                leftZ = v0[2];
                leftZD = zDelta01;

                rightX = v0[0];
                rightXD = xDelta02;
                rightZ = v0[2];
                rightZD = zDelta02;
            }
            else
            {
                leftX = v0[0];
                leftXD = xDelta02;
                leftZ = v0[2];
                leftZD = zDelta02;

                rightX = v0[0];
                rightXD = xDelta01;
                rightZ = v0[2];
                rightZD = zDelta01;
            }
        }
        else
        {
            needsInit = true;
            yMid = -1;
            double invYDelta12 = 1.0 / yDelta12;
            double invYDelta02 = 1.0 / yDelta02;
            double dir12 = xDelta12 * invYDelta12;
            double slope12 = zDelta12 * invYDelta12;
            double dir02 = xDelta02 * invYDelta02;
            double slope02 = zDelta02 * invYDelta02;

            isRightLong = dir02 < dir12;

            if (isRightLong)
            {
                leftX = v1[0];
                leftXD = dir12;
                leftZ = v1[2];
                leftZD = slope12;
                leftYRes = (double)yStart - v1[1];

                rightX = v0[0];
                rightXD = dir02;
                rightZ = v0[2];
                rightZD = slope02;
                rightYRes = (double)yStart - v0[1];
            }
            else
            {
                leftX = v0[0];
                leftXD = dir02;
                leftZ = v0[2];
                leftZD = slope02;
                leftYRes = (double)yStart - v0[1];

                rightX = v1[0];
                rightXD = dir12;
                rightZ = v1[2];
                rightZD = slope12;
                rightYRes = (double)yStart - v1[1];
            }
        }
    }
    else
    {
        needsInit = true;
        double invYDelta01 = 1.0 / yDelta01;
        double invYDelta02 = 1.0 / yDelta02;
        double dir01 = xDelta01 * invYDelta01;
        double dir02 = xDelta02 * invYDelta02;
        double slope01 = zDelta01 * invYDelta01;
        double slope02 = zDelta02 * invYDelta02;
        isRightLong = dir01 < dir02;

        leftYRes = rightYRes = (double)yStart - v0[1];

        if (isRightLong)
        {
            leftX = v0[0];
            leftXD = dir01;
            leftZ = v0[2];
            leftZD = slope01;

            rightX = v0[0];
            rightXD = dir02;
            rightZ = v0[2];
            rightZD = slope02;
        }
        else
        {
            leftX = v0[0];
            leftXD = dir02;
            leftZ = v0[2];
            leftZD = slope02;

            rightX = v0[0];
            rightXD = dir01;
            rightZ = v0[2];
            rightZD = slope01;
        }

        lastYRes = (double)yMid - v1[1];
        double invYDelta12 = 1.0 / yDelta12;
        lastXD = xDelta12 * invYDelta12;
        lastZD = zDelta12 * invYDelta12;
    }
    double Nx = zDelta01 * yDelta02 - yDelta01 * zDelta02;
    double Nz = xDelta01 * yDelta02 - yDelta01 * xDelta02;
    const double EPSILON = 1e-5;
    if (std::abs(Nz) < EPSILON) {
        dZ = 0;
    }else{
        dZ = Nx / Nz;
    }
    // if (isRightLong)
    //     printf("%d %f\n\n\n", yMid, lastXD);

    // if (std::abs(rightXD - leftXD) > 1e-2){
    //     interX = ((rightX - yDelta02*rightXD) - (leftX - yDelta02*leftXD));
    //     interZ = ((rightZ - yDelta02*rightZD) - (leftZ - yDelta02*leftZD));
    //     dZ = interZ / interX;
    // }else{
    //     dZ = 0;
    // }

    // if (std::abs(dZ) > 10){
    //     printf("(%f, %f, %f)\n(%f, %f, %f)\n(%f, %f, %f)\n\n",
    //                    v0[0], v0[1], v0[2],
    //                    v1[0], v1[1], v1[2],
    //                    v2[0], v2[1], v2[2]);
    //     printf("%f %f %f %f %f %f\n\n", dZ, interZ, rightX, leftX, rightXD, leftXD);
    //     // exit(1);
    // }

    // if (rightX - leftX > 0.1){
    //     interX = (rightX - leftX);
    //     interZ = (rightZ - leftZ);
    //     dZ = interZ / interX;
    // }else{

    // }
}

void ScanTri::increment(int y)
{
    if (needsInit)
        init();

    if (y == yMid)
        swap();

    leftX += leftXD;
    rightX += rightXD;
    leftZ += leftZD;
    rightZ += rightZD;

    // y = y + 1;
    // double yRatio = (y - v[0][1]) / (v[2][1] - v[0][1]);
    // yRatio = std::fmin(1, yRatio);
    // double xl = v[0][0] + (v[2][0] - v[0][0]) * yRatio;
    // double zl = v[0][2] + (v[2][2] - v[0][2]) * yRatio;
    // double xr, zr;
    // if (y > v[1][1])
    // {
    //     yRatio = (y - v[1][1]) / (v[2][1] - v[1][1]);
    //     yRatio = std::fmin(1, yRatio);
    //     xr = v[1][0] + (v[2][0] - v[1][0]) * yRatio;
    //     zr = v[1][2] + (v[2][2] - v[1][2]) * yRatio;
    // }
    // else
    // {
    //     yRatio = (y - v[0][1]) / (v[1][1] - v[0][1]);
    //     yRatio = std::fmin(1, yRatio);
    //     xr = v[0][0] + (v[1][0] - v[0][0]) * yRatio;
    //     zr = v[0][2] + (v[1][2] - v[0][2]) * yRatio;
    // }

    // if (xl > xr)
    // {
    //     std::swap(xl, xr);
    //     std::swap(zl, zr);
    // }

    // if (y < v[2][1] &&
    //     ((int)leftX != int(xl) || (int)rightX != int(xr)))
    // {
    //     printf("(%f, %f, %f)\n(%f, %f, %f)\n(%f, %f, %f)\n\n",
    //            v[0][0], v[0][1], v[0][2],
    //            v[1][0], v[1][1], v[1][2],
    //            v[2][0], v[2][1], v[2][2]);
    //     printf("%f %f \n%f %f \n\n", xl, leftX, xr, rightX);
    //     printf("%d \n%f %f %f\n\n", yMid, leftYRes, leftXD, rightXD);
    //     exit(0);
    // }

    // leftX = xl;
    // leftZ = zl;
    // rightX = xr;
    // rightZ = zr;
}

void ScanTri::init()
{
    leftX = leftX + leftXD * leftYRes;
    leftZ = leftZ + leftZD * leftYRes;
    rightX = rightX + rightXD * rightYRes;
    rightZ = rightZ + rightZD * rightYRes;
    needsInit = false;
}

void ScanTri::swap()
{
    if (isRightLong)
    {
        leftX = leftX + (lastXD - leftXD) * lastYRes;
        leftZ = leftZ + (lastZD - leftZD) * lastYRes;
        leftXD = lastXD;
        leftZD = lastZD;
    }
    else
    {
        rightX = rightX + (lastXD - rightXD) * lastYRes;
        rightZ = rightZ + (lastZD - rightZD) * lastYRes;
        rightXD = lastXD;
        rightZD = lastZD;
    }
}
