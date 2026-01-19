#include "triangle.h"
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <limits>

Triangle::Triangle()
{
}

Triangle::Triangle(double* v0, double* v1, double* v2,
    int yStart, int yEnd,
    int xStart, int xEnd) : yStart(yStart), yEnd(yEnd),
    xStart(xStart), xEnd(xEnd)
{
    minZ = std::fmin(std::fmin(v0[2], v1[2]), v2[2]);
    maxZ = std::fmax(std::fmax(v0[2], v1[2]), v2[2]);

    // 暂存 v0, v1, v2 到空闲成员变量中
    // v0[0], v0[1], v0[2] -> leftX, leftXD, leftZ
    // v1[0], v1[1], v1[2] -> leftZD, rightX, rightXD
    // v2[0], v2[1], v2[2] -> rightZ, rightZD, leftYRes
    leftX = v0[0];
    leftXD = v0[1];
    leftZ = v0[2];
    leftZD = v1[0];
    rightX = v1[1];
    rightXD = v1[2];
    rightZ = v2[0];
    rightZD = v2[1];
    leftYRes = v2[2];
}

void Triangle::init()
{
    // 从成员变量中取出暂存的 v0, v1, v2
    double v0[3] = {leftX, leftXD, leftZ};
    double v1[3] = {leftZD, rightX, rightXD};
    double v2[3] = {rightZ, rightZD, leftYRes};

    yMid = v1[1];

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
    if (std::abs(Nz) < 1e-5) {
        dZ = 0;
    }
    else {
        dZ = Nx / Nz;
    }
}

void Triangle::increment(int y)
{
    if (needsInit)
        reachStart();

    if (y == yMid)
        reachMid();

    leftX += leftXD;
    rightX += rightXD;
    leftZ += leftZD;
    rightZ += rightZD;
}

void Triangle::reachStart()
{
    leftX = leftX + leftXD * leftYRes;
    leftZ = leftZ + leftZD * leftYRes;
    rightX = rightX + rightXD * rightYRes;
    rightZ = rightZ + rightZD * rightYRes;
    needsInit = false;
}

void Triangle::reachMid()
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
