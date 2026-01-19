#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>

class Camera {
public:
    int width = 800;
    int height = 600;
    double fovx = 60.0;
    double eye[3] = {0, 0, 0};
    double lookat[3] = {0, 0, -1};
    double up[3] = {0, 1, 0};

    // 计算视图矩阵，输出3x4矩阵（只需要前3行）
    void getViewMatrix(double m[3][4]) const {
        // forward = normalize(lookat - eye)
        double fx = lookat[0] - eye[0];
        double fy = lookat[1] - eye[1];
        double fz = lookat[2] - eye[2];
        double flen = std::sqrt(fx*fx + fy*fy + fz*fz);
        fx /= flen; fy /= flen; fz /= flen;

        // right = normalize(forward × up)
        double rx = fy * up[2] - fz * up[1];
        double ry = fz * up[0] - fx * up[2];
        double rz = fx * up[1] - fy * up[0];
        double rlen = std::sqrt(rx*rx + ry*ry + rz*rz);
        rx /= rlen; ry /= rlen; rz /= rlen;

        // newUp = right × forward
        double ux = ry * fz - rz * fy;
        double uy = rz * fx - rx * fz;
        double uz = rx * fy - ry * fx;

        // 构建视图矩阵
        m[0][0] = rx;  m[0][1] = ry;  m[0][2] = rz;
        m[1][0] = ux;  m[1][1] = uy;  m[1][2] = uz;
        m[2][0] = -fx; m[2][1] = -fy; m[2][2] = -fz;
        m[0][3] = -(rx * eye[0] + ry * eye[1] + rz * eye[2]);
        m[1][3] = -(ux * eye[0] + uy * eye[1] + uz * eye[2]);
        m[2][3] = fx * eye[0] + fy * eye[1] + fz * eye[2];
    }
};

#endif
