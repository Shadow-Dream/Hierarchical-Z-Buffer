#include "zbuffer_renderer.h"
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cmath>

static const double NEAR_PLANE = 0.01;

enum { SV_X = 0, SV_Y = 1, SV_Z = 2, SV_VX = 3, SV_VY = 4, SV_VZ = 5 };

ZBufferRenderer::ZBufferRenderer()
    : width_(0), height_(0), zbuffer_(nullptr), screenVerts_(nullptr) {
}

ZBufferRenderer::~ZBufferRenderer() {
    free(zbuffer_);
    free(screenVerts_);
}

void ZBufferRenderer::init(int width, int height, int vertexCount, int triangleCount) {
    width_ = width;
    height_ = height;
    zbuffer_ = (double*)malloc(width * height * sizeof(double));
    screenVerts_ = (double*)malloc(vertexCount * 6 * sizeof(double));
}

void ZBufferRenderer::render(const double* vertices, int vertexCount,
                              const int* triangles, int triangleCount,
                              const Camera& camera) {
    const double fovx = camera.fovx;

    double viewMatrix[3][4];
    camera.getViewMatrix(viewMatrix);

    // 清空 zbuffer
    int zbufSize = width_ * height_;
    for (int i = 0; i < zbufSize; i++) {
        zbuffer_[i] = std::numeric_limits<double>::max();
    }

    transformVertices(vertices, vertexCount, viewMatrix, fovx, NEAR_PLANE);

    for (int i = 0; i < triangleCount; i++) {
        const int* tri = triangles + i * 3;
        int v0 = tri[0], v1 = tri[1], v2 = tri[2];

        if (v0 < 0 || v0 >= vertexCount ||
            v1 < 0 || v1 >= vertexCount ||
            v2 < 0 || v2 >= vertexCount) {
            continue;
        }

        const double* sv0 = screenVerts_ + v0 * 6;
        const double* sv1 = screenVerts_ + v1 * 6;
        const double* sv2 = screenVerts_ + v2 * 6;

        int insideCount = (sv0[SV_VZ] < -NEAR_PLANE ? 1 : 0)
                        + (sv1[SV_VZ] < -NEAR_PLANE ? 1 : 0)
                        + (sv2[SV_VZ] < -NEAR_PLANE ? 1 : 0);

        if (insideCount == 0) {
            continue;
        } else if (insideCount == 3) {
            rasterizeTriangle(sv0, sv1, sv2);
        } else {
            double clipped[4][6];
            int numVerts = clipTriangleNearPlane(sv0, sv1, sv2, NEAR_PLANE, fovx, clipped);
            if (numVerts >= 3) {
                rasterizeTriangle(clipped[0], clipped[1], clipped[2]);
            }
            if (numVerts == 4) {
                rasterizeTriangle(clipped[0], clipped[2], clipped[3]);
            }
        }
    }
}

void ZBufferRenderer::transformVertices(const double* vertices, int vertexCount,
                                         const double viewMatrix[3][4],
                                         double fovx, double nearPlane) {
    const double aspect = (double)width_ / height_;
    const double fovxRad = fovx * M_PI / 180.0;
    const double tanHalfFovX = std::tan(fovxRad * 0.5);
    const double tanHalfFovY = tanHalfFovX / aspect;
    const double halfW = width_ * 0.5;
    const double halfH = height_ * 0.5;

    for (int i = 0; i < vertexCount; i++) {
        const double* vert = vertices + i * 3;
        double* sv = screenVerts_ + i * 6;

        double vx = viewMatrix[0][0] * vert[0] + viewMatrix[0][1] * vert[1] + viewMatrix[0][2] * vert[2] + viewMatrix[0][3];
        double vy = viewMatrix[1][0] * vert[0] + viewMatrix[1][1] * vert[1] + viewMatrix[1][2] * vert[2] + viewMatrix[1][3];
        double vz = viewMatrix[2][0] * vert[0] + viewMatrix[2][1] * vert[1] + viewMatrix[2][2] * vert[2] + viewMatrix[2][3];

        sv[SV_VX] = vx;
        sv[SV_VY] = vy;
        sv[SV_VZ] = vz;

        if (vz < -nearPlane) {
            double invZ = -1.0 / vz;
            sv[SV_X] = (vx * invZ / tanHalfFovX + 1.0) * halfW;
            sv[SV_Y] = (vy * invZ / tanHalfFovY + 1.0) * halfH;
            sv[SV_Z] = -vz;
        } else {
            sv[SV_X] = 0;
            sv[SV_Y] = 0;
            sv[SV_Z] = 0;
        }
    }
}

void ZBufferRenderer::clipEdge(const double* inside, const double* outside,
                                double nearZ, double fovx, double* result) {
    double t = (-nearZ - inside[SV_VZ]) / (outside[SV_VZ] - inside[SV_VZ]);

    result[SV_VX] = inside[SV_VX] + t * (outside[SV_VX] - inside[SV_VX]);
    result[SV_VY] = inside[SV_VY] + t * (outside[SV_VY] - inside[SV_VY]);
    result[SV_VZ] = -nearZ;

    const double aspect = (double)width_ / height_;
    const double fovxRad = fovx * M_PI / 180.0;
    const double tanHalfFovX = std::tan(fovxRad * 0.5);
    const double tanHalfFovY = tanHalfFovX / aspect;

    double invZ = 1.0 / nearZ;
    result[SV_X] = (result[SV_VX] * invZ / tanHalfFovX + 1.0) * width_ * 0.5;
    result[SV_Y] = (result[SV_VY] * invZ / tanHalfFovY + 1.0) * height_ * 0.5;
    result[SV_Z] = nearZ;
}

int ZBufferRenderer::clipTriangleNearPlane(const double* sv0, const double* sv1, const double* sv2,
                                            double nearZ, double fovx,
                                            double outVerts[4][6]) {
    bool in0 = sv0[SV_VZ] < -nearZ;
    bool in1 = sv1[SV_VZ] < -nearZ;
    bool in2 = sv2[SV_VZ] < -nearZ;
    int insideCount = (in0 ? 1 : 0) + (in1 ? 1 : 0) + (in2 ? 1 : 0);

    if (insideCount == 3) {
        for (int j = 0; j < 6; j++) {
            outVerts[0][j] = sv0[j];
            outVerts[1][j] = sv1[j];
            outVerts[2][j] = sv2[j];
        }
        return 3;
    }

    if (insideCount == 0) return 0;

    if (insideCount == 1) {
        const double *inside, *out1, *out2;
        if (in0) { inside = sv0; out1 = sv1; out2 = sv2; }
        else if (in1) { inside = sv1; out1 = sv2; out2 = sv0; }
        else { inside = sv2; out1 = sv0; out2 = sv1; }

        for (int j = 0; j < 6; j++) outVerts[0][j] = inside[j];
        clipEdge(inside, out1, nearZ, fovx, outVerts[1]);
        clipEdge(inside, out2, nearZ, fovx, outVerts[2]);
        return 3;
    }

    const double *outside, *in1ptr, *in2ptr;
    if (!in0) { outside = sv0; in1ptr = sv1; in2ptr = sv2; }
    else if (!in1) { outside = sv1; in1ptr = sv2; in2ptr = sv0; }
    else { outside = sv2; in1ptr = sv0; in2ptr = sv1; }

    for (int j = 0; j < 6; j++) outVerts[0][j] = in1ptr[j];
    clipEdge(in1ptr, outside, nearZ, fovx, outVerts[1]);
    clipEdge(in2ptr, outside, nearZ, fovx, outVerts[2]);
    for (int j = 0; j < 6; j++) outVerts[3][j] = in2ptr[j];
    return 4;
}

void ZBufferRenderer::rasterizeTriangle(const double* v0, const double* v1, const double* v2) {
    const double* sorted[3] = {v0, v1, v2};
    if (sorted[0][SV_Y] > sorted[1][SV_Y]) std::swap(sorted[0], sorted[1]);
    if (sorted[1][SV_Y] > sorted[2][SV_Y]) std::swap(sorted[1], sorted[2]);
    if (sorted[0][SV_Y] > sorted[1][SV_Y]) std::swap(sorted[0], sorted[1]);

    const double* top = sorted[0];
    const double* mid = sorted[1];
    const double* bot = sorted[2];

    int yStart = std::max(0, (int)std::ceil(top[SV_Y]));
    int yEnd = std::min(height_ - 1, (int)std::floor(bot[SV_Y]));
    if (yStart > yEnd) return;

    double invDy02 = (bot[SV_Y] - top[SV_Y]) > 0.0001 ? 1.0 / (bot[SV_Y] - top[SV_Y]) : 0;
    double invDy01 = (mid[SV_Y] - top[SV_Y]) > 0.0001 ? 1.0 / (mid[SV_Y] - top[SV_Y]) : 0;
    double invDy12 = (bot[SV_Y] - mid[SV_Y]) > 0.0001 ? 1.0 / (bot[SV_Y] - mid[SV_Y]) : 0;
    int yMid = (int)std::ceil(mid[SV_Y]);

    for (int y = yStart; y <= yEnd; y++) {
        double t02 = (y - top[SV_Y]) * invDy02;
        double x1 = top[SV_X] + t02 * (bot[SV_X] - top[SV_X]);
        double z1 = top[SV_Z] + t02 * (bot[SV_Z] - top[SV_Z]);

        double x2, z2;
        if (y < yMid) {
            double t01 = (y - top[SV_Y]) * invDy01;
            x2 = top[SV_X] + t01 * (mid[SV_X] - top[SV_X]);
            z2 = top[SV_Z] + t01 * (mid[SV_Z] - top[SV_Z]);
        } else {
            double t12 = (y - mid[SV_Y]) * invDy12;
            x2 = mid[SV_X] + t12 * (bot[SV_X] - mid[SV_X]);
            z2 = mid[SV_Z] + t12 * (bot[SV_Z] - mid[SV_Z]);
        }

        if (x1 > x2) { std::swap(x1, x2); std::swap(z1, z2); }

        int xStart = std::max(0, (int)std::ceil(x1));
        int xEnd = std::min(width_ - 1, (int)std::floor(x2));
        if (xStart > xEnd) continue;

        double spanWidth = x2 - x1;
        double dzx = (spanWidth > 0.0001) ? (z2 - z1) / spanWidth : 0;
        double z = z1 + (xStart - x1) * dzx;

        double* zbufRow = zbuffer_ + y * width_;
        for (int x = xStart; x <= xEnd; x++) {
            if (z > 0 && z < zbufRow[x]) {
                zbufRow[x] = z;
            }
            z += dzx;
        }
    }
}
