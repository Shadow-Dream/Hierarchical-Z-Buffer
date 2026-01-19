#ifndef ZBUFFER_RENDERER_H
#define ZBUFFER_RENDERER_H

#include "camera.h"

class ZBufferRenderer {
public:
    ZBufferRenderer();
    ~ZBufferRenderer();

    void init(int width, int height, int vertexCount, int triangleCount);

    void render(const double* vertices, int vertexCount,
                const int* triangles, int triangleCount,
                const Camera& camera);

    const double* getZBuffer() const { return zbuffer_; }
    int getWidth() const { return width_; }
    int getHeight() const { return height_; }

private:
    int width_;
    int height_;
    double* zbuffer_;
    double* screenVerts_;

    void transformVertices(const double* vertices, int vertexCount,
                           const double viewMatrix[3][4],
                           double fovx, double nearPlane);
    void rasterizeTriangle(const double* sv0, const double* sv1, const double* sv2);
    int clipTriangleNearPlane(const double* sv0, const double* sv1, const double* sv2,
                               double nearZ, double fovx,
                               double outVerts[4][6]);
    void clipEdge(const double* inside, const double* outside,
                  double nearZ, double fovx, double* result);
};

#endif
