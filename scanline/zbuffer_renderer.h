#ifndef ZBUFFER_RENDERER_H
#define ZBUFFER_RENDERER_H

#include "camera.h"
#include "triangle.h"
#include <vector>

// 性能报告函数
void printPerfReport();

struct ScanNode {
    ScanTri data;
    ScanNode *next;
};

struct ScanBucket {
    ScanNode *prevOfHead;  // head = prevOfHead->next
    ScanNode *tail;        // tail == nullptr 表示 bucket 为空
    int nextY;             // 链表中下一个 bucket 的 y 值，-1 表示无效
};

class ZBufferRenderer
{
public:
    ZBufferRenderer();
    ~ZBufferRenderer();

    void init(int width, int height, int vertexCount, int triangleCount);

    void render(const double *vertices, int vertexCount,
                const int *triangles, int triangleCount,
                const Camera &camera);

    const double *getZBuffer() const { return zbuffer_; }
    int getWidth() const { return width_; }
    int getHeight() const { return height_; }

private:
    int width_;
    int height_;
    double aspect_;
    double tanHalfFovX_;
    double tanHalfFovY_;
    double halfW_;
    double halfH_;

    // 预分配缓冲区
    double *zbuffer_;
    double *screenVerts_;
    std::vector<Triangle> candidateTriangles_;
    std::vector<ScanNode> nodeBuffer_;
    std::vector<ScanBucket> buckets_;
    std::vector<double> clipVertBuffer_;
    int clipVertIdx_;

    void transformVertices(const double *vertices, int vertexCount,
                           const double viewMatrix[3][4], double fovx);

    void clipTriangleOne(double *sv0, double *sv1, double *sv2,
                         bool in0, bool in1, bool in2, double fovx);
    void clipTriangleTwo(double *sv0, double *sv1, double *sv2,
                         bool in0, bool in1, bool in2, double fovx);
    void clipEdge(const double *inside, const double *outside,
                  double fovx, double *result);
    void constructTriangle(double *v0, double *v1, double *v2);

    void scanLine();
};

#endif
