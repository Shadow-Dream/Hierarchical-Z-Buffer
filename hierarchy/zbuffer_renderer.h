#ifndef ZBUFFER_RENDERER_H
#define ZBUFFER_RENDERER_H

#include "camera.h"
#include "triangle.h"
#include <vector>

// 混合策略参数
static const int MAX_QUAD_DEPTH = 6;  // 四叉树最大递归深度
static const int FAIL_THRESHOLD = 1000;  // 失败次数触发更新
static const int TRIANGLE_THRESHOLD = 1000;  // 三角形数量触发更新
static const bool SORT_BEFORE_FILL = true;

struct Node {
    int qIdx, bIdx, left, right, top, bottom, depth;
};

struct Range {
    double xMin, xMax, zMin, xZMin, xZMax;
};

class ZBufferRenderer
{
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
    double aspect_;
    double tanHalfFovX_;
    double tanHalfFovY_;
    double halfW_;
    double halfH_;

    // 预分配缓冲区
    double* zbuffer_;
    double* screenVerts_;
    std::vector<Triangle> triangles_;
    std::vector<double> clipVertBuffer_;
    int clipVertIdx_;

    double* treeZ_;
    int treeSize_;
    double* zs, * xStarts, * xEnds, * xZMins, * xZMaxs;

    // 扫描线专属缓存（按y坐标索引，用于混合策略）
    double* scanXStarts_;   // 每行的xStart
    double* scanXEnds_;     // 每行的xEnd
    double* scanZs_;        // 每行的起始z值

    // 懒更新相关成员
    bool* dirtyNodes_;          // 标记需要更新的节点
    int failCount_;             // 扫描线失败计数
    int triangleCount_;         // 已处理三角形计数
    int lastUpdateTriangle_;    // 上次更新时的三角形计数

    void transformVertices(const double* vertices, int vertexCount,
        const double viewMatrix[3][4], double fovx);

    void clipTriangleOne(double* sv0, double* sv1, double* sv2,
        bool in0, bool in1, bool in2, double fovx);
    void clipTriangleTwo(double* sv0, double* sv1, double* sv2,
        bool in0, bool in1, bool in2, double fovx);
    void clipEdge(const double* inside, const double* outside,
        double fovx, double* result);
    void constructTriangle(double* v0, double* v1, double* v2);

    void fillTriangles();

    void fillTriangle(Triangle& tri, Node& node);

    Node findContain(Triangle& tri, int maxDepth = MAX_QUAD_DEPTH);

    // 混合策略：扫描线填充
    void scanlineFill(Triangle& tri, Node& node);

    // 懒更新相关函数
    void markNodeDirty(int qIdx);
    bool shouldUpdateTree();
    double updateNode(int qIdx, int left, int right, int top, int bottom, int depth);
    void updateDirtyNodes();

    Range precomputeTriangle(Triangle& tri, int top, int bottom, int yIdx, int depth);
};

#endif
