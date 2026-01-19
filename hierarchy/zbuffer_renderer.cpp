#include "zbuffer_renderer.h"
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <immintrin.h>

static const double NEAR_PLANE = 0.01;

int quadIndexToBinIndex(int i)
{
    int x = 3 * i + 1;
    int L = ((31 - __builtin_clz(x)) >> 1);
    int pow4 = 1 << (L << 1);
    int startQ = (pow4 - 1) / 3;
    int startB = (1 << L) - 1;

    int p = i - startQ;
    int q = _pext_u32(p, 0xAAAAAAAA);
    return startB + q;
}

enum
{
    SV_X = 0,
    SV_Y = 1,
    SV_Z = 2,
    SV_VX = 3,
    SV_VY = 4,
    SV_VZ = 5
};

ZBufferRenderer::ZBufferRenderer()
    : width_(0), height_(0), zbuffer_(nullptr), screenVerts_(nullptr),
      clipVertIdx_(0), scanXStarts_(nullptr), scanXEnds_(nullptr), scanZs_(nullptr),
      dirtyNodes_(nullptr), failCount_(0), triangleCount_(0), lastUpdateTriangle_(0)
{
}

ZBufferRenderer::~ZBufferRenderer()
{
    free(zbuffer_);
    free(screenVerts_);
    free(scanXStarts_);
    free(scanXEnds_);
    free(scanZs_);
    free(dirtyNodes_);
}

void ZBufferRenderer::transformVertices(const double *vertices, int vertexCount,
                                        const double viewMatrix[3][4], double fovx)
{
    for (int i = 0; i < vertexCount; i++)
    {
        const double *vert = vertices + i * 3;
        double *sv = screenVerts_ + i * 6;

        double vx = viewMatrix[0][0] * vert[0] + viewMatrix[0][1] * vert[1] + viewMatrix[0][2] * vert[2] + viewMatrix[0][3];
        double vy = viewMatrix[1][0] * vert[0] + viewMatrix[1][1] * vert[1] + viewMatrix[1][2] * vert[2] + viewMatrix[1][3];
        double vz = viewMatrix[2][0] * vert[0] + viewMatrix[2][1] * vert[1] + viewMatrix[2][2] * vert[2] + viewMatrix[2][3];

        sv[SV_VX] = vx;
        sv[SV_VY] = vy;
        sv[SV_VZ] = vz;

        if (vz < -NEAR_PLANE)
        {
            double invZ = -1.0 / vz;
            sv[SV_X] = (vx * invZ / tanHalfFovX_ + 1.0) * halfW_;
            sv[SV_Y] = (vy * invZ / tanHalfFovY_ + 1.0) * halfH_;
            sv[SV_Z] = -vz;
        }
        else
        {
            sv[SV_X] = 0;
            sv[SV_Y] = 0;
            sv[SV_Z] = 0;
        }
    }
}

void ZBufferRenderer::clipEdge(const double *inside, const double *outside,
                               double fovx, double *result)
{
    double t = (-NEAR_PLANE - inside[SV_VZ]) / (outside[SV_VZ] - inside[SV_VZ]);

    double svx = inside[SV_VX] + t * (outside[SV_VX] - inside[SV_VX]);
    double svy = inside[SV_VY] + t * (outside[SV_VY] - inside[SV_VY]);

    double invZ = 1.0 / NEAR_PLANE;
    result[SV_X] = (svx * invZ / tanHalfFovX_ + 1.0) * halfW_;
    result[SV_Y] = (svy * invZ / tanHalfFovY_ + 1.0) * halfH_;
    result[SV_Z] = NEAR_PLANE;
}

void ZBufferRenderer::clipTriangleOne(double *sv0, double *sv1, double *sv2,
                                      bool in0, bool in1, bool in2, double fovx)
{
    double *inside, *out1, *out2;
    if (in0)
    {
        inside = sv0;
        out1 = sv1;
        out2 = sv2;
    }
    else if (in1)
    {
        inside = sv1;
        out1 = sv2;
        out2 = sv0;
    }
    else
    {
        inside = sv2;
        out1 = sv0;
        out2 = sv1;
    }

    double *newVerts = &clipVertBuffer_[clipVertIdx_];
    clipVertIdx_ += 6;
    clipEdge(inside, out1, fovx, newVerts + 0);
    clipEdge(inside, out2, fovx, newVerts + 3);
    constructTriangle(inside, newVerts + 0, newVerts + 3);
}

void ZBufferRenderer::clipTriangleTwo(double *sv0, double *sv1, double *sv2,
                                      bool in0, bool in1, bool in2, double fovx)
{
    double *outside, *in1ptr, *in2ptr;
    if (!in0)
    {
        outside = sv0;
        in1ptr = sv1;
        in2ptr = sv2;
    }
    else if (!in1)
    {
        outside = sv1;
        in1ptr = sv2;
        in2ptr = sv0;
    }
    else
    {
        outside = sv2;
        in1ptr = sv0;
        in2ptr = sv1;
    }

    double *newVerts = &clipVertBuffer_[clipVertIdx_];
    clipVertIdx_ += 6;
    clipEdge(in1ptr, outside, fovx, newVerts + 0);
    clipEdge(in2ptr, outside, fovx, newVerts + 3);
    constructTriangle(in1ptr, newVerts + 0, newVerts + 3);
    constructTriangle(in1ptr, newVerts + 3, in2ptr);
}

void ZBufferRenderer::constructTriangle(double *v0, double *v1, double *v2)
{
    if (v0[1] < v1[1])
    {
        if (v1[1] > v2[1])
        {
            if (v0[1] < v2[1])
                std::swap(v1, v2);
            else
            {
                double *tmp = v0;
                v0 = v2;
                v2 = v1;
                v1 = tmp;
            }
        }
    }
    else
    {
        if (v0[1] < v2[1])
            std::swap(v0, v1);
        else
        {
            if (v1[1] < v2[1])
            {
                double *tmp = v0;
                v0 = v1;
                v1 = v2;
                v2 = tmp;
            }
            else
                std::swap(v0, v2);
        }
    }

    int yStart = v0[1];
    if (yStart >= height_)
        return;
    int yEnd = v2[1];
    if (yEnd < 0)
        return;

    int xStart, xEnd;

    if (yStart < 0)
    {
        double *newVerts = &clipVertBuffer_[clipVertIdx_];
        clipVertIdx_ += 6;

        if (v1[1] < 0)
        {
            double t02 = -v0[1] / (v2[1] - v0[1]);
            newVerts[0] = v0[0] + (v2[0] - v0[0]) * t02;
            newVerts[1] = 0.0;
            newVerts[2] = v0[2] + (v2[2] - v0[2]) * t02;

            double t12 = -v1[1] / (v2[1] - v1[1]);
            newVerts[3] = v1[0] + (v2[0] - v1[0]) * t12;
            newVerts[4] = 0.0;
            newVerts[5] = v1[2] + (v2[2] - v1[2]) * t12;

            xStart = std::min(newVerts[0], std::min(newVerts[3], v2[0]));
            if (xStart >= width_)
                return;
            xEnd = std::max(newVerts[0], std::max(newVerts[3], v2[0]));
            if (xEnd < 0)
                return;
            triangles_.emplace_back(newVerts, newVerts + 3, v2, 0, yEnd, xStart, xEnd);
        }
        else
        {
            double t01 = -v0[1] / (v1[1] - v0[1]);
            newVerts[0] = v0[0] + (v1[0] - v0[0]) * t01;
            newVerts[1] = 0.0;
            newVerts[2] = v0[2] + (v1[2] - v0[2]) * t01;

            double t02 = -v0[1] / (v2[1] - v0[1]);
            newVerts[3] = v0[0] + (v2[0] - v0[0]) * t02;
            newVerts[4] = 0.0;
            newVerts[5] = v0[2] + (v2[2] - v0[2]) * t02;

            xStart = std::min(newVerts[0], std::min(newVerts[3], v2[0]));
            if (xStart < width_)
            {
                xEnd = std::max(newVerts[0], std::max(newVerts[3], v2[0]));
                if (xEnd >= 0)
                    triangles_.emplace_back(
                        newVerts, newVerts + 3, v2, 0, yEnd, xStart, xEnd);
            }

            xStart = std::min(newVerts[0], std::min(v1[0], v2[0]));
            if (xStart >= width_)
                return;
            xEnd = std::max(newVerts[0], std::max(v1[0], v2[0]));
            if (xEnd < 0)
                return;
            triangles_.emplace_back(newVerts, v1, v2, 0, yEnd, xStart, xEnd);
        }
    }
    else
    {
        xStart = std::min(v0[0], std::min(v1[0], v2[0]));
        if (xStart >= width_)
            return;
        xEnd = std::max(v0[0], std::max(v1[0], v2[0]));
        if (xEnd < 0)
            return;
        triangles_.emplace_back(v0, v1, v2, yStart, yEnd, xStart, xEnd);
    }
}

void ZBufferRenderer::init(int width, int height, int vertexCount, int triangleCount)
{
    width_ = width;
    height_ = height;

    // zbuffer
    zbuffer_ = (double *)malloc(width * height * sizeof(double));

    // screenVerts
    screenVerts_ = (double *)malloc(vertexCount * 6 * sizeof(double));

    // 裁剪顶点缓冲区 (每个三角形最多产生2个新顶点，每个6个double，再加上constructTriangle中的裁剪)
    clipVertBuffer_.resize(triangleCount * 4 * 6);

    // 候选三角形 (最坏情况：每个三角形裁剪成2个)
    triangles_.reserve(triangleCount * 2);

    int maxAxis = std::max(width, height);

    int n, m;

    n = 33 - (int)__builtin_clz(maxAxis - 1);
    m = 1 << n;
    m = m - 1;
    n = n << 1;
    n = 1 << n;
    n = (n - 1) / 3;
    // printf("%d %d", m, n);
    // exit(0);

    treeZ_ = (double *)malloc(n * sizeof(double));
    zs = (double *)malloc(n * sizeof(double));
    xStarts = (double *)malloc(m * sizeof(double));
    xEnds = (double *)malloc(m * sizeof(double));
    xZMins = (double *)malloc(m * sizeof(double));
    xZMaxs = (double *)malloc(m * sizeof(double));
    dirtyNodes_ = (bool *)malloc(n * sizeof(bool));

    // 扫描线专属缓存
    scanXStarts_ = (double *)malloc(height * sizeof(double));
    scanXEnds_ = (double *)malloc(height * sizeof(double));
    scanZs_ = (double *)malloc(height * sizeof(double));

    treeSize_ = n;
}

void ZBufferRenderer::render(const double *vertices, int vertexCount,
                             const int *triangles, int triangleCount,
                             const Camera &camera)
{
    const double fovx = camera.fovx;

    aspect_ = (double)width_ / height_;
    const double fovxRad = fovx * 3.14159265358979323846 / 180.0;
    tanHalfFovX_ = std::tan(fovxRad * 0.5);
    tanHalfFovY_ = tanHalfFovX_ / aspect_;
    halfW_ = width_ * 0.5;
    halfH_ = height_ * 0.5;

    double viewMatrix[3][4];
    camera.getViewMatrix(viewMatrix);

    std::fill(zbuffer_, zbuffer_ + width_ * height_, std::numeric_limits<double>::max());
    std::fill(treeZ_, treeZ_ + treeSize_, std::numeric_limits<double>::max());
    clipVertIdx_ = 0;
    triangles_.clear();
    std::fill(dirtyNodes_, dirtyNodes_ + treeSize_, false);
    failCount_ = 0;
    triangleCount_ = 0;
    lastUpdateTriangle_ = 0;

    transformVertices(vertices, vertexCount, viewMatrix, fovx);

    // 3. 裁剪三角形
    for (int i = 0; i < triangleCount; i++)
    {
        const int *tri = triangles + i * 3;
        int v0 = tri[0], v1 = tri[1], v2 = tri[2];

        double *sv0 = screenVerts_ + v0 * 6;
        double *sv1 = screenVerts_ + v1 * 6;
        double *sv2 = screenVerts_ + v2 * 6;

        bool in0 = sv0[SV_VZ] < -NEAR_PLANE;
        bool in1 = sv1[SV_VZ] < -NEAR_PLANE;
        bool in2 = sv2[SV_VZ] < -NEAR_PLANE;

        switch (in0 + in1 + in2)
        {
        case 1:
            clipTriangleOne(sv0, sv1, sv2, in0, in1, in2, fovx);
            break;
        case 2:
            clipTriangleTwo(sv0, sv1, sv2, in0, in1, in2, fovx);
            break;
        case 3:
            constructTriangle(sv0, sv1, sv2);
            break;
        }
    }

    fillTriangles();
}

Node ZBufferRenderer::findContain(Triangle &tri, int maxDepth)
{
    int xStart = std::max(tri.xStart, 0);
    int yStart = std::max(tri.yStart, 0);
    int xEnd = std::min(tri.xEnd, width_ - 1);
    int yEnd = std::min(tri.yEnd + 1, height_ - 1);
    int left = 0, right = width_, top = 0, bottom = height_;

    int xRange = (width_ >> 1), yRange = (height_ >> 1);

    int nodeIdx = 0;
    int depth = 0;
    while (depth < maxDepth)
    {
        double z = treeZ_[nodeIdx];
        if (tri.minZ > z)
        {
            return {-1, 0, 0, 0, 0, 0, 0};
        }
        int childIds = (nodeIdx << 2) + 1;
        bool childContain = false;

        for (int i = 0; i < 4; i++)
        {
            int childLeft, childRight, childTop, childBottom;

            if (i & 1)
            {
                childLeft = left + xRange;
                childRight = right;
            }
            else
            {
                childLeft = left;
                childRight = left + xRange;
            }

            if (i & 2)
            {
                childTop = top + yRange;
                childBottom = bottom;
            }
            else
            {
                childTop = top;
                childBottom = top + yRange;
            }

            bool inChild = (xStart >= childLeft && xEnd <= childRight && yStart >= childTop && yEnd <= childBottom);

            if (inChild)
            {
                childContain = true;
                left = childLeft;
                right = childRight;
                top = childTop;
                bottom = childBottom;
                nodeIdx = childIds + i;
                break;
            }
        }

        if (!childContain)
        {
            break;
        }

        xRange = right - left;
        yRange = bottom - top;
        xRange = xRange >> 1;
        yRange = yRange >> 1;
        depth++; // 进入下一层
    }

    int bIdx = quadIndexToBinIndex(nodeIdx);

    return {nodeIdx, bIdx, left, right, top, bottom, depth};
}

void ZBufferRenderer::fillTriangle(Triangle &tri, Node &node)
{
    double xStart = xStarts[node.bIdx];
    double xEnd = xEnds[node.bIdx];
    double maxZ = treeZ_[node.qIdx];
    int xStartInt = std::max(int(xStart), 0);
    int xEndInt = std::min(int(xEnd), width_);

    if (node.left >= xEndInt || node.right <= xStartInt)
        return;

    double minZ;
    if (tri.dZ > 0)
    {
        minZ = zs[node.bIdx] + tri.dZ * std::fmax(node.left - xZMins[node.bIdx], 0);
    }
    else
    {
        minZ = zs[node.bIdx] + tri.dZ * (std::fmin(node.right, xZMaxs[node.bIdx]) - xZMins[node.bIdx]);
    }

    if (minZ > maxZ)
        return;

    if (node.depth == MAX_QUAD_DEPTH)
    {
        markNodeDirty(node.qIdx);
        scanlineFill(tri, node);
        return;
    }

    int xRange = (node.right - node.left) >> 1;
    int yRange = (node.bottom - node.top) >> 1;
    int childQIds = (node.qIdx << 2) + 1;
    int childBIds = (node.bIdx << 1) + 1;
    for (int i = 0; i < 4; i++)
    {
        Node childNode;
        childNode.depth = node.depth + 1;
        childNode.qIdx = childQIds + i;

        if (i & 2)
        {
            childNode.bIdx = childBIds + 1;
            childNode.top = node.top + yRange;
            childNode.bottom = node.bottom;
        }
        else
        {
            childNode.bIdx = childBIds;
            childNode.top = node.top;
            childNode.bottom = node.top + yRange;
        }

        if (childNode.top > tri.yEnd || childNode.bottom <= tri.yStart)
        {
            continue;
        }

        if (i & 1)
        {
            childNode.left = node.left + xRange;
            childNode.right = node.right;
        }
        else
        {
            childNode.left = node.left;
            childNode.right = node.left + xRange;
        }
        fillTriangle(tri, childNode);
    }
}

void ZBufferRenderer::scanlineFill(Triangle &tri, Node &node)
{
    int yStart = std::max(node.top, tri.yStart);
    int yEnd = std::min(node.bottom, tri.yEnd + 1);

    for (int y = yStart; y < yEnd; y++)
    {
        double xStart = scanXStarts_[y];
        double xEnd = scanXEnds_[y];
        double z = scanZs_[y];

        int xStartInt = std::max((int)xStart, node.left);
        int xEndInt = std::min((int)xEnd, node.right);

        double xOffset = xStartInt - xStart;
        xOffset = xOffset > 0 ? xOffset : 0;
        double curZ = z + tri.dZ * xOffset;

        double *row = zbuffer_ + y * width_;

        for (int x = xStartInt; x < xEndInt; x++, curZ += tri.dZ)
        {
            if (curZ < row[x])
            {
                row[x] = curZ;
            }
            else
            {
                failCount_++;
            }
        }
    }
}

Range ZBufferRenderer::precomputeTriangle(Triangle &tri, int top, int bottom, int yIdx, int depth)
{
    if (bottom <= tri.yStart || top > tri.yEnd)
    {
        return {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::min(),
            std::numeric_limits<double>::max(),
            0, 0};
    }
    int yRange = (bottom - top) >> 1;

    // 达到深度限制或叶节点，填充扫描线缓存
    if (depth == MAX_QUAD_DEPTH || !yRange)
    {
        double xMin = std::numeric_limits<double>::max();
        double xMax = std::numeric_limits<double>::min();
        double zMin = std::numeric_limits<double>::max();
        double xZMin = 0, xZMax = 0;

        int yStart = std::max(tri.yStart, top);
        int yEnd = std::min(tri.yEnd + 1, bottom);

        for (int y = yStart; y < yEnd; y++)
        {
            double xStart = tri.leftX;
            double xEnd = tri.rightX;
            double z = tri.leftZ;

            scanXStarts_[y] = xStart;
            scanXEnds_[y] = xEnd;
            scanZs_[y] = z;

            xMin = std::min(xMin, xStart);
            xMax = std::max(xMax, xEnd);
            if (z < zMin)
            {
                zMin = z;
                xZMin = xStart;
                xZMax = xEnd;
            }
            tri.increment(y);
        }

        xStarts[yIdx] = xMin;
        xEnds[yIdx] = xMax;
        zs[yIdx] = zMin;
        xZMins[yIdx] = xZMin;
        xZMaxs[yIdx] = xZMax;

        return {xMin, xMax, zMin, xZMin, xZMax};
    }

    Range range1 = precomputeTriangle(tri, top, top + yRange, (yIdx << 1) + 1, depth + 1);
    Range range2 = precomputeTriangle(tri, top + yRange, bottom, (yIdx << 1) + 2, depth + 1);
    double xStart = std::min(range1.xMin, range2.xMin);
    double xEnd = std::max(range1.xMax, range2.xMax);
    double zMin, xZMin, xZMax;
    if (range1.zMin < range2.zMin)
    {
        zMin = range1.zMin;
        xZMin = range1.xZMin;
        xZMax = range1.xZMax;
    }
    else
    {
        zMin = range2.zMin;
        xZMin = range2.xZMin;
        xZMax = range2.xZMax;
    }
    xStarts[yIdx] = xStart;
    xEnds[yIdx] = xEnd;
    zs[yIdx] = zMin;
    xZMins[yIdx] = xZMin;
    xZMaxs[yIdx] = xZMax;
    return {xStart, xEnd, zMin, xZMin, xZMax};
}

// 标记节点为脏（向上传播到根）
void ZBufferRenderer::markNodeDirty(int qIdx)
{
    while (qIdx >= 0)
    {
        if (dirtyNodes_[qIdx])
            break; // 已标记则停止（祖先也已标记）
        dirtyNodes_[qIdx] = true;
        qIdx = (qIdx - 1) >> 2; // 父节点索引
    }
}

// 判断是否需要更新四叉树
bool ZBufferRenderer::shouldUpdateTree()
{
    // 策略1：失败次数阈值
    if (failCount_ >= FAIL_THRESHOLD)
    {
        return true;
    }

    // 策略2：处理三角形数量阈值
    if (triangleCount_ - lastUpdateTriangle_ >= TRIANGLE_THRESHOLD)
    {
        return true;
    }

    return false;
}

// 自顶向下递归更新脏节点
double ZBufferRenderer::updateNode(int qIdx, int left, int right, int top, int bottom, int depth)
{
    if (!dirtyNodes_[qIdx])
    {
        return treeZ_[qIdx];
    }

    double maxZ = 0;
    if (depth == MAX_QUAD_DEPTH)
    {
        double curZ = treeZ_[qIdx];
        for (int y = top; y < bottom; y++)
        {
            double *row = zbuffer_ + y * width_;
            for (int x = left; x < right; x++)
            {
                double z = row[x];
                maxZ = maxZ > z ? maxZ : z;
                if (maxZ >= curZ)
                {
                    dirtyNodes_[qIdx] = false;
                    return curZ;
                }
            }
        }
        treeZ_[qIdx] = maxZ;
        dirtyNodes_[qIdx] = false;
        return maxZ;
    }

    int xRange = (right - left) >> 1;
    int yRange = (bottom - top) >> 1;
    int childBase = (qIdx << 2) + 1;

    for (int i = 0; i < 4; i++)
    {
        int childLeft, childRight, childTop, childBottom;

        if (i & 1)
        {
            childLeft = left + xRange;
            childRight = right;
        }
        else
        {
            childLeft = left;
            childRight = left + xRange;
        }

        if (i & 2)
        {
            childTop = top + yRange;
            childBottom = bottom;
        }
        else
        {
            childTop = top;
            childBottom = top + yRange;
        }
        double childZ = updateNode(childBase + i, childLeft, childRight, childTop, childBottom, depth + 1);
        maxZ = maxZ > childZ ? maxZ : childZ;
    }

    treeZ_[qIdx] = maxZ;
    dirtyNodes_[qIdx] = false;
    return maxZ;
}

// 执行懒更新
void ZBufferRenderer::updateDirtyNodes()
{
    if (dirtyNodes_[0])
    {

        updateNode(0, 0, width_, 0, height_, 0);
    }

    failCount_ = 0;
    lastUpdateTriangle_ = triangleCount_;
}

void ZBufferRenderer::fillTriangles()
{
    if (SORT_BEFORE_FILL)
    {
        std::sort(triangles_.begin(), triangles_.end(),
                  [](const Triangle &a, const Triangle &b)
                  {
                      return a.minZ < b.minZ;
                  });
    }

    for (auto &tri : triangles_)
    {
        // 检查是否需要更新四叉树
        if (shouldUpdateTree())
        {
            updateDirtyNodes();
        }

        Node node = findContain(tri, MAX_QUAD_DEPTH);
        if (node.qIdx == -1)
        {
            continue;
        }
        tri.init();
        precomputeTriangle(tri, node.top, node.bottom, node.bIdx, node.depth);
        fillTriangle(tri, node);

        triangleCount_++;
    }
}
