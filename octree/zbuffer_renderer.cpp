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
static const double FAR_PLANE = 11;

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
    : width_(0), height_(0), zbuffer_(nullptr), screenVerts_(nullptr), clipVertIdx_(0),
      scanXStarts_(nullptr), scanXEnds_(nullptr), scanZs_(nullptr),
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

void ZBufferRenderer::transformVertices(const double* vertices, int vertexCount,
    const double viewMatrix[3][4], double fovx)
{
    for (int i = 0; i < vertexCount; i++)
    {
        const double* vert = vertices + i * 3;
        double* sv = screenVerts_ + i * 6;

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

void ZBufferRenderer::clipEdge(const double* inside, const double* outside,
    double fovx, double* result)
{
    double t = (-NEAR_PLANE - inside[SV_VZ]) / (outside[SV_VZ] - inside[SV_VZ]);

    double svx = inside[SV_VX] + t * (outside[SV_VX] - inside[SV_VX]);
    double svy = inside[SV_VY] + t * (outside[SV_VY] - inside[SV_VY]);

    double invZ = 1.0 / NEAR_PLANE;
    result[SV_X] = (svx * invZ / tanHalfFovX_ + 1.0) * halfW_;
    result[SV_Y] = (svy * invZ / tanHalfFovY_ + 1.0) * halfH_;
    result[SV_Z] = NEAR_PLANE;
}

void ZBufferRenderer::clipTriangleOne(double* sv0, double* sv1, double* sv2,
    bool in0, bool in1, bool in2, double fovx)
{
    double* inside, * out1, * out2;
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

    double* newVerts = &clipVertBuffer_[clipVertIdx_];
    clipVertIdx_ += 6;
    clipEdge(inside, out1, fovx, newVerts + 0);
    clipEdge(inside, out2, fovx, newVerts + 3);
    constructTriangle(inside, newVerts + 0, newVerts + 3);
}

void ZBufferRenderer::clipTriangleTwo(double* sv0, double* sv1, double* sv2,
    bool in0, bool in1, bool in2, double fovx)
{
    double* outside, * in1ptr, * in2ptr;
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

    double* newVerts = &clipVertBuffer_[clipVertIdx_];
    clipVertIdx_ += 6;
    clipEdge(in1ptr, outside, fovx, newVerts + 0);
    clipEdge(in2ptr, outside, fovx, newVerts + 3);
    constructTriangle(in1ptr, newVerts + 0, newVerts + 3);
    constructTriangle(in1ptr, newVerts + 3, in2ptr);
}

void ZBufferRenderer::constructTriangle(double* v0, double* v1, double* v2)
{
    if (v0[1] < v1[1])
    {
        if (v1[1] > v2[1])
        {
            if (v0[1] < v2[1])
                std::swap(v1, v2);
            else
            {
                double* tmp = v0;
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
                double* tmp = v0;
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
        double* newVerts = &clipVertBuffer_[clipVertIdx_];
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
    zbuffer_ = (double*)malloc(width * height * sizeof(double));

    // screenVerts
    screenVerts_ = (double*)malloc(vertexCount * 6 * sizeof(double));

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

    treeZ_ = (double*)malloc(n * sizeof(double));
    zs = (double*)malloc(n * sizeof(double));
    xStarts = (double*)malloc(m * sizeof(double));
    xEnds = (double*)malloc(m * sizeof(double));
    xZMins = (double*)malloc(m * sizeof(double));
    xZMaxs = (double*)malloc(m * sizeof(double));
    treeSize_ = n;

    // 扫描线专属缓存
    scanXStarts_ = (double*)malloc(height * sizeof(double));
    scanXEnds_ = (double*)malloc(height * sizeof(double));
    scanZs_ = (double*)malloc(height * sizeof(double));

    // 懒更新
    dirtyNodes_ = (bool*)malloc(n * sizeof(bool));

    // 单数组存储三角形索引和临时缓冲区
    triIndices_ = (int*)malloc(triangleCount * 2 * sizeof(int));
    triIndices2_ = (int*)malloc(triangleCount * 2 * sizeof(int));
    triChildIdx_ = (int8_t*)malloc(triangleCount * 2 * sizeof(int8_t));
    int maxOctreeNode = std::min(triangleCount * 2, width * height * 2);
    octreeNodes = (OctreeNode*)malloc(maxOctreeNode * sizeof(OctreeNode));
}

void ZBufferRenderer::render(const double* vertices, int vertexCount,
    const int* triangles, int triangleCount,
    const Camera& camera)
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
    std::fill(dirtyNodes_, dirtyNodes_ + treeSize_, false);
    clipVertIdx_ = 0;
    triangles_.clear();
    failCount_ = 0;
    triangleCount_ = 0;
    lastUpdateTriangle_ = 0;
    octreeSize = 0;

    transformVertices(vertices, vertexCount, viewMatrix, fovx);

    // 3. 裁剪三角形
    for (int i = 0; i < triangleCount; i++)
    {
        const int* tri = triangles + i * 3;
        int v0 = tri[0], v1 = tri[1], v2 = tri[2];

        double* sv0 = screenVerts_ + v0 * 6;
        double* sv1 = screenVerts_ + v1 * 6;
        double* sv2 = screenVerts_ + v2 * 6;

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

    buildOctree();
    Node root = { 0, 0, 0, width_, 0, height_, 0 };  // depth=0
    fillOctree(0, root, 0, width_, 0, height_, NEAR_PLANE, FAR_PLANE, 1);
}

Node ZBufferRenderer::findContain(Triangle& tri)
{
    int xStart = std::max(tri.xStart, 0);
    int yStart = std::max(tri.yStart, 0);
    int xEnd = std::min(tri.xEnd, width_ - 1);
    int yEnd = std::min(tri.yEnd + 1, height_ - 1);
    int left = 0, right = width_, top = 0, bottom = height_;

    int xRange = (width_ >> 1), yRange = (height_ >> 1);

    int nodeIdx = 0;
    while (xRange | yRange)
    {
        double z = treeZ_[nodeIdx];
        if (tri.minZ > z)
        {
            return { -1, 0, 0, 0, 0 };
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
    }
    int bIdx = quadIndexToBinIndex(nodeIdx);

    return { nodeIdx, bIdx, left, right, top, bottom, 0 };  // 添加depth=0
}

void ZBufferRenderer::fillTriangle(Triangle& tri, Node& node)
{
    if (node.right == node.left || node.top == node.bottom)
    {
        return;
    }

    double xStart = xStarts[node.bIdx];
    double xEnd = xEnds[node.bIdx];
    double maxZ = treeZ_[node.qIdx];
    int xStartInt = std::max(int(xStart), 0);
    int xEndInt = std::min(int(xEnd), width_);

    if (node.left >= xEndInt || node.right <= xStartInt)
    {
        return;
    }

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
    {
        return;
    }

    // 深度限制：切换到扫描线
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
        childNode.qIdx = childQIds + i;
        childNode.depth = node.depth + 1;

        if (i & 2)
        {
            if (yRange)
            {
                childNode.bIdx = childBIds + 1;
            }
            else
            {
                childNode.bIdx = node.bIdx;
            }

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

void ZBufferRenderer::buildOctree()
{
    // 初始化triIndices_
    int n = (int)triangles_.size();
    for (int i = 0; i < n; i++)
    {
        triIndices_[i] = i;
    }

    // 只创建根节点，标记为未分区（children[0] = -1）
    OctreeNode& root = octreeNodes[octreeSize++];
    root.triStart = 0;
    root.triEnd = n;
    root.minZ = NEAR_PLANE;
    for (int i = 0; i < 8; i++)
        root.children[i] = -1;  // -1 表示未初始化
}

// 延迟初始化：分区单个节点（不递归）
void ZBufferRenderer::partitionNode(int nodeIdx,
    int minX, int maxX, int minY, int maxY,
    double minZ, double maxZ, int xyDepth)
{
    OctreeNode& node = octreeNodes[nodeIdx];
    int start = node.triStart;
    int end = node.triEnd;
    int triCount = end - start;

    int xRange = (maxX - minX) >> 1;
    int yRange = (maxY - minY) >> 1;
    int middleX = minX + xRange;
    int middleY = minY + yRange;
    double middleZ = (minZ + maxZ) / 2;

    // 叶节点：不再分区
    if (triCount <= octreeThreshold || xyDepth == MAX_QUAD_DEPTH)
    {
        for (int i = 0; i < 8; i++)
            node.children[i] = 0;
        return;
    }

    // 第一遍：计数并保存归属
    int counts[9] = {0};
    for (int j = start; j < end; j++)
    {
        Triangle& tri = triangles_[triIndices_[j]];
        int c;
        if (tri.maxX() <= middleX && tri.maxY() < middleY && tri.maxZ <= middleZ)
            c = 0;
        else if (tri.minX() >= middleX && tri.maxY() < middleY && tri.maxZ <= middleZ)
            c = 1;
        else if (tri.maxX() <= middleX && tri.minY() >= middleY && tri.maxZ <= middleZ)
            c = 2;
        else if (tri.minX() >= middleX && tri.minY() >= middleY && tri.maxZ <= middleZ)
            c = 3;
        else if (tri.maxX() <= middleX && tri.maxY() < middleY && tri.minZ >= middleZ)
            c = 4;
        else if (tri.minX() >= middleX && tri.maxY() < middleY && tri.minZ >= middleZ)
            c = 5;
        else if (tri.maxX() <= middleX && tri.minY() >= middleY && tri.minZ >= middleZ)
            c = 6;
        else if (tri.minX() >= middleX && tri.minY() >= middleY && tri.minZ >= middleZ)
            c = 7;
        else
            c = 8;
        triChildIdx_[j] = c;
        counts[c]++;
    }

    // 计算每个子节点的起始位置
    int offsets[10];
    offsets[0] = start;
    for (int i = 0; i < 9; i++)
        offsets[i + 1] = offsets[i] + counts[i];

    // 第二遍：直接复制到正确位置
    int cursors[9];
    for (int i = 0; i < 9; i++)
        cursors[i] = offsets[i];

    for (int j = start; j < end; j++)
    {
        int c = triChildIdx_[j];
        triIndices2_[cursors[c]++] = triIndices_[j];
    }

    // 复制回原数组
    memcpy(triIndices_ + start, triIndices2_ + start, triCount * sizeof(int));

    // 当前节点的三角形在最后
    node.triStart = offsets[8];
    node.triEnd = end;

    // 创建子节点（但不分区，标记为 -1 表示未初始化）
    for (int c = 0; c < 8; c++)
    {
        if (counts[c] > 0)
        {
            node.children[c] = octreeSize;
            OctreeNode& child = octreeNodes[octreeSize++];
            child.triStart = offsets[c];
            child.triEnd = offsets[c + 1];

            // 计算子节点的 minZ
            if (c & 4)
                child.minZ = middleZ;
            else
                child.minZ = minZ;

            // 标记为未分区
            for (int i = 0; i < 8; i++)
                child.children[i] = -1;
        }
        else
        {
            node.children[c] = 0;
        }
    }
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
            double* row = zbuffer_ + y * width_;
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

// 扫描线填充
void ZBufferRenderer::scanlineFill(Triangle& tri, Node& node)
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

        double* row = zbuffer_ + y * width_;

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

Range ZBufferRenderer::precomputeTriangle(Triangle& tri, int top, int bottom, int yIdx, int depth)
{
    if (bottom <= tri.yStart || top > tri.yEnd)
    {
        return {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::min(),
            std::numeric_limits<double>::max(),
            0, 0 };
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

        return { xMin, xMax, zMin, xZMin, xZMax };
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
    return { xStart, xEnd, zMin, xZMin, xZMax };
}

void ZBufferRenderer::fillOctree(int octreeIdx, Node& node,
    int minX, int maxX, int minY, int maxY,
    double minZ, double maxZ, int xyDepth)
{
    OctreeNode& octree = octreeNodes[octreeIdx];
    if (octree.minZ > treeZ_[node.qIdx])
    {
        return;
    }

    // 延迟初始化：如果节点未分区，先分区
    if (octree.children[0] == -1)
    {
        partitionNode(octreeIdx, minX, maxX, minY, maxY, minZ, maxZ, xyDepth);
    }

    for (int i = octree.triStart; i < octree.triEnd; i++)
    {
        // 检查是否需要更新四叉树
        if (shouldUpdateTree())
        {
            updateDirtyNodes();
        }

        Triangle& tri = triangles_[triIndices_[i]];
        tri.init();
        precomputeTriangle(tri, node.top, node.bottom, node.bIdx, node.depth);
        fillTriangle(tri, node);

        triangleCount_++;  // 计数
    }

    int childQIds = (node.qIdx << 2) + 1;
    int childBIds = (node.bIdx << 1) + 1;
    int xRange = (node.right - node.left) >> 1;
    int yRange = (node.bottom - node.top) >> 1;

    int xRangeSpace = (maxX - minX) >> 1;
    int yRangeSpace = (maxY - minY) >> 1;
    int middleX = minX + xRangeSpace;
    int middleY = minY + yRangeSpace;
    double middleZ = (minZ + maxZ) / 2;

    for (int i = 0; i < 8; i++)
    {
        int childIdx = octree.children[i];
        if (childIdx <= 0)  // 0 表示空，-1 不应该出现（已分区）
            continue;

        Node childNode;
        if (xRange | yRange)
        {
            childNode.qIdx = childQIds + (i & 3);
            childNode.depth = node.depth + 1;
        }
        else
        {
            childNode.qIdx = node.qIdx;
            childNode.depth = node.depth;
        }

        if (i & 2)
        {
            if (yRange)
            {
                childNode.bIdx = childBIds + 1;
            }
            else
            {
                childNode.bIdx = node.bIdx;
            }

            childNode.top = node.top + yRange;
            childNode.bottom = node.bottom;
        }
        else
        {
            if (!yRange)
            {
                continue;
            }
            childNode.bIdx = childBIds;
            childNode.top = node.top;
            childNode.bottom = node.top + yRange;
        }

        if (i & 1)
        {
            childNode.left = node.left + xRange;
            childNode.right = node.right;
        }
        else
        {
            if (!xRange)
            {
                continue;
            }
            childNode.left = node.left;
            childNode.right = node.left + xRange;
        }

        // 计算子节点的空间边界
        int cMinX, cMaxX, cMinY, cMaxY;
        double cMinZ, cMaxZ;

        if (i & 1) { cMinX = middleX; cMaxX = maxX; }
        else { cMinX = minX; cMaxX = middleX; }

        if (i & 2) { cMinY = middleY; cMaxY = maxY; }
        else { cMinY = minY; cMaxY = middleY; }

        if (i & 4) { cMinZ = middleZ; cMaxZ = maxZ; }
        else { cMinZ = minZ; cMaxZ = middleZ; }

        fillOctree(childIdx, childNode, cMinX, cMaxX, cMinY, cMaxY, cMinZ, cMaxZ, xyDepth + 1);
    }
}