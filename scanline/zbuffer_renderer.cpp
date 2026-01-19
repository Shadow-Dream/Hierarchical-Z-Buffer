#include "zbuffer_renderer.h"
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cstdio>

// 性能测量变量
static double t_clearZbuf = 0;
static double t_transform = 0;
static double t_clip = 0;
static double t_sort = 0;
static double t_scanLine = 0;
static double t_clearBuckets = 0;
static double t_buildActiveList = 0;
static double t_pixelFill = 0;
static double t_removeNodes = 0;
static double t_triIncrement = 0;
static int perfRunCount = 0;

void printPerfReport()
{
    if (perfRunCount == 0)
        return;
    printf("\n=== Performance Report (avg of %d runs) ===\n", perfRunCount);
    printf("render() breakdown:\n");
    printf("  1. clearZbuffer:      %.4f ms\n", t_clearZbuf / perfRunCount);
    printf("  2. transformVertices: %.4f ms\n", t_transform / perfRunCount);
    printf("  3. clipTriangles:     %.4f ms\n", t_clip / perfRunCount);
    printf("  4. sort:              %.4f ms\n", t_sort / perfRunCount);
    printf("  5. scanLine total:    %.4f ms\n", t_scanLine / perfRunCount);
    printf("\nscanLine() breakdown:\n");
    printf("  5a. clearBuckets:     %.4f ms\n", t_clearBuckets / perfRunCount);
    printf("  5b. buildActiveList:  %.4f ms\n", t_buildActiveList / perfRunCount);
    printf("  5c. pixelFill:        %.4f ms\n", t_pixelFill / perfRunCount);
    printf("  5d. triIncrement:     %.4f ms\n", t_triIncrement / perfRunCount);
    printf("  5e. removeNodes:      %.4f ms\n", t_removeNodes / perfRunCount);
    printf("==========================================\n");
}

static const double NEAR_PLANE = 0.01;

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
    : width_(0), height_(0), zbuffer_(nullptr), screenVerts_(nullptr), clipVertIdx_(0)
{
}

ZBufferRenderer::~ZBufferRenderer()
{
    free(zbuffer_);
    free(screenVerts_);
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
    candidateTriangles_.reserve(triangleCount * 2);

    // 节点缓冲区
    nodeBuffer_.resize(triangleCount * 2);

    // 桶数组
    buckets_.resize(height);
}

void ZBufferRenderer::render(const double *vertices, int vertexCount,
                             const int *triangles, int triangleCount,
                             const Camera &camera)
{

    const double fovx = camera.fovx;

    aspect_ = (double)width_ / height_;
    const double fovxRad = fovx * M_PI / 180.0;
    tanHalfFovX_ = std::tan(fovxRad * 0.5);
    tanHalfFovY_ = tanHalfFovX_ / aspect_;
    halfW_ = width_ * 0.5;
    halfH_ = height_ * 0.5;

    double viewMatrix[3][4];
    camera.getViewMatrix(viewMatrix);

    std::fill(zbuffer_, zbuffer_ + width_ * height_, std::numeric_limits<double>::max());

    // 重置裁剪缓冲区索引
    clipVertIdx_ = 0;

    // 清空候选三角形
    candidateTriangles_.clear();

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

    
    std::sort(candidateTriangles_.begin(), candidateTriangles_.end());
    

    scanLine();


    perfRunCount++;
}

void ZBufferRenderer::scanLine()
{

    if (candidateTriangles_.empty())
        return;

    int bufIdx = 0;

    // 5a. 清空桶
    std::memset(buckets_.data(), 0, buckets_.size() * sizeof(ScanBucket));
    

    ScanNode dummy;
    dummy.next = nullptr;
    int dummyY = -1;  // dummy->next 所属的 bucket 的 y 值，-1 表示链表为空

    int y = 0;
    size_t i = 0;

    double localBuildTime = 0;
    double localPixelTime = 0;
    double localIncrTime = 0;
    double localRemoveTime = 0;

    while (y < height_)
    {
        // 5b. 构建活动列表
        
        for (; i < candidateTriangles_.size() && candidateTriangles_[i].yStart <= y; i++)
        {
            int yEnd = std::min(candidateTriangles_[i].yEnd, height_ - 1);

            ScanNode *node = &nodeBuffer_[bufIdx++];
            node->data = ScanTri(candidateTriangles_[i]);

            ScanBucket &bkt = buckets_[yEnd];
            if (bkt.tail == nullptr)
            {
                node->next = dummy.next;
                dummy.next = node;
                bkt.tail = node;
                bkt.prevOfHead = &dummy;
                bkt.nextY = dummyY;
                // 更新原来紧跟dummy的bucket的prevOfHead
                if (dummyY >= 0)
                {
                    buckets_[dummyY].prevOfHead = node;
                }
                dummyY = yEnd;
            }
            else
            {
                node->next = bkt.prevOfHead->next;
                bkt.prevOfHead->next = node;
            }
        }
        
        
        for (ScanNode *node = dummy.next; node; node = node->next)
        {
            ScanTri &tri = node->data;

            int xStart = std::max(0, int(tri.leftX));
            int xEnd = std::min(width_, int(tri.rightX));
            // printf("%d: %d %d\n", y, xStart, xEnd); 
            double dz = tri.dZ;
            // double refDz = (tri.rightZ - tri.leftZ) / (tri.rightX - tri.leftX);
            // printf("%d %f %f\n", y, dz, refDz); 
            // if (std::abs(dzRef - dz) > 0.01)
            // {
            //     printf("(%f, %f, %f)\n(%f, %f, %f)\n(%f, %f, %f)\n\n(%f, %f, %f, %f, %f, %f, %f)\n",
            //            tri.v0[0], tri.v0[1], tri.v0[2],
            //            tri.v1[0], tri.v1[1], tri.v1[2],
            //            tri.v2[0], tri.v2[1], tri.v2[2],
            //            tri.interRatio, tri.interZ, tri.interX, dz, dzRef, tri.rightZ - tri.leftZ, tri.rightX - tri.leftX);
            //     exit(0);
            // }
            double z = tri.leftZ + dz * std::fmax(xStart - tri.leftX,0);
            double *zbufRow = zbuffer_ + y * width_;
            for (int x = xStart; x < xEnd; z += dz, x++)
            {
                
                if (z < zbufRow[x])
                {
                    zbufRow[x] = z;
                    // printf("%d %d: %f\n\n", x, y, z);
                }
            }
            node->data.increment(y);
        }
        
        // 5d. 三角形增量更新
        
        for (ScanNode *node = dummy.next; node; node = node->next)
        {
            
        }
        
        
        ScanBucket &bkt = buckets_[y];
        if (bkt.tail != nullptr)
        {
            bkt.prevOfHead->next = bkt.tail->next;
            // 通过 nextY 链找到下一个有效的 bucket（跳过已删除的）
            int nextY = bkt.nextY;
            while (nextY >= 0 && buckets_[nextY].tail == nullptr)
            {
                nextY = buckets_[nextY].nextY;
            }
            if (nextY >= 0)
            {
                buckets_[nextY].prevOfHead = bkt.prevOfHead;
            }
            // 如果删除的是链表第一个 bucket，更新 dummyY
            if (y == dummyY)
            {
                dummyY = nextY;
            }
            bkt.tail = nullptr;
        }
        
        if (dummy.next == nullptr)
            y = i < candidateTriangles_.size() ? candidateTriangles_[i].yStart : height_;
        else
            y++;
    }

    t_buildActiveList += localBuildTime;
    t_pixelFill += localPixelTime;
    t_triIncrement += localIncrTime;
    t_removeNodes += localRemoveTime;
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

    int yStart = (int)v0[1];
    if (yStart >= height_)
        return;
    int yEnd = (int)v2[1];
    if (yEnd < 0)
        return;
    int xStart = (int)std::min(v0[0], std::min(v1[0], v2[0]));
    if (xStart >= width_)
        return;
    int xEnd = (int)std::max(v0[0], std::max(v1[0], v2[0]));
    if (xEnd < 0)
        return;

    if (v0[1] < 0)
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

            candidateTriangles_.emplace_back(newVerts, newVerts + 3, v2, 0, yEnd);
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

            candidateTriangles_.emplace_back(newVerts, newVerts + 3, v2, 0, yEnd);
            candidateTriangles_.emplace_back(newVerts, v1, v2, 0, yEnd);
        }
    }
    else
    {
        candidateTriangles_.emplace_back(v0, v1, v2, yStart, yEnd);
    }
}
