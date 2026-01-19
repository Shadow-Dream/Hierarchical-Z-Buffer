#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#include <string>

class OBJParser {
public:
    OBJParser();
    ~OBJParser();

    bool parse(const std::string& filename);

    // vertices: double[vertexCount * 3], 每个顶点 {x, y, z}
    const double* getVertices() const { return vertices_; }
    // triangles: int[triangleCount * 3], 每个三角形 {v0, v1, v2}
    const int* getTriangles() const { return triangles_; }
    int getVertexCount() const { return vertexCount_; }
    int getTriangleCount() const { return triangleCount_; }

private:
    double* vertices_;      // [x0,y0,z0, x1,y1,z1, ...]
    int* triangles_;        // [v0,v1,v2, v0,v1,v2, ...]
    int vertexCount_;
    int triangleCount_;
    int vertexCapacity_;
    int triangleCapacity_;

    void addVertex(double x, double y, double z);
    void addTriangle(int v0, int v1, int v2);
    int parseVertexIndex(const char* token);
};

#endif
