#include "obj_parser.h"
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iostream>

OBJParser::OBJParser()
    : vertices_(nullptr), triangles_(nullptr),
      vertexCount_(0), triangleCount_(0),
      vertexCapacity_(0), triangleCapacity_(0) {
}

OBJParser::~OBJParser() {
    free(vertices_);
    free(triangles_);
}

void OBJParser::addVertex(double x, double y, double z) {
    if (vertexCount_ >= vertexCapacity_) {
        vertexCapacity_ = vertexCapacity_ == 0 ? 1024 : vertexCapacity_ * 2;
        vertices_ = (double*)realloc(vertices_, vertexCapacity_ * 3 * sizeof(double));
    }
    double* v = vertices_ + vertexCount_ * 3;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    vertexCount_++;
}

void OBJParser::addTriangle(int v0, int v1, int v2) {
    if (triangleCount_ >= triangleCapacity_) {
        triangleCapacity_ = triangleCapacity_ == 0 ? 1024 : triangleCapacity_ * 2;
        triangles_ = (int*)realloc(triangles_, triangleCapacity_ * 3 * sizeof(int));
    }
    int* t = triangles_ + triangleCount_ * 3;
    t[0] = v0;
    t[1] = v1;
    t[2] = v2;
    triangleCount_++;
}

int OBJParser::parseVertexIndex(const char* token) {
    int idx = atoi(token);
    return idx - 1;
}

bool OBJParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open OBJ file: " << filename << std::endl;
        return false;
    }

    vertexCount_ = 0;
    triangleCount_ = 0;

    char line[1024];
    while (file.getline(line, sizeof(line))) {
        char* ptr = line;
        while (*ptr == ' ' || *ptr == '\t') ptr++;

        if (ptr[0] == 'v' && ptr[1] == ' ') {
            double x, y, z;
            sscanf(ptr + 2, "%lf %lf %lf", &x, &y, &z);
            addVertex(x, y, z);
        } else if (ptr[0] == 'f' && ptr[1] == ' ') {
            ptr += 2;
            int indices[16];
            int count = 0;

            char* token = strtok(ptr, " \t\r\n");
            while (token && count < 16) {
                indices[count++] = parseVertexIndex(token);
                token = strtok(nullptr, " \t\r\n");
            }

            for (int i = 1; i + 1 < count; i++) {
                addTriangle(indices[0], indices[i], indices[i + 1]);
            }
        }
    }

    return true;
}
