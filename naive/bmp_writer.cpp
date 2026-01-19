#include "bmp_writer.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cstdint>
#include <cstdlib>

bool BMPWriter::write(const std::string& filename,
                      const double* zbuffer,
                      int width, int height) {
    // Find depth range
    double minZ = std::numeric_limits<double>::max();
    double maxZ = std::numeric_limits<double>::lowest();

    int pixelCount = width * height;
    for (int i = 0; i < pixelCount; i++) {
        if (zbuffer[i] < std::numeric_limits<double>::max() * 0.9) {
            if (zbuffer[i] < minZ) minZ = zbuffer[i];
            if (zbuffer[i] > maxZ) maxZ = zbuffer[i];
        }
    }

    // 限制可视化的最大深度为10米，提高近处物体的对比度
    // 但如果整个场景都在10米之外，则使用原始范围
    const double MAX_VIS_DEPTH = 10.0;
    double visMaxZ = (minZ < MAX_VIS_DEPTH) ? std::min(maxZ, MAX_VIS_DEPTH) : maxZ;

    // BMP headers
    #pragma pack(push, 1)
    struct BMPHeader {
        uint16_t bfType = 0x4D42;
        uint32_t bfSize;
        uint16_t bfReserved1 = 0;
        uint16_t bfReserved2 = 0;
        uint32_t bfOffBits = 54;
    };

    struct DIBHeader {
        uint32_t biSize = 40;
        int32_t biWidth;
        int32_t biHeight;
        uint16_t biPlanes = 1;
        uint16_t biBitCount = 24;
        uint32_t biCompression = 0;
        uint32_t biSizeImage;
        int32_t biXPelsPerMeter = 2835;
        int32_t biYPelsPerMeter = 2835;
        uint32_t biClrUsed = 0;
        uint32_t biClrImportant = 0;
    };
    #pragma pack(pop)

    int rowSize = ((width * 3 + 3) / 4) * 4;
    int imageSize = rowSize * height;

    BMPHeader bmpHeader;
    bmpHeader.bfSize = 54 + imageSize;

    DIBHeader dibHeader;
    dibHeader.biWidth = width;
    dibHeader.biHeight = height;
    dibHeader.biSizeImage = imageSize;

    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create BMP file: " << filename << std::endl;
        return false;
    }

    file.write(reinterpret_cast<char*>(&bmpHeader), sizeof(bmpHeader));
    file.write(reinterpret_cast<char*>(&dibHeader), sizeof(dibHeader));

    // Write pixels
    unsigned char* row = (unsigned char*)malloc(rowSize);
    double range = visMaxZ - minZ;
    if (range < 0.0001) range = 1.0;
    double invRange = 1.0 / range;

    for (int y = 0; y < height; y++) {
        const double* zbufRow = zbuffer + y * width;
        for (int x = 0; x < width; x++) {
            double z = zbufRow[x];
            unsigned char gray;
            if (z >= std::numeric_limits<double>::max() * 0.9) {
                gray = 0;
            } else {
                // 超过可视化最大深度的显示为最暗
                double clampedZ = std::min(z, visMaxZ);
                double normalized = (visMaxZ - clampedZ) * invRange;
                gray = (unsigned char)(normalized * 255);
            }
            row[x * 3 + 0] = gray;
            row[x * 3 + 1] = gray;
            row[x * 3 + 2] = gray;
        }
        // Padding
        for (int x = width * 3; x < rowSize; x++) {
            row[x] = 0;
        }
        file.write(reinterpret_cast<char*>(row), rowSize);
    }

    free(row);
    file.close();
    return true;
}
