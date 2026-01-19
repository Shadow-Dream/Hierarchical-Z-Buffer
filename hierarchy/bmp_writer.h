#ifndef BMP_WRITER_H
#define BMP_WRITER_H

#include <string>

class BMPWriter {
public:
    static bool write(const std::string& filename,
                      const double* zbuffer,
                      int width, int height);
};

#endif
