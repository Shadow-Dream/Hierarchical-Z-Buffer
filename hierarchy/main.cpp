#include <iostream>
#include <string>
#include <chrono>
#include "obj_parser.h"
#include "xml_parser.h"
#include "zbuffer_renderer.h"
#include "bmp_writer.h"

int main(int argc, char* argv[]) {
    std::string sceneFolder;
    if (argc == 2) {
        sceneFolder = argv[1];
    }
    else {
        sceneFolder = "scenes/living-room";
    }
    
    if (!sceneFolder.empty() && sceneFolder.back() == '/') {
        sceneFolder.pop_back();
    }

    std::string objFile = sceneFolder + "/scene.obj";
    std::string xmlFile = sceneFolder + "/scene.xml";
    std::string bmpFile = sceneFolder + "/scene.bmp";

    OBJParser objParser;
    if (!objParser.parse(objFile)) {
        return 1;
    }

    Camera camera;
    XMLParser xmlParser;
    if (!xmlParser.parse(xmlFile, camera)) {
        return 1;
    }

    ZBufferRenderer renderer;
    renderer.init(camera.width, camera.height,
                  objParser.getVertexCount(), objParser.getTriangleCount());

    // 多次运行以获得稳定的性能数据
    const int NUM_RUNS = 10;
    double totalTime = 0;

    for (int run = 0; run < NUM_RUNS; run++) {
        auto startTime = std::chrono::high_resolution_clock::now();
        renderer.render(objParser.getVertices(), objParser.getVertexCount(),
                        objParser.getTriangles(), objParser.getTriangleCount(),
                        camera);
        auto endTime = std::chrono::high_resolution_clock::now();
        double renderTime = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        totalTime += renderTime;
        std::cout << "Run " << (run + 1) << ": " << renderTime << " ms" << std::endl;
    }

    std::cout << "\nAverage render time: " << (totalTime / NUM_RUNS) << " ms" << std::endl;

    BMPWriter::write(bmpFile, renderer.getZBuffer(),
                     renderer.getWidth(), renderer.getHeight());

    return 0;
}
