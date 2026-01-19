#include "xml_parser.h"
#include <fstream>
#include <iostream>

bool XMLParser::parse(const std::string& filename, Camera& camera) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open XML file: " << filename << std::endl;
        return false;
    }

    content_ = std::string((std::istreambuf_iterator<char>(file)),
                            std::istreambuf_iterator<char>());

    std::string cameraSection = getSection("camera");
    camera.width = static_cast<int>(getAttribute(cameraSection, "width"));
    camera.height = static_cast<int>(getAttribute(cameraSection, "height"));
    camera.fovx = getAttribute(cameraSection, "fovx");

    std::string eyeSection = getSection("eye");
    camera.eye[0] = getAttribute(eyeSection, "x");
    camera.eye[1] = getAttribute(eyeSection, "y");
    camera.eye[2] = getAttribute(eyeSection, "z");

    std::string lookatSection = getSection("lookat");
    camera.lookat[0] = getAttribute(lookatSection, "x");
    camera.lookat[1] = getAttribute(lookatSection, "y");
    camera.lookat[2] = getAttribute(lookatSection, "z");

    std::string upSection = getSection("up");
    camera.up[0] = getAttribute(upSection, "x");
    camera.up[1] = getAttribute(upSection, "y");
    camera.up[2] = getAttribute(upSection, "z");

    return true;
}

std::string XMLParser::getSection(const std::string& tagName) {
    std::string startTag = "<" + tagName;
    size_t pos = content_.find(startTag);
    if (pos == std::string::npos) return "";

    size_t endPos = content_.find(">", pos);
    if (endPos == std::string::npos) return "";

    return content_.substr(pos, endPos - pos + 1);
}

double XMLParser::getAttribute(const std::string& section, const std::string& attr) {
    std::string search = attr + "=\"";
    size_t pos = section.find(search);
    if (pos == std::string::npos) return 0.0;

    pos += search.length();
    size_t end = section.find("\"", pos);
    if (end == std::string::npos) return 0.0;

    return std::stod(section.substr(pos, end - pos));
}
