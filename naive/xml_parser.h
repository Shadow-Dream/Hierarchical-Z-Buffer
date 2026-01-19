#ifndef XML_PARSER_H
#define XML_PARSER_H

#include <string>
#include "camera.h"

class XMLParser {
public:
    XMLParser() = default;

    bool parse(const std::string& filename, Camera& camera);

private:
    std::string content_;

    double getAttribute(const std::string& section, const std::string& attr);
    std::string getSection(const std::string& tagName);
};

#endif // XML_PARSER_H
