#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "read_mesh.h"

std::unique_ptr<CoordinatesAndConnectivity> read_mesh(const std::string &filename) {
    FILE *f = fopen(filename.c_str(), "r");
    if (f == NULL) {
        throw "read_mesh: cannot open file";
    }

    const int line_max = 2048;
    char line[line_max];

    if (fgets(line, line_max, f) == NULL) {
        fclose(f);
        throw "read_mesh: cannot read any line in the file";
    }

    auto mesh = std::unique_ptr<CoordinatesAndConnectivity>{
        new CoordinatesAndConnectivity{
            std::vector<double>(),
            std::vector<size_t>(),
        }};

    return mesh;
}