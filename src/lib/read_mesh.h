#pragma once

#include <memory>
#include <string>
#include <vector>

struct CoordinatesAndConnectivity {
    std::vector<double> coordinates;
    std::vector<size_t> connectivity;
};

std::unique_ptr<CoordinatesAndConnectivity> read_mesh(const std::string &filename);
