#include <cmath>
#include <iostream>
#include <map>

#include "../src/libfem2d.h"
#include "laclib.h"

using namespace std;

void run(int argc, char **argv) {
    // load the mesh
    auto pps = string("1800"); // points
    auto ccs = string("3387"); // cells
    auto home = string(std::getenv("HOME"));
    auto fn_mesh = home + string("/Downloads/meshes/quarter_ring2d_" + pps + "points_" + ccs + "cells.msh");
    auto mesh = read_mesh(fn_mesh);

    // parameters
    auto ncell = mesh->connectivity.size() / 3;
    auto param_young = vector<double>(ncell, 1000.0);
    auto param_poisson = vector<double>(ncell, 0.25);
    auto param_cross_area = vector<double>{};

    // essential boundary conditions
    auto npoint = mesh->coordinates.size() / 2;
    map<node_dof_pair_t, double> essential_bcs{};
    for (size_t p = 0; p < npoint; p++) {
        auto x = mesh->coordinates[p * 2];
        auto y = mesh->coordinates[p * 2 + 1];
        if (fabs(x) < 1e-11) {
            // fix left edge horizontally
            essential_bcs[{p, AlongX}] = 0.0;
        }
        if (fabs(y) < 1e-11) {
            // fix bottom edge vertically
            essential_bcs[{p, AlongY}] = 0.0;
        }
    }

    // natural boundary conditions
    // TODO

    // options
    // auto solid_triangle = true;
    // auto plane_stress = false;
    // auto thickness = 1.0;
    // auto use_expanded_bdb = true;
    // auto use_expanded_bdb_full = false;
}

MAIN_FUNCTION(run)