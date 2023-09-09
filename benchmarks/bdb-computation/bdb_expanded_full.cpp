#include <cmath>
#include <iostream>
#include <map>

#include "../../src/libfem2d.h"
#include "laclib.h"

using namespace std;

void run(int argc, char **argv) {
    // load the mesh
    auto pps = string("1648167"); // points
    auto ccs = string("3291387"); // cells
    auto home = string(std::getenv("HOME"));
    auto fn_mesh = home + string("/Downloads/meshes/quarter_ring2d_" + pps + "points_" + ccs + "cells.msh");
    auto mesh = read_mesh(fn_mesh);

    // parameters
    auto ncell = mesh->connectivity.size() / 3;
    auto param_young = vector<double>(ncell, 1000.0);
    auto param_poisson = vector<double>(ncell, 0.25);
    auto param_cross_area = vector<double>{};

    // boundary conditions
    map<node_dof_pair_t, double> essential_bcs{};
    map<node_dof_pair_t, double> natural_bcs{};

    // options
    auto solid_triangle = true;
    auto plane_stress = false;
    auto thickness = 1.0;
    auto use_expanded_bdb = true;
    auto use_expanded_bdb_full = true;

    // allocate fem
    auto fem = Fem2d::make_new(solid_triangle,
                               plane_stress,
                               thickness,
                               use_expanded_bdb,
                               use_expanded_bdb_full,
                               mesh->coordinates,
                               mesh->connectivity,
                               param_young,
                               param_poisson,
                               param_cross_area,
                               essential_bcs,
                               natural_bcs);

    // calculate kk for all elements
    auto stopwatch = Stopwatch::make_new();
    for (size_t e = 0; e < ncell; e++) {
        fem->calculate_element_stiffness(e);
    }
    stopwatch.stop("all elements: ", true);
}

MAIN_FUNCTION(run)