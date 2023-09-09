#include <cmath>
#include <iostream>
#include <map>

#include "../../src/libfem2d.h"
#include "laclib.h"

using namespace std;

void run(int argc, char **argv) {
    // get arguments from command line
    vector<string> defaults{
        "expanded", // {expanded, expanded_full, classical}
        "7",        // number of runs
    };
    auto args = extract_arguments_or_use_defaults(argc, argv, defaults);

    // method to compute Bᵀ ⋅ D ⋅ B
    bool use_expanded_bdb;
    bool use_expanded_bdb_full;
    auto method = string(args[0]);
    auto label = string();
    if (method == "expanded") {
        label = "     expanded: ";
        use_expanded_bdb = true;
        use_expanded_bdb_full = false;
    } else if (method == "expanded_full") {
        label = "expanded_full: ";
        use_expanded_bdb = true;
        use_expanded_bdb_full = true;
    } else if (args[0] == "classical") {
        label = "    classical: ";
        use_expanded_bdb = false;
        use_expanded_bdb_full = false;
    } else {
        throw "method to compute Bᵀ ⋅ D ⋅ B must be one of {expanded, expanded_full, classical}";
    }

    // number of runs
    size_t number_of_runs = std::atoi(args[1].c_str());

    // load the mesh
    // auto pps = string("1800"); // points
    // auto ccs = string("3387"); // cells
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
    for (size_t run = 0; run < number_of_runs; run++) {
        for (size_t e = 0; e < ncell; e++) {
            fem->calculate_element_stiffness(e);
        }
    }
    stopwatch.stop(label, true);
}

MAIN_FUNCTION(run)
