#include <iostream>

#include "../src/libfem2d.h"
#include "laclib.h"

using std::cout;
using std::endl;

void run(int argc, char **argv) {

    // Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
    // Element Method, Wiley, Fifth Edition, 664p
    //
    // TEST GOAL
    //
    // This test verifies a plane-strain simulation with Tri3 elements
    //
    // MESH
    //
    //               1.0 kN/m²
    //         ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    //  0.0  ▷0---------1---------2
    //         |       ,'|       ,'|   E = 1e6 kN/m²
    //         |  0  ,'  |  2  ,'  |   ν = 0.3
    //         |   ,'    |   ,'    |
    //         | ,'   1  | ,'  3   |   connectivity:
    // -0.5  ▷3'--------4'--------5     0 : 1 0 3
    //         |       ,'|       ,'|     1 : 3 4 1
    //         |  4  ,'  |  6  ,'  |     2 : 2 1 4
    //         |   ,'    |   ,'    |     3 : 4 5 2
    //         | ,'   5  | ,'   7  |     4 : 4 3 6
    // -1.0  ▷6'--------7'--------8     5 : 6 7 4
    //         △        △        △    6 : 5 4 7
    //                                   7 : 7 8 5
    //        0.0       0.5       1.0
    //
    // BOUNDARY CONDITIONS
    //
    // Fix left edge horizontally
    // Fix bottom edge vertically
    // Distributed load Qn = -1.0 on top edge
    //
    // CONFIGURATION AND PARAMETERS
    //
    // Static simulation
    // Young = 1e6
    // Poisson = 0.3
    // Plane-strain

    // options
    auto solid_triangle = true;
    auto plane_stress = false;
    auto thickness = 1.0;
    auto use_expanded_bdb = true;
    auto use_expanded_bdb_full = false;

    // nodes
    auto coordinates = vector<double>{
        0.0, 0.0,   // 0
        0.5, 0.0,   // 1
        1.0, 0.0,   // 2
        0.0, -0.5,  // 3
        0.5, -0.5,  // 4
        1.0, -0.5,  // 5
        0.0, -1.0,  // 6
        0.5, -1.0,  // 7
        1.0, -1.0}; // 8

    // elements
    auto connectivity = vector<size_t>{
        1, 0, 3,  // 0
        3, 4, 1,  // 1
        2, 1, 4,  // 2
        4, 5, 2,  // 3
        4, 3, 6,  // 4
        6, 7, 4,  // 5
        5, 4, 7,  // 6
        7, 8, 5}; // 7

    // parameters
    auto param_young = vector<double>(8, 1e6);
    auto param_poisson = vector<double>(8, 0.3);
    auto param_cross_area = vector<double>{};

    // boundary conditions
    map<node_dof_pair_t, double> essential_bcs{
        {{0, AlongX}, 0.0},
        {{3, AlongX}, 0.0},
        {{6, AlongX}, 0.0},
        {{6, AlongY}, 0.0},
        {{7, AlongY}, 0.0},
        {{8, AlongY}, 0.0}};
    map<node_dof_pair_t, double> natural_bcs{
        {{0, AlongY}, -0.25},
        {{1, AlongY}, -0.50},
        {{2, AlongY}, -0.25}};

    // get FEM solution
    auto fem = Fem2d::make_new(solid_triangle,
                               plane_stress,
                               thickness,
                               use_expanded_bdb,
                               use_expanded_bdb_full,
                               coordinates,
                               connectivity,
                               param_young,
                               param_poisson,
                               param_cross_area,
                               essential_bcs,
                               natural_bcs);
    fem->solve();
    print_vector("uu(smith_plane_strain_5dot2)", fem->uu);

    // check solution
    auto correct_uu = vector<double>{
        0.000000000000000e+00, -9.100000000000005e-07, // 0
        1.950000000000001e-07, -9.100000000000006e-07, // 1
        3.900000000000002e-07, -9.100000000000000e-07, // 2
        0.000000000000000e+00, -4.550000000000002e-07, // 3
        1.950000000000002e-07, -4.550000000000004e-07, // 4
        3.900000000000004e-07, -4.549999999999999e-07, // 5
        0.000000000000000e+00, 0.000000000000000e+00,  // 6
        1.950000000000004e-07, 0.000000000000000e+00,  // 7
        3.900000000000004e-07, 0.000000000000000e+00}; // 8

    // check
    if (!equal_vectors_tol(fem->uu, correct_uu, 1e-15)) {
        throw "smith_plane_strain_5dot2 failed";
    }
}

MAIN_FUNCTION(run)
