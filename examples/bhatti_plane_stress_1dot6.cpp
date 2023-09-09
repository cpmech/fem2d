#include <iostream>

#include "../src/libfem2d.h"
#include "laclib.h"

using namespace std;

void run(int argc, char **argv) {

    // Bhatti's Example 1.6 on page 32
    //
    // Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    //
    // TEST GOAL
    //
    // This test verifies the equilibrium of a thin bracket modelled by assuming plane-stress
    //
    // MESH
    //
    // 2.0  fixed 1'-,_load                connectivity:
    //            |     '-,_      load      eid : vertices
    // 1.5 - - -  |        ,'3-,__            0 :  0, 2, 3
    //            |  1   ,'  |    '-,_        1 :  3, 1, 0
    // 1.0 - - -  |    ,'    |  3   ,-'5      2 :  2, 4, 5
    //            |  ,'  0   |   ,-'   |      3 :  5, 3, 2
    //            |,'        |,-'   2  |
    // 0.0  fixed 0----------2---------4   constraints:
    //           0.0        2.0       4.0   fixed on x and y
    //
    // BOUNDARY CONDITIONS
    //
    // Fully fixed @ points 0 and 1
    // Distributed load along edges (1,3) and (3,5) with Qn = -20
    //
    // CONFIGURATION AND PARAMETERS
    //
    // Static simulation
    // Young = 10,000
    // Poisson = 0.2
    // Plane-stress with thickness = 0.25

    // options
    auto solid_triangle = true;
    auto plane_stress = true;
    auto thickness = 0.25;
    auto use_expanded_bdb = true;
    auto use_expanded_bdb_full = false;

    // nodes
    auto coordinates = vector<double>{
        0.0, 0.0,  // 0
        0.0, 2.0,  // 1
        2.0, 0.0,  // 2
        2.0, 1.5,  // 3
        4.0, 0.0,  // 4
        4.0, 1.0}; // 5

    // bars
    auto connectivity = vector<size_t>{
        0, 2, 3,  // 0
        3, 1, 0,  // 1
        2, 4, 5,  // 2
        5, 3, 2}; // 3

    // parameters
    auto param_young = vector<double>{10000.0, 10000.0, 10000.0, 10000.0};
    auto param_poisson = vector<double>{0.2, 0.2, 0.2, 0.2};
    auto param_cross_area = vector<double>{};

    // boundary conditions
    map<node_dof_pair_t, double> essential_bcs{
        {{0, AlongX}, 0.0},
        {{0, AlongY}, 0.0},
        {{1, AlongX}, 0.0},
        {{1, AlongY}, 0.0}};
    map<node_dof_pair_t, double> natural_bcs{
        {{1, AlongX}, -1.25},
        {{1, AlongY}, -5.0},
        {{3, AlongX}, -2.5},
        {{3, AlongY}, -10.0},
        {{5, AlongX}, -1.25},
        {{5, AlongY}, -5.0}};

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
    print_vector("uu(bhatti_plane_stress_1dot6)", fem->uu);

    // bhatti's solution
    auto correct_uu = vector<double>{
        0.000000000000000e+00, 0.000000000000000e+00,   // 0
        0.000000000000000e+00, 0.000000000000000e+00,   // 1
        -1.035527877607004e-02, -2.552969847657423e-02, // 2
        4.727650463081949e-03, -2.473565538172127e-02,  // 3
        -1.313941349422282e-02, -5.549310752960183e-02, // 4
        8.389015766816341e-05, -5.556637423271112e-02}; // 5

    // check
    if (!equal_vectors_tol(fem->uu, correct_uu, 1e-15)) {
        throw "bhatti_plane_stress_1dot6 failed";
    }
}

MAIN_FUNCTION(run)
