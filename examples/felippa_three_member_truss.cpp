#include <iostream>

#include "../src/libfem2d.h"
#include "laclib.h"

using std::cout;
using std::endl;

void run(int argc, char **argv) {

    // GEOMETRY
    //
    //                      fy=1 ↑
    // ---                       2 →
    //  ↑                      ,'| fx=2
    //  |                    ,'  |
    //  |                  ,'    |
    //  |       EA=200√2 ,'      |
    // 10          (2) ,'        | EA=50
    //  |            ,'          | (1)
    //  |          ,'            |
    //  |        ,'              |
    //  |      ,'    EA=100      |
    //  ↓    ,'       (0)        |
    // ---  0--------------------1
    //     | |                  | |
    //      ⇊ uy=-0.5     uy=0.4 ⇈
    //
    //      |←------- 10 -------→|
    //
    // BOUNDARY CONDITIONS
    //
    // node 0: x-fixed with a vertical displacement: uy = -0.5
    // node 1: x-fixed with a vertical displacement: uy = 0.4
    // node 2: fx = 2.0 and fy = 1.0
    //
    // EXPECTED RESULTS
    //
    // kk * uu = ff
    //
    // correct_uu = {0.0, -0.5, 0.0, 0.4, -0.5, 0.2}
    // correct_ff = {-2.0, -2.0, 0.0, 1.0, 2.0, 1.0}
    //
    // Reference
    // Carlos Felippa I-FEM Page 3-12 Chapter 3 The Direct Stiffness Method II

    // options
    auto solid_triangle = false;
    auto plane_stress = false;
    auto thickness = 1.0;
    auto use_expanded_bdb = true;
    auto use_expanded_bdb_full = false;

    // nodes
    auto coordinates = vector<double>{0.0, 0.0, 10.0, 0.0, 10.0, 10.0};

    // elements
    auto connectivity = vector<size_t>{0, 1, 1, 2, 2, 0};

    // parameters
    auto param_young = vector<double>{100.0, 50.0, 200.0};
    auto param_poisson = vector<double>{};
    auto param_cross_area = vector<double>{1.0, 1.0, SQRT_2};

    // boundary conditions
    map<node_dof_pair_t, double> essential_bcs{
        {{0, AlongX}, 0.0},
        {{0, AlongY}, -0.5},
        {{1, AlongY}, 0.4}};
    map<node_dof_pair_t, double> natural_bcs{
        {{2, AlongX}, 2.0},
        {{2, AlongY}, 1.0},
    };

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
    print_vector("uu(felippa_three_member_truss)", fem->uu);

    // Felippa I-FEM page 3-13
    auto correct_uu = vector<double>{0.0, -0.5, 0.0, 0.4, -0.5, 0.2};

    // check
    if (!equal_vectors_tol(fem->uu, correct_uu, 1e-15)) {
        throw "felippa_three_member_truss failed";
    }
}

MAIN_FUNCTION(run)
