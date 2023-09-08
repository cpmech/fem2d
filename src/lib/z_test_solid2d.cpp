#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <map>
#include <vector>

#include "../util/doctest.h"
#include "constants.h"
#include "fem2d.h"
#include "laclib.h"

using namespace std;

#define _SUBCASE(name) if (false)

TEST_CASE("solid2d") {
    SUBCASE("plane-stress bracket (Bhatti example 1.6)") {
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

        // input data
        auto solid_triangle = true;
        auto plane_stress = true;
        auto thickness = 0.25;
        auto coordinates = vector<double>{
            0.0, 0.0,  // 0
            0.0, 2.0,  // 1
            2.0, 0.0,  // 2
            2.0, 1.5,  // 3
            4.0, 0.0,  // 4
            4.0, 1.0}; // 5
        auto connectivity = vector<size_t>{
            0, 2, 3,  // 0
            3, 1, 0,  // 1
            2, 4, 5,  // 2
            5, 3, 2}; // 3
        auto param_young = vector<double>{10000.0, 10000.0, 10000.0, 10000.0};
        auto param_poisson = vector<double>{0.2, 0.2, 0.2, 0.2};
        auto param_cross_area = vector<double>{};
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

        // allocate fem solver
        auto fem = Fem2d::make_new(solid_triangle,
                                   plane_stress,
                                   thickness,
                                   coordinates,
                                   connectivity,
                                   param_young,
                                   param_poisson,
                                   param_cross_area,
                                   essential_bcs,
                                   natural_bcs);

        // check element stiffness
        fem->calculate_element_stiffness(0);
        auto correct_kk0 = Matrix::from_row_major(
            6, 6,
            vector<double>{
                9.765625000000001e+02, 0.000000000000000e+00, -9.765625000000001e+02, 2.604166666666667e+02, 0.000000000000000e+00, -2.604166666666667e+02,  // 0
                0.000000000000000e+00, 3.906250000000000e+02, 5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02, 0.000000000000000e+00,  // 1
                -9.765625000000001e+02, 5.208333333333334e+02, 1.671006944444445e+03, -7.812500000000000e+02, -6.944444444444445e+02, 2.604166666666667e+02, // 2
                2.604166666666667e+02, -3.906250000000000e+02, -7.812500000000000e+02, 2.126736111111111e+03, 5.208333333333334e+02, -1.736111111111111e+03, // 3
                0.000000000000000e+00, -5.208333333333334e+02, -6.944444444444445e+02, 5.208333333333334e+02, 6.944444444444445e+02, 0.000000000000000e+00,  // 4
                -2.604166666666667e+02, 0.000000000000000e+00, 2.604166666666667e+02, -1.736111111111111e+03, 0.000000000000000e+00, 1.736111111111111e+03}  // 5
        );
        CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk0->data, 1e-12));

        fem->calculate_element_stiffness(1);
        auto correct_kk1 = Matrix::from_row_major(
            6, 6,
            vector<double>{
                1.302083333333333e+03, 0.000000000000000e+00, -9.765625000000001e+02, 2.604166666666667e+02, -3.255208333333334e+02, -2.604166666666667e+02, // 0
                0.000000000000000e+00, 5.208333333333334e+02, 5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02, -1.302083333333333e+02, // 1
                -9.765625000000001e+02, 5.208333333333334e+02, 1.253255208333333e+03, -5.859375000000000e+02, -2.766927083333334e+02, 6.510416666666666e+01, // 2
                2.604166666666667e+02, -3.906250000000000e+02, -5.859375000000000e+02, 1.595052083333333e+03, 3.255208333333333e+02, -1.204427083333333e+03, // 3
                -3.255208333333334e+02, -5.208333333333334e+02, -2.766927083333334e+02, 3.255208333333333e+02, 6.022135416666667e+02, 1.953125000000000e+02, // 4
                -2.604166666666667e+02, -1.302083333333333e+02, 6.510416666666666e+01, -1.204427083333333e+03, 1.953125000000000e+02, 1.334635416666667e+03} // 5
        );
        // fem->kk_element->print();
        CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk1->data, 1e-12));

        fem->calculate_element_stiffness(2);
        auto correct_kk2 = Matrix::from_row_major(
            6, 6,
            vector<double>{
                6.510416666666667e+02, 0.000000000000000e+00, -6.510416666666667e+02, 2.604166666666667e+02, 0.000000000000000e+00, -2.604166666666667e+02,  // 0
                0.000000000000000e+00, 2.604166666666667e+02, 5.208333333333334e+02, -2.604166666666667e+02, -5.208333333333334e+02, 0.000000000000000e+00,  // 1
                -6.510416666666667e+02, 5.208333333333334e+02, 1.692708333333333e+03, -7.812500000000000e+02, -1.041666666666667e+03, 2.604166666666667e+02, // 2
                2.604166666666667e+02, -2.604166666666667e+02, -7.812500000000000e+02, 2.864583333333333e+03, 5.208333333333334e+02, -2.604166666666667e+03, // 3
                0.000000000000000e+00, -5.208333333333334e+02, -1.041666666666667e+03, 5.208333333333334e+02, 1.041666666666667e+03, 0.000000000000000e+00,  // 4
                -2.604166666666667e+02, 0.000000000000000e+00, 2.604166666666667e+02, -2.604166666666667e+03, 0.000000000000000e+00, 2.604166666666667e+03}  // 5
        );
        CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk2->data, 1e-12));

        fem->calculate_element_stiffness(3);
        auto correct_kk3 = Matrix::from_row_major(
            6, 6,
            vector<double>{
                9.765625000000001e+02, 0.000000000000000e+00, -6.510416666666667e+02, 2.604166666666667e+02, -3.255208333333334e+02, -2.604166666666667e+02, // 0
                0.000000000000000e+00, 3.906250000000000e+02, 5.208333333333334e+02, -2.604166666666667e+02, -5.208333333333334e+02, -1.302083333333333e+02, // 1
                -6.510416666666667e+02, 5.208333333333334e+02, 1.128472222222222e+03, -5.208333333333334e+02, -4.774305555555556e+02, 0.000000000000000e+00, // 2
                2.604166666666667e+02, -2.604166666666667e+02, -5.208333333333334e+02, 1.909722222222222e+03, 2.604166666666667e+02, -1.649305555555555e+03, // 3
                -3.255208333333334e+02, -5.208333333333334e+02, -4.774305555555556e+02, 2.604166666666667e+02, 8.029513888888889e+02, 2.604166666666667e+02, // 4
                -2.604166666666667e+02, -1.302083333333333e+02, 0.000000000000000e+00, -1.649305555555555e+03, 2.604166666666667e+02, 1.779513888888889e+03} // 5
        );
        CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk3->data, 1e-12));

        fem->calculate_rhs_and_global_stiffness();
        auto kk = fem->kk_coo->as_matrix();
        // kk->print();

        // solve the linear system
        fem->solve();
        print_vector("uu", fem->uu);

        // check solution
        auto correct_uu = vector<double>{
            0.000000000000000e+00, 0.000000000000000e+00,   // 0
            0.000000000000000e+00, 0.000000000000000e+00,   // 1
            -1.035527877607004e-02, -2.552969847657423e-02, // 2
            4.727650463081949e-03, -2.473565538172127e-02,  // 3
            -1.313941349422282e-02, -5.549310752960183e-02, // 4
            8.389015766816341e-05, -5.556637423271112e-02}; // 5
        CHECK(equal_vectors_tol(fem->uu, correct_uu, 1e-15));
    }
}
