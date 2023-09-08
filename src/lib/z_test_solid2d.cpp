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

        // correct values
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
        auto correct_uu = vector<double>{
            0.000000000000000e+00, 0.000000000000000e+00,   // 0
            0.000000000000000e+00, 0.000000000000000e+00,   // 1
            -1.035527877607004e-02, -2.552969847657423e-02, // 2
            4.727650463081949e-03, -2.473565538172127e-02,  // 3
            -1.313941349422282e-02, -5.549310752960183e-02, // 4
            8.389015766816341e-05, -5.556637423271112e-02}; // 5

        SUBCASE("classical Bᵀ ⋅ D ⋅ B approach") {
            // allocate fem solver
            auto use_expanded_bdb = false;
            auto use_expanded_bdb_full = false;
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

            // check element stiffness
            fem->calculate_element_stiffness(0);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk0->data, 1e-12));

            fem->calculate_element_stiffness(1);
            // fem->kk_element->print();
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk1->data, 1e-12));

            fem->calculate_element_stiffness(2);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk2->data, 1e-12));

            fem->calculate_element_stiffness(3);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk3->data, 1e-12));

            // fem->calculate_rhs_and_global_stiffness();
            // auto kk = fem->kk_coo->as_matrix();
            // kk->print();

            // solve the linear system
            fem->solve();
            // print_vector("uu", fem->uu);

            // check solution
            CHECK(equal_vectors_tol(fem->uu, correct_uu, 1e-15));
        }

        SUBCASE("expanded (full) Bᵀ ⋅ D ⋅ B") {
            // allocate fem solver
            auto use_expanded_bdb = false;
            auto use_expanded_bdb_full = true;
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

            // check element stiffness
            fem->calculate_element_stiffness(0);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk0->data, 1e-12));
            // fem->kk_element->print();

            fem->calculate_element_stiffness(1);
            // fem->kk_element->print();
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk1->data, 1e-12));

            fem->calculate_element_stiffness(2);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk2->data, 1e-12));

            fem->calculate_element_stiffness(3);
            CHECK(equal_vectors_tol(fem->kk_element->data, correct_kk3->data, 1e-12));

            // solve the linear system
            fem->solve();

            // check solution
            CHECK(equal_vectors_tol(fem->uu, correct_uu, 1e-15));
        }

        SUBCASE("expanded (upper triangle) Bᵀ ⋅ D ⋅ B") {
            // allocate fem solver
            auto use_expanded_bdb = true;
            auto use_expanded_bdb_full = false;
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

            // check element stiffness
            fem->calculate_element_stiffness(0);
            auto kk0_upper = correct_kk0->get_copy();
            for (size_t i = 1; i < 6; i++) {
                for (size_t j = 0; j < i; j++) {
                    kk0_upper->set(i, j, 0.0);
                }
            }
            CHECK(equal_vectors_tol(fem->kk_element->data, kk0_upper->data, 1e-12));
            fem->kk_element->print();

            fem->calculate_element_stiffness(1);
            auto kk1_upper = correct_kk1->get_copy();
            for (size_t i = 1; i < 6; i++) {
                for (size_t j = 0; j < i; j++) {
                    kk1_upper->set(i, j, 0.0);
                }
            }
            CHECK(equal_vectors_tol(fem->kk_element->data, kk1_upper->data, 1e-12));

            fem->calculate_element_stiffness(2);
            auto kk2_upper = correct_kk2->get_copy();
            for (size_t i = 1; i < 6; i++) {
                for (size_t j = 0; j < i; j++) {
                    kk2_upper->set(i, j, 0.0);
                }
            }
            CHECK(equal_vectors_tol(fem->kk_element->data, kk2_upper->data, 1e-12));

            fem->calculate_element_stiffness(3);
            auto kk3_upper = correct_kk3->get_copy();
            for (size_t i = 1; i < 6; i++) {
                for (size_t j = 0; j < i; j++) {
                    kk3_upper->set(i, j, 0.0);
                }
            }
            CHECK(equal_vectors_tol(fem->kk_element->data, kk3_upper->data, 1e-12));

            // solve the linear system
            fem->solve();

            // check solution
            CHECK(equal_vectors_tol(fem->uu, correct_uu, 1e-15));
        }
    }

    SUBCASE("plane-strain domain (Smith's Example 5.2)") {
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

        // input data
        auto solid_triangle = true;
        auto plane_stress = false;
        auto thickness = 1.0;
        auto use_expanded_bdb = false;
        auto use_expanded_bdb_full = false;
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
        auto connectivity = vector<size_t>{
            1, 0, 3,  // 0
            3, 4, 1,  // 1
            2, 1, 4,  // 2
            4, 5, 2,  // 3
            4, 3, 6,  // 4
            6, 7, 4,  // 5
            5, 4, 7,  // 6
            7, 8, 5}; // 7
        auto param_young = vector<double>(8, 1e6);
        auto param_poisson = vector<double>(8, 0.3);
        auto param_cross_area = vector<double>{};
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

        // allocate fem solver
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

        // solve the linear system
        fem->solve();
        // print_vector("uu", fem->uu);

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
        CHECK(equal_vectors_tol(fem->uu, correct_uu, 1e-15));
    }
}
