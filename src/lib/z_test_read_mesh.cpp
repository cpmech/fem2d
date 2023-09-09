#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "../util/doctest.h"
#include "laclib.h"
#include "read_mesh.h"
#include <vector>

#ifndef DATA_DIR
#define DATA_DIR "data"
#endif

using namespace std;

#define _SUBCASE(name) if (false)

TEST_CASE("read_mesh") {

    auto data_path = string(DATA_DIR) + "/meshes/";

    SUBCASE("read_mesh works (felippa_three_member_truss)") {
        auto mesh = read_mesh(data_path + "felippa_three_member_truss.msh");
        auto correct_coo = vector<double>{
            0.0, 0.0,    // 0
            10.0, 0.0,   // 1
            10.0, 10.0}; // 2
        auto correct_con = vector<size_t>{
            0, 1,  // 0
            1, 2,  // 1
            2, 0}; // 2
        CHECK(equal_vectors_tol(mesh->coordinates, correct_coo, 1e-15));
        CHECK(equal_vectors(mesh->connectivity, correct_con));
    }

    SUBCASE("read_mesh works (smith_plane_strain_5dot2)") {
        auto mesh = read_mesh(data_path + "smith_plane_strain_5dot2.msh");
        auto correct_coo = vector<double>{
            0.0, 0.0,   // 0
            0.5, 0.0,   // 1
            1.0, 0.0,   // 2
            0.0, -0.5,  // 3
            0.5, -0.5,  // 4
            1.0, -0.5,  // 5
            0.0, -1.0,  // 6
            0.5, -1.0,  // 7
            1.0, -1.0}; // 8
        auto correct_con = vector<size_t>{
            1, 0, 3,  // 0
            3, 4, 1,  // 1
            2, 1, 4,  // 2
            4, 5, 2,  // 3
            4, 3, 6,  // 4
            6, 7, 4,  // 5
            5, 4, 7,  // 6
            7, 8, 5}; // 7
        CHECK(equal_vectors_tol(mesh->coordinates, correct_coo, 1e-15));
        CHECK(equal_vectors(mesh->connectivity, correct_con));
    }
}
