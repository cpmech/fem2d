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

    SUBCASE("read_mesh works") {
        auto mesh = read_mesh(data_path + "felippa_three_member_truss.msh");

        print_vector("coordinates", mesh->coordinates);
        auto correct_coo = vector<double>{
            0.0, 0.0,   // 0
            10.0, 0.0,  // 1
            10.0, 10.0, // 2
        };
        auto correct_con = vector<size_t>{
            0, 1, // 0
            1, 2, // 1
            2, 0, // 2
        };
        CHECK(equal_vectors_tol(mesh->coordinates, correct_coo, 1e-15));
        CHECK(equal_vectors(mesh->connectivity, correct_con));
    }
}
