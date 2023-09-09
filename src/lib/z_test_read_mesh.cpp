#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "../util/doctest.h"
#include "laclib.h"
#include "read_mesh.h"

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
    }
}
