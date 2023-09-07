#pragma once

#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "laclib.h"

using namespace std;

/// @brief Defines the index of a local DOF (0 or 1)
enum LocalDOF {
    AlongX = 0,
    AlongY = 1,
};

/// @brief Holds the pair (node_number, dof_number)
typedef std::tuple<size_t, LocalDOF> node_dof_pair_t;

/// @brief Implements a finite element solver for trusses in 2D
struct Fem2d {
    /// @brief Simulate linear elastic solid triangles with 3 nodes instead of linear elastic rods with 2nodes
    bool solid_triangle;

    /// @brief If solid_triangle, simulate plane-stress instead of plane-strain
    bool plane_stress;

    /// @brief Thickness of the plate if plane-stress and solid_triangle (will be set to 1.0 if plane-stress = false)
    double thickness;

    /// @brief Holds the number of nodes = coordinates.size() / 2
    size_t number_of_nodes;

    /// @brief Holds the number of elements = connectivity.size() / 2
    size_t number_of_elements;

    /// @brief Holds the total number of DOFs = 2 * number_of_nodes
    size_t total_ndof;

    /// @brief Coordinates x0 y0  x1 y1  ...  xnn ynn (size = 2 * number_of_nodes)
    std::vector<double> coordinates;

    /// @brief Connectivity 0 1  0 2  1 2  (size = 2 * number_of_elements)
    std::vector<size_t> connectivity;

    /// @brief Holds all Young's modulus (size = number_of_elements)
    std::vector<double> param_young;

    /// @brief Holds all Poisson coefficients (solid_triangle only) (size = number_of_elements)
    std::vector<double> param_poisson;

    /// @brief Holds all cross-sectional areas (rod element only) (size = number_of_elements)
    std::vector<double> param_cross_area;

    /// @brief Essential (displacement) prescribed? (size = total_ndof)
    std::vector<bool> essential_prescribed;

    /// @brief Essential (displacement) boundary conditions (size = total_ndof)
    std::vector<double> essential_boundary_conditions;

    /// @brief Natural (force) boundary conditions (size = total_ndof)
    std::vector<double> natural_boundary_conditions;

    /// @brief Holds the strain-displacement B matrix (4 x 6; used with solid_triangle)
    std::unique_ptr<Matrix> bb;

    /// @brief Holds the elastic D matrix (4 x 4; used with solid_triangle)
    std::unique_ptr<Matrix> dd;

    /// @brief Holds the result of transpose(B) * D (6 x 4; used with solid_triangle)
    std::unique_ptr<Matrix> bb_t_dd;

    /// @brief Element stiffness matrix (6 x 6 if solid_triangle; 4 x 4 otherwise)
    std::unique_ptr<Matrix> kk_element;

    /// @brief maps local to global DOF
    std::vector<size_t> m;

    /// @brief Global displacements (size = total_ndof)
    std::vector<double> uu;

    /// @brief Right-hand side vector = global forces, corrected for prescribed displacements (size = total_ndof)
    std::vector<double> rhs;

    /// @brief Global stiffness matrix in COO format (nnz = (10 or 21) * number_of_elements)
    std::unique_ptr<CooMatrix> kk_coo;

    /// @brief Global stiffness matrix in CSR format (nnz = (10 or 21) * number_of_elements)
    std::unique_ptr<CsrMatrixMkl> kk_csr;

    /// @brief Holds the linear system solver
    std::unique_ptr<SolverDss> lin_sys_solver;

    /// @brief Allocates a new Truss2D structure
    /// @param solid_triangle Plane-stress or plane-strain analysis with triangles instead of frames in 2D
    /// @param thickness Out-of-plane thickness if solid-triangle and plane-stress
    /// @param plane_stress If solid-triangle, simulate plane-stress instead of plane-strain
    /// @param coordinates x0 y0  x1 y1  ...  xnn ynn (size = 2 * number_of_nodes)
    /// @param connectivity 0 1 (2)  0 2 (3)  1 2 (4)  (size = (2 or 3) * number_of_elements)
    /// @param param_young All Young's modulus (size = number_of_elements)
    /// @param param_poisson All Poisson coefficients (solid_triangle only) (size = number_of_elements)
    /// @param param_cross_area All cross-sectional areas (rod element only) (size = number_of_elements)
    /// @param essential_bcs prescribed boundary conditions. maps (node_number,dof_number) => value
    /// @param natural_bcs natural boundary conditions. maps (node_number,dof_number) => value
    inline static std::unique_ptr<Fem2d> make_new(bool solid_triangle,
                                                  bool plane_stress,
                                                  double thickness,
                                                  const std::vector<double> &coordinates,
                                                  const std::vector<size_t> &connectivity,
                                                  const std::vector<double> &param_young,
                                                  const std::vector<double> &param_poisson,
                                                  const std::vector<double> &param_cross_area,
                                                  const std::map<node_dof_pair_t, double> &essential_bcs,
                                                  const std::map<node_dof_pair_t, double> &natural_bcs) {

        size_t element_num_node = solid_triangle ? 3 : 2;
        auto number_of_nodes = coordinates.size() / 2;
        auto number_of_elements = connectivity.size() / element_num_node;
        auto total_ndof = 2 * number_of_nodes;

        // The number sum_band below corresponds to the number of values in the
        // element stiffness matrix on the diagonal and above the diagonal (upper triangle)
        // The actual number of non-zeros is less than 10 * number_of_elements, so
        // this could be optimized by doing an assembly first
        auto sum_band = solid_triangle ? 21 : 10;
        auto nnz_max = sum_band * number_of_elements;

        auto essential_prescribed = std::vector<bool>(total_ndof, false);
        auto essential_boundary_conditions = std::vector<double>(total_ndof, 0.0);
        auto natural_boundary_conditions = std::vector<double>(total_ndof, 0.0);

        for (const auto &[key, value] : essential_bcs) {
            const auto [node, dof] = key;
            auto global_dof = node * 2 + dof;
            essential_prescribed[global_dof] = true;
            essential_boundary_conditions[global_dof] = value;
        }

        for (const auto &[key, value] : natural_bcs) {
            const auto [node, dof] = key;
            auto global_dof = node * 2 + dof;
            natural_boundary_conditions[global_dof] = value;
        }

        StoredLayout layout = UPPER_TRIANGULAR;

        auto options = DssOptions::make_new();
        options->symmetric = true;
        options->positive_definite = true;

        return std::unique_ptr<Fem2d>{new Fem2d{
            solid_triangle,
            plane_stress,
            plane_stress ? thickness : 1.0,
            number_of_nodes,
            number_of_elements,
            total_ndof,
            coordinates,
            connectivity,
            param_young,
            param_poisson,
            param_cross_area,
            essential_prescribed,
            essential_boundary_conditions,
            natural_boundary_conditions,
            solid_triangle ? Matrix::make_new(4, 6) : NULL, // bb
            solid_triangle ? Matrix::make_new(4, 4) : NULL, // dd
            solid_triangle ? Matrix::make_new(6, 4) : NULL, // bb_t_dd
            solid_triangle ? Matrix::make_new(6, 6) : Matrix::make_new(4, 4),
            std::vector<size_t>(solid_triangle ? 6 : 4),
            std::vector<double>(total_ndof),
            std::vector<double>(total_ndof),
            CooMatrix::make_new(layout, total_ndof, nnz_max),
            NULL,
            SolverDss::make_new(options),
        }};
    }

    /// @brief Calculates the element stiffness (Elastic Rod)
    /// @param e index of element (rod) in 0 <= e < number_of_elements
    void calculate_element_stiffness_elastic_rod(size_t e);

    /// @brief Calculates the element stiffness (Solid Triangle)
    /// @param e index of element (rod) in 0 <= e < number_of_elements
    void calculate_element_stiffness_solid_triangle(size_t e);

    /// @brief Calculates the element stiffness
    /// @param e index of element (rod) in 0 <= e < number_of_elements
    inline void calculate_element_stiffness(size_t e) {
        if (solid_triangle) {
            calculate_element_stiffness_solid_triangle(e);
        } else {
            calculate_element_stiffness_elastic_rod(e);
        }
    }

    /// @brief Calculates the global stiffness
    void calculate_rhs_and_global_stiffness();

    /// @brief Solves the mechanical problem
    void solve();
};
