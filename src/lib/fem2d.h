#pragma once

#include "laclib.h"
#include <map>
#include <memory>
#include <tuple>
#include <vector>

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

    /// @brief Properties = E*A (size = number_of_elements)
    std::vector<double> properties;

    /// @brief Essential (displacement) prescribed? (size = total_ndof)
    std::vector<bool> essential_prescribed;

    /// @brief Essential (displacement) boundary conditions (size = total_ndof)
    std::vector<double> essential_boundary_conditions;

    /// @brief Natural (force) boundary conditions (size = total_ndof)
    std::vector<double> natural_boundary_conditions;

    /// @brief Element stiffness matrix (4 x 4)
    std::unique_ptr<Matrix> kk_element;

    /// @brief Global displacements (size = total_ndof)
    std::vector<double> uu;

    /// @brief Right-hand side vector = global forces, corrected for prescribed displacements (size = total_ndof)
    std::vector<double> rhs;

    /// @brief Global stiffness matrix in COO format (nnz = 10 * number_of_elements)
    std::unique_ptr<CooMatrix> kk_coo;

    /// @brief Global stiffness matrix in CSR format (nnz = 10 * number_of_elements)
    std::unique_ptr<CsrMatrixMkl> kk_csr;

    /// @brief Holds the linear system solver
    std::unique_ptr<SolverDss> lin_sys_solver;

    /// @brief Allocates a new Truss2D structure
    /// @param coordinates x0 y0  x1 y1  ...  xnn ynn (size = 2 * number_of_nodes)
    /// @param connectivity 0 1  0 2  1 2  (size = 2 * number_of_elements)
    /// @param properties E*A (size = number_of_elements)
    /// @param essential_bcs prescribed boundary conditions. maps (node_number,dof_number) => value
    /// @param natural_bcs natural boundary conditions. maps (node_number,dof_number) => value
    inline static std::unique_ptr<Fem2d> make_new(
        const std::vector<double> &coordinates,
        const std::vector<size_t> &connectivity,
        const std::vector<double> &properties,
        const std::map<node_dof_pair_t, double> &essential_bcs,
        const std::map<node_dof_pair_t, double> &natural_bcs) {

        auto number_of_nodes = coordinates.size() / 2;
        auto number_of_elements = connectivity.size() / 2;
        auto total_ndof = 2 * number_of_nodes;

        // The number 10 below corresponds to the number of values in the
        // element stiffness matrix on the diagonal and above the diagonal (upper triangle)
        // The actual number of non-zeros is less than 10 * number_of_elements, so
        // this could be optimized by doing an assembly first
        auto nnz_max = 10 * number_of_elements;

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
            number_of_nodes,
            number_of_elements,
            total_ndof,
            coordinates,
            connectivity,
            properties,
            essential_prescribed,
            essential_boundary_conditions,
            natural_boundary_conditions,
            Matrix::make_new(4, 4),
            std::vector<double>(total_ndof),
            std::vector<double>(total_ndof),
            CooMatrix::make_new(layout, total_ndof, nnz_max),
            NULL,
            SolverDss::make_new(options),
        }};
    }

    /// @brief Returns the distance between node a and b
    /// @param a index of node in 0 <= a < number_of_nodes
    /// @param b index of node in 0 <= b < number_of_nodes
    inline double calculate_length(size_t a, size_t b) const {
        if (a >= coordinates.size()) {
            throw "cannot calculate the length because the a-th node index is out-of-range";
        }
        if (b >= coordinates.size()) {
            throw "cannot calculate the length because the b-th node index is out-of-range";
        }
        auto xa = coordinates[a * 2];
        auto ya = coordinates[a * 2 + 1];
        auto xb = coordinates[b * 2];
        auto yb = coordinates[b * 2 + 1];
        auto dx = xb - xa;
        auto dy = yb - ya;
        return sqrt(dx * dx + dy * dy);
    }

    /// @brief Calculates the element stiffness
    /// @param e index of element (rod) in 0 <= e < number_of_elements
    void calculate_element_stiffness(size_t e);

    /// @brief Calculates the global stiffness
    void calculate_rhs_and_global_stiffness();

    /// @brief Solves the mechanical problem
    void solve();
};
