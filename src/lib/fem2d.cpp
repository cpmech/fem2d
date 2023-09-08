#include "fem2d.h"
#include "constants.h"
#include "fem2d.h"
#include "laclib.h"
#include "linear_elasticity.h"

void Fem2d::calculate_element_stiffness_elastic_rod(size_t e) {
    size_t a = connectivity[e * 2];
    size_t b = connectivity[e * 2 + 1];
    double xa = coordinates[a * 2];
    double ya = coordinates[a * 2 + 1];
    double xb = coordinates[b * 2];
    double yb = coordinates[b * 2 + 1];
    double dx = xb - xa;
    double dy = yb - ya;
    double l = sqrt(dx * dx + dy * dy);
    double c = (xb - xa) / l;
    double s = (yb - ya) / l;
    double p = param_young[e] * param_cross_area[e] / l;

    // computing upper triangle only
    //      _                   _
    //     |  c*c c*s -c*c -c*s  | 0
    // E A |   .  s*s -c*s -s*s  | 1
    // --- |   .   .   c*c  c*s  | 2
    //  L  |_  .   .    .   s*s _| 3
    //         0   1    2    3
    kk_element->set(0, 0, p * c * c);
    kk_element->set(0, 1, p * c * s);
    kk_element->set(0, 2, -p * c * c);
    kk_element->set(0, 3, -p * c * s);

    kk_element->set(1, 1, p * s * s);
    kk_element->set(1, 2, -p * c * s);
    kk_element->set(1, 3, -p * s * s);

    kk_element->set(2, 2, p * c * c);
    kk_element->set(2, 3, p * c * s);

    kk_element->set(3, 3, p * s * s);
}

void Fem2d::calculate_element_stiffness_solid_triangle(size_t e) {
    if (e >= number_of_elements) {
        throw "cannot calculate element stiffness because the element index is out-of-range";
    }

    // calculate elastic modulus
    linear_elasticity_modulus(dd, param_young[e], param_poisson[e], plane_stress);

    // auxiliary data
    size_t a = connectivity[e * 3];
    size_t b = connectivity[e * 3 + 1];
    size_t c = connectivity[e * 3 + 2];
    double x0 = coordinates[a * 2];
    double y0 = coordinates[a * 2 + 1];
    double x1 = coordinates[b * 2];
    double y1 = coordinates[b * 2 + 1];
    double x2 = coordinates[c * 2];
    double y2 = coordinates[c * 2 + 1];
    double a0 = y1 - y2;
    double a1 = y2 - y0;
    double a2 = y0 - y1;
    double b0 = x2 - x1;
    double b1 = x0 - x2;
    double b2 = x1 - x0;
    double f0 = x1 * y2 - x2 * y1;
    double f1 = x2 * y0 - x0 * y2;
    double f2 = x0 * y1 - x1 * y0;
    double area = (f0 + f1 + f2) / 2.0;
    double r = 2.0 * area;
    double s = r * SQRT_2;

    // element B-matrix (plane-strain and plane-stress)
    bb->set(0, 0, a0 / r);
    bb->set(0, 1, 0.0);
    bb->set(0, 2, a1 / r);
    bb->set(0, 3, 0.0);
    bb->set(0, 4, a2 / r);
    bb->set(0, 5, 0.0);
    bb->set(1, 0, 0.0);
    bb->set(1, 1, b0 / r);
    bb->set(1, 2, 0.0);
    bb->set(1, 3, b1 / r);
    bb->set(1, 4, 0.0);
    bb->set(1, 5, b2 / r);
    bb->set(2, 0, 0.0);
    bb->set(2, 1, 0.0);
    bb->set(2, 2, 0.0);
    bb->set(2, 3, 0.0);
    bb->set(2, 4, 0.0);
    bb->set(2, 5, 0.0);
    bb->set(3, 0, b0 / s);
    bb->set(3, 1, a0 / s);
    bb->set(3, 2, b1 / s);
    bb->set(3, 3, a1 / s);
    bb->set(3, 4, b2 / s);
    bb->set(3, 5, a2 / s);

    // K = Bᵀ ⋅ D ⋅ B ⋅ thickness ⋅ area
    mat_t_mat_mul(bb_t_dd, 1.0, bb, dd);
    mat_mat_mul(kk_element, thickness * area, bb_t_dd, bb);
}

void Fem2d::calculate_rhs_and_global_stiffness() {
    // The linear system is partitioned into unknown (1) and
    // prescribed (2) sub-matrices and sub-vectors
    //
    //  _             _
    // |  [K11] [K12]  | / {u1} \   / {f1} \  << unknown displacements
    // |               | |      | = |      |
    // |_ [K21] [K22] _| \ {u2} /   \ {f2} /  << prescribed displacements
    //       ^     ^
    // unknown     prescribed
    //
    // Note that:
    //
    // [K11]{u1} = {f1} - [K12]{u2}
    //
    // We need to solve the modified system
    //
    //  _             _
    // |  [K11]  [0]   | / {?1} \   / {f1} - [K12]{u2} \  << unknown displacements
    // |               | |      | = |                  |
    // |_  [0]   [1]  _| \ {?2} /   \       {u2}       /  << prescribed displacements
    //
    // {rhs1} = {f1} - [K12]{u2}
    // {rhs2} = {u2}

    // initialize uu and right-hand side vector
    // also, put ones on the diagonal of the global stiffness matrix
    kk_coo->pos = 0; // reset position
    for (size_t i = 0; i < total_ndof; ++i) {
        if (essential_prescribed[i]) {
            uu[i] = essential_boundary_conditions[i];  // {u2}: needed to correct RHS vector later on
            rhs[i] = essential_boundary_conditions[i]; // {rhs2}: because diagonal(K;prescribed) = 1
            kk_coo->put(i, i, 1.0);                    // [K22]: set diagonal(K;prescribed) = 1
        } else {
            uu[i] = 0.0;                             // {?1}: irrelevant, actually
            rhs[i] = natural_boundary_conditions[i]; // {rhs1} := {f1}, external forces
        }
    }

    // number of rows = number of columns in the element matrix
    size_t nrow = solid_triangle ? 6 : 4;

    // fix RHS vector and assemble stiffness
    for (size_t e = 0; e < number_of_elements; ++e) {
        calculate_element_stiffness(e);
        if (solid_triangle) {
            size_t a = connectivity[e * 3];
            size_t b = connectivity[e * 3 + 1];
            size_t c = connectivity[e * 3 + 2];
            m[0] = a * 2;
            m[1] = a * 2 + 1;
            m[2] = b * 2;
            m[3] = b * 2 + 1;
            m[4] = c * 2;
            m[5] = c * 2 + 1;
        } else {
            size_t a = connectivity[e * 2];
            size_t b = connectivity[e * 2 + 1];
            m[0] = a * 2;
            m[1] = a * 2 + 1;
            m[2] = b * 2;
            m[3] = b * 2 + 1;
        }
        for (size_t i = 0; i < nrow; ++i) {
            if (!essential_prescribed[m[i]]) {
                // {rhs1} -= [K12]{u2}, correct RHS
                for (size_t j = 0; j < nrow; ++j) {
                    if (essential_prescribed[m[j]]) {
                        if (j >= i) {
                            rhs[m[i]] -= kk_element->get(i, j) * uu[m[j]];
                        } else {
                            // must get (i,j) from upper triangle
                            rhs[m[i]] -= kk_element->get(j, i) * uu[m[j]];
                        }
                    }
                }
                // [K11]: assemble upper triangle into global stiffness
                for (size_t j = i; j < nrow; ++j) { // j = i => local upper triangle
                    if (!essential_prescribed[m[j]]) {
                        if (m[j] >= m[i]) {
                            kk_coo->put(m[i], m[j], kk_element->get(i, j));
                        } else {
                            // must go to the global upper triangle
                            kk_coo->put(m[j], m[i], kk_element->get(i, j));
                        }
                    }
                }
            }
        }
    }

    // convert COO to CSR
    if (kk_csr == NULL) {
        kk_csr = CsrMatrixMkl::from(kk_coo);
    } else {
        kk_csr.reset();
        kk_csr = CsrMatrixMkl::from(kk_coo);
    }
}

void Fem2d::solve() {
    if (kk_csr == NULL) {
        calculate_rhs_and_global_stiffness();
    }
    lin_sys_solver->analyze(kk_csr);
    lin_sys_solver->factorize(kk_csr);
    lin_sys_solver->solve(uu, rhs); // uu = inv(kk) * ff
}
