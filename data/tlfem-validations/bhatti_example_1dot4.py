import numpy as np
from tlfem.mesh import Mesh
from tlfem.solver import Solver
from tlfem.util import CompareArrays

# Bhatti's Example 1.4 on page 25
#
# Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
#
# TEST GOAL
#
# This test verifies a 2D frame with rod elements and concentrated forces
#
# MESH
#
#               (3)
#               [2]
#     2----------------------3
#     |'.  (4)           _.-'
#     |  '.[3]       _.-'
#     |    '.    _.-'  (1)
# (2) |      '1-'      [1]
# [2] |      /
#     |     /
#     |    / (0)   The lines are ROD (Lin2) elements
#     |   /  [1]
#     |  /
#     | /    (#) indicates cell id
#     0'     [#] indicates attribute id
#
# BOUNDARY CONDITIONS
#
# Fully fixed @ points 0 and 3
# Concentrated load @ point 1 with Fy = -150,000
#
# CONFIGURATION AND PARAMETERS
#
# Static simulation
# Attribute 1: Area = 4,000; Young = 200,000
# Attribute 2: Area = 3,000; Young = 200,000
# Attribute 3: Area = 2,000; Young =  70,000

# vertices
V = [
    [0, -100, 0.0, 0.0],
    [1, -200, 1500, 3500],
    [2, 0, 0.0, 5000],
    [3, -300, 5000, 5000],
]
# cells/bars
C = [
    [0, -1, [0, 1]],
    [1, -1, [1, 3]],
    [2, -2, [0, 2]],
    [3, -2, [2, 3]],
    [4, -3, [1, 2]],
]
m = Mesh(V, C)
# m.draw()
# m.show()

# create dictionary with all parameters
p = {
    -1: {"type": "EelasticRod", "E": 200000, "A": 4000},
    -2: {"type": "EelasticRod", "E": 200000, "A": 3000},
    -3: {"type": "EelasticRod", "E": 70000, "A": 2000},
}

# allocate fem solver object
s = Solver(m, p)

# set boundary conditions
vertex_bcs = {
    -100: {"ux": 0.0, "uy": 0.0},
    -200: {"fy": -150000},
    -300: {"ux": 0.0, "uy": 0.0},
}
s.set_bcs(vb=vertex_bcs)

# solve steady
o = s.solve_steady(reactions=True)

# output
o.print_res("U")
o.print_res("R", True)

# bhatti's solution
Ubh = np.array(
    [
        [0.000000000000000e00, 0.000000000000000e00],
        [5.389536380057676e-01, -9.530613006371175e-01],
        [2.647036149579491e-01, -2.647036149579490e-01],
        [0.000000000000000e00, 0.000000000000000e00],
    ]
)
ux_fem = [o.Uout["ux"][n][-1] for n in range(m.nv)]
uy_fem = [o.Uout["uy"][n][-1] for n in range(m.nv)]
CompareArrays(ux_fem, Ubh[:, 0], tol=1e-15)
CompareArrays(uy_fem, Ubh[:, 1], tol=1e-15)
