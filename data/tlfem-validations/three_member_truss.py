import numpy as np
from tlfem.mesh import Mesh
from tlfem.solver import Solver
from tlfem.util import CompareArrays

# GEOMETRY
#
#                      fy=1 ↑
# ---                       2 →
#  ↑                      ,'| fx=2
#  |                    ,'  |
#  |                  ,'    |
#  |       EA=200√2 ,'      |
# 10          (2) ,'        | EA=50
#  |            ,'          | (1)
#  |          ,'            |
#  |        ,'              |
#  |      ,'    EA=100      |
#  ↓    ,'       (0)        |
# ---  0--------------------1
#     | |                  | |
#      ⇊ uy=-0.5     uy=0.4 ⇈
#
#      |←------- 10 -------→|
#
# BOUNDARY CONDITIONS
#
# node 0: x-fixed with a vertical displacement: uy = -0.5
# node 1: x-fixed with a vertical displacement: uy = 0.4
# node 2: fx = 2.0 and fy = 1.0
#
# EXPECTED RESULTS
#
# kk * uu = ff
#
# correct_uu = [0.0, -0.5, 0.0, 0.4, -0.5, 0.2]
# correct_ff = [-2.0, -2.0, 0.0, 1.0, 2.0, 1.0]

# vertices
V = [
    [0, -100, 0.0, 0.0],
    [1, -200, 10.0, 0.0],
    [2, -300, 10.0, 10.0],
]
# elements (cells)
C = [
    [0, -1, [0, 1]],
    [1, -2, [1, 2]],
    [2, -3, [0, 2]],
]
m = Mesh(V, C)
# m.draw()
# m.show()

# parameters
p = {
    -1: {"type": "EelasticRod", "E": 100.0, "A": 1.0},
    -2: {"type": "EelasticRod", "E": 50.0, "A": 1.0},
    -3: {"type": "EelasticRod", "E": 200.0, "A": np.sqrt(2.0)},
}

# allocate fem solver object
s = Solver(m, p)

# set boundary conditions
vertex_bcs = {
    -100: {"ux": 0.0, "uy": -0.5},
    -200: {"ux": 0.0, "uy": 0.4},
    -300: {"fx": 2.0, "fy": 1.0},
}
s.set_bcs(vb=vertex_bcs)

# solve steady
o = s.solve_steady(reactions=True, output_F=True)

# output
o.print_res("U")
o.print_res("R", True)

# check solution
uu = np.array(
    [(o.Uout["ux"][n][-1], o.Uout["uy"][n][-1]) for n in range(m.nv)]
).flatten()
uu_exp = [0.0, -0.5, 0.0, 0.4, -0.5, 0.2]
ff_exp = [-2.0, -2.0, 0.0, 1.0, 2.0, 1.0]
CompareArrays(uu, uu_exp, tol=1e-15)
CompareArrays(o.F, ff_exp, tol=1e-15)
