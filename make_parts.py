#!/usr/bin/env /bin/python3
import numpy as np
from sfoda_utils import Grid

L = 124000
H = 957


def get_depth(x):
    a1 = 7.900666e+02
    b1 = 2.717665e-02
    c1 = -6.374170e-01
    a2 = 3.754591e+02
    b2 = 5.443289e-02
    c2 = 1.232468e+00
    a3 = 1.180416e+02
    b3 = 9.870790e-02
    c3 = 1.833616e+00
    a4 = 6.888052e+01
    b4 = 1.406913e-01
    c4 = 2.538620e+00
    a5 = 9.968110e+01
    b5 = 2.027221e-01
    c5 = 4.999398e+00
    a6 = 1.130123e+02
    b6 = 1.917498e-01
    c6 = 2.483625e+00
    a7 = 9.585347e+00
    b7 = 2.925321e-01
    c7 = 3.719423e+00
    a8 = 5.146615e+00
    b8 = 4.044217e-01
    c8 = 8.558507e-01

    z = a1 * np.sin(b1 *x/ 1000.0 + c1) + a2 * np.sin(b2 * x/ 1000.0 + c2) + a3 * np.sin(b3 * x/ 1000.0 + c3) + a4 * np.sin(b4 * x/ 1000.0 + c4) + a5 * np.sin(b5 * x/ 1000.0 + c5) + a6 * np.sin(b6 * x/ 1000.0 + c6) + a7 * np.sin(b7 * x/ 1000.0 + c7) + a8 * np.sin(b8 * x/ 1000.0 + c8)

    return 957 - z


parts_x = 50
parts_z = 20

margins_x = 2000
margins_z = 50

parts = np.zeros((parts_x*parts_z, 3), dtype=np.float64)

real_L = L - margins_x*2

step_x = real_L/parts_x

sungrid = Grid("data", VERBOSE=True, Nk=0)
sungrid.loadBathy("data/depth.dat-voro")
y1 = np.min(sungrid.yv)
y2 = np.max(sungrid.yv)

for i in range(parts_x):
    for j in range(parts_z):
        num = i*parts_z + j
        x = i*step_x + margins_x
        H = get_depth(x)
        if H > margins_z*2:
            real_H = H - margins_z*2
        else:
            continue
        step_z = real_H/parts_z
        parts[num][0] = x
        parts[num][1] = (y1 + y2) / 2
        parts[num][2] = j*step_z + margins_z

np.savetxt("data/2d_particles.dat", parts[parts[:, 0] != 0], fmt="%.6e %.6e %.6e")
