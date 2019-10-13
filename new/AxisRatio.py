import numpy as np
import h5py
import json
import sys
sys.path.append('F:\Linux')
import illustris_python as il
import inertia_tensors



def MassTensorEigenVals(coor, mas, half_r):
    '''
    Return eigenvalues of the mass tensor, sorted by M1 < M2 < M3
    '''
    r = coor - coor[0]
    r[r > 37500] -= 75000
    r[r < -37500] += 75000
    dis = np.linalg.norm(r, axis=1)
    inside = dis < (half_r * 2)
    r = r[inside]
    mas = mas[inside]

    M_x = ((mas * r[:, 0]**2).sum())**0.5 / mas.sum()**0.5
    M_y = ((mas * r[:, 1]**2).sum())**0.5 / mas.sum()**0.5
    M_z = ((mas * r[:, 2]**2).sum())**0.5 / mas.sum()**0.5
    M = np.array([M_x, M_y, M_z])
    M.sort()
    return M