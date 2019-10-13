import numpy as np
import h5py
import sys
import illustris_python as il
    
# tngdisk = np.load('f:/Linux/localRUN/diskID_TNG.npy')
# il1disk = np.load('f:/Linux/localRUN/diskID_il1.npy')

# tng_barred = np.load('f:/Linux/localRUN/barredID_TNG.npy')
# il1_barred = np.load('f:/Linux/localRUN/barredID_il1.npy')

# il1_subCM = il.func.loadSubhalos('il1', 135, 'SubhaloCM')
# StellarHalfmassRads = (il.func.loadSubhalos('il1',135,'SubhaloHalfmassRadType'))[:,4]

def angular(subhaloID, simu='il1', snapnum=135, pt=4, rad=2):
    if pt == 4:
        masses = il.func.loadgalaxy(simu, snapnum, subhaloID, partType=pt, fields='Masses')
        coor = il.func.loadgalaxy(simu, snapnum, subhaloID, partType=pt, fields='Coordinates')
        vel = il.func.loadgalaxy(simu, snapnum, subhaloID, partType=pt, fields='Velocities')
        CM = il.func.loadSubhalos(simu, snapnum, 'SubhaloCM')[subhaloID]
        half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[subhaloID, pt]

    if pt == 4:
        sf_time = il.func.loadgalaxy(simu, snapnum, subhaloID, partType=pt, fields='GFM_StellarFormationTime')
        is_star = (sf_time >= 0.0)  # don't use wind particles

        vel = vel[is_star]
        masses = masses[is_star]
        coor = coor[is_star]

    r = coor - CM
    r[r > 37500] -= 75000
    r[r < -37500] += 75000
    if rad != 0:
        dis = ((r**2).sum(1))**0.5
        inside = dis < (half_r*rad)

        vel = vel[inside]
        masses = masses[inside]
        r = r[inside]

    V = np.sum(vel * masses[:, np.newaxis], 0) / masses.sum()
    vel = vel - V
    L = np.sum(np.cross(r, vel * masses[:, np.newaxis]), axis=0)

    return L


def ang(masses, coor, vel, cm):
    r = coor - cm

    r[r > 37500] -= 75000
    r[r < -37500] += 75000

    V = np.sum(vel * masses[:, np.newaxis], 0) / masses.sum()
    vel = vel - V
    L = np.sum(np.cross(r, vel * masses[:, np.newaxis]), axis=0)

    return L
