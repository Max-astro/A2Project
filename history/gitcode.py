basePath = '/Raid1/Illustris/TNG/'
import numpy as np
from illustris_python.snapshot import loadSubhalo
from illustris_python.groupcat import loadSubhalos

def specific_angular_momentum(x, v, m):
    """
    specific angular momentum of a group of particles
    
    Parameters
    ----------
    x : array_like
        array particle positions of shape (Nptcl, ndim)
    v : array_like
        array of particle velcities wth shape (Nptcl, ndim)
    m : array_like
        array of particle masses of shape (Nptcl,)
    Returns
    -------
    L : nump.array
        specific angular momentum vector
    """
    return (m[:,np.newaxis]*np.cross(x,v)).sum(axis=0)

def galaxy_ang_mom(gal_id, basePath, snapNum, reduced=True):
    """
    Parameters
    ----------
    gal_id : int
    basepath : string
    snapNum : int
    Lbox : array_like
    reduced : bool
    Returns
    -------
    eig_vals, eig_vecs
    """

    # load galaxy position (most bound particle)
    gal_positions = loadSubhalos(basePath, snapNum, fields=['SubhaloPos'])/1000.0
    gal_position = gal_positions[gal_id]

    # half mass radius
    gal_rhalfs = loadSubhalos(basePath, snapNum, fields=['SubhaloHalfmassRadType'])[:,4]/1000.0
    gal_rhalf = gal_rhalfs[gal_id]

    # load stellar particles
    ptcl_coords = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Coordinates'])/1000.0
    ptcl_masses = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Masses'])*10.0**10
    ptcl_vels = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Velocities'])
    sf_time = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['GFM_StellarFormationTime'])
    is_a_star = (sf_time>=0.0) # don't use wind particles

    # account for PBCs
    dx = ptcl_coords[:,0] - gal_position[0]
    dy = ptcl_coords[:,1] - gal_position[1]
    dz = ptcl_coords[:,2] - gal_position[2]

    ptcl_coords = np.vstack((dx,dy,dz)).T

    r = np.sqrt(np.sum(ptcl_coords**2, axis=1))/gal_rhalf
    mask = (r<=10.0) & (is_a_star)

    L = specific_angular_momentum(ptcl_coords[mask], ptcl_vels[mask], ptcl_masses[mask])

    mag_L = np.sqrt(np.sum(L**2,axis=-1))

    return L, mag_L, L/mag_L