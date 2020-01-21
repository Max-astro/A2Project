import numpy as np
import h5py
import illustris_python as il
import matplotlib.pyplot as plt

#TNG slice

maxid = 0
minid = 0

for i in range(448):
    with h5py.File('/Raid1/Illustris/TNG-100/Snapshot/snap_99/snap_099.%d.hdf5'%i,'r') as f:
        ids = np.array(f['PartType1']['ParticleIDs'])
        maxid = ids.max() if (ids.max() > maxid) else maxid
        minid = ids.min() if (ids.min() < minid) else minid


print("maxid = ",maxid)
print("minid = ",minid)