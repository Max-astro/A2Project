import numpy as np
import sys
import h5py
import time
import illustris_python as il

f = h5py.File('/Raid1/Illustris/Illustris-1/snap_ics.hdf5', 'r')

ids = np.array(f['PartType1']['ParticleIDs'])
id_dict = dict(zip(ids, range(len(ids))))


data = np.load('/Raid0/zhouzb/match/tng_to_il1_GID.npy', allow_pickle=True).item()
il1_h5 = h5py.File('/Raid0/zhouzb/match/il1_haloPtcl.hdf5', 'a')
il1_h5.create_group('MatchID')
il1_h5.create_group('diskGID')
for haloID in data.keys():
