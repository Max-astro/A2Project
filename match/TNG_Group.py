import numpy as np
import sys
import h5py
import illustris_python as il

f = h5py.File('/Raid1/Illustris/TNG/snap_ics.hdf5', 'r')
ids = np.array(f['PartType1']['ParticleIDs'])
id_dict = dict(zip(ids, range(len(ids))))

matchedID = np.load('/Raid0/zhouzb/match/tng_MatchedGID.npy', allow_pickle=True)
GID = np.load('/Raid0/zhouzb/match/tng_diskGIDs.npy', allow_pickle=True)

groupLen_tng = il.func.loadhalos('TNG', 99, 'GroupLenType')[:, 1]
groupLen_il1 = il.func.loadhalos('il1', 135, 'GroupLenType')[:, 1]

il1_h5 = h5py.File('/Raid0/zhouzb/match/tng_haloPtcl.hdf5', 'a')

il1_h5.create_group('diskGID')
for haloID in GID:
    if groupLen[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('TNG', 99, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, 10000)
    pt_Index = []
    for i in PtcIDs:
        pt_Index.append(id_dict[i])
    il1_h5['diskGID']['%d'%haloID] = pt_Index

il1_h5.create_group('MatchID')
for haloID in matchedID:
    if groupLen_il1[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('TNG', 99, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, 10000)
    pt_Index = []
    for i in PtcIDs:
        pt_Index.append(id_dict[i])
    il1_h5['MatchID']['%d'%haloID] = pt_Index

il1_h5.close()