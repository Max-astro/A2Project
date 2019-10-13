import numpy as np
import sys
import h5py
import illustris_python as il

f = h5py.File('/Raid1/Illustris/Illustris-1/snap_ics.hdf5', 'r')
ids = np.array(f['PartType1']['ParticleIDs'])
id_dict = dict(zip(ids, range(len(ids))))

matchedID = np.load('/Raid0/zhouzb/match/il1_MatchedGID.npy', allow_pickle=True)
# match_dict = np.load('/Raid0/zhouzb/match/TNG_match_dict.npy', allow_pickle=True).item()
GID = np.load('/Raid0/zhouzb/match/il1_diskGIDs.npy', allow_pickle=True)

groupLen = il.func.loadhalos('il1', 135, 'GroupLenType')[:, 1]

il1_h5 = h5py.File('/Raid0/zhouzb/match/il1_haloPtcl.hdf5', 'a')

il1_h5.create_group('diskGID')
for haloID in GID:
    if groupLen[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('il1', 135, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, 10000)
    pt_Index = []
    for i in PtcIDs:
        pt_Index.append(id_dict[i])
    il1_h5['diskGID']['%d'%haloID] = pt_Index



il1_h5.create_group('MatchID')
for haloID in matchedID:
    if groupLen[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('il1', 135, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, 10000)
    pt_Index = []
    for i in PtcIDs:
        pt_Index.append(id_dict[i])
    il1_h5['MatchID']['%d'%haloID] = pt_Index

il1_h5.close()