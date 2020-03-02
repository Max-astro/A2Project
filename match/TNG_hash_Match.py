import numpy as np
import sys
import h5py
import illustris_python as il

f_tng = h5py.File('/Raid1/Illustris/TNG/snap_ics.hdf5', 'r')
f_il1 = h5py.File('/Raid1/Illustris/Illustris-1/snap_ics.hdf5', 'r')
ids_tng = np.array(f_tng['PartType1']['ParticleIDs'])
ids_il1 = np.array(f_il1['PartType1']['ParticleIDs'])
id_dict = dict(zip(ids_tng, ids_il1))
ids_tng = 0
ids_il1 = 0

matchedID = np.load('/Raid0/zhouzb/match/TNG_match_dict.npy', allow_pickle=True).item()

groupLen = il.func.loadhalos('TNG', 99, 'GroupLenType')[:, 1]

Pid_Matched = {}


for haloID in matchedID.keys():
    if groupLen[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('TNG', 99, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, int(len(PtcIDs) / 10))
    
    for i in PtcIDs:
        pt_Index.append(id_dict[i])


il1_h5.create_group('MatchID')
for haloID in matchedID:
    if groupLen[haloID] < 100000:
        continue
    PtcIDs = il.func.loadhalo('TNG', 99, haloID, 1, 'ParticleIDs')
    PtcIDs = np.random.choice(PtcIDs, 10000)
    pt_Index = []
    for i in PtcIDs:
        pt_Index.append(id_dict[i])
    il1_h5['MatchID']['%d'%haloID] = pt_Index

il1_h5.close()