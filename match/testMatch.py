import numpy as np
import sys
import h5py
import illustris_python as il

f_tng = h5py.File('/Raid1/Illustris/TNG-100/snap_ics.hdf5', 'r')
f_il1 = h5py.File('/Raid1/Illustris/Illustris-1/snap_ics.hdf5', 'r')
ids_tng = f_tng['PartType1']['ParticleIDs']
ids_il1 = f_il1['PartType1']['ParticleIDs']

matchedID = np.load('/Raid0/zhouzb/match/TNG_match_dict.npy', allow_pickle=True).item()

groupLen = il.func.loadhalos('TNG', 99, 'GroupLenType')[:, 1]


il1_pf = h5py.File('/Raid0/zhouzb/match/il1_haloPtcl.hdf5', 'r')
tng_pf = h5py.File('/Raid0/zhouzb/match/tng_haloPtcl.hdf5', 'r')
print('Load data successed.')

Pid_Matched = {}
for haloID in matchedID.keys():
    tng_PidIndex = tng_pf['diskGID']['%s' % haloID]
    Pid_Matched[haloID] = {} 
    for il1_halo in matchedID[haloID]:
        il1_haloPtcl = il.func.loadhalo('il1', 135, il1_halo, 1, fields='ParticleIDs')
        #count matched particle number
        count = 0
        for pid in tng_PidIndex:
            if ids_il1[pid] in il1_haloPtcl:
                count += 1
        Pid_Matched[haloID][il1_halo] = count / len(tng_PidIndex)
    print('Halo %d finisded'%haloID)

np.save('/Raid0/zhouzb/match/PidMatch_tng_to_il1.npy', Pid_Matched)
