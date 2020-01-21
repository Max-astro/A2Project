#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import h5py
import sys
sys.path.append('F:\Linux')
import illustris_python as il
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt

def find_id(all_id):
    il1_id = np.zeros(len(all_id), dtype=int)
    for i in range(len(all_id)):
        il1_id[i] = -1
        for j in range(5):
            if all_id[i][j] != -1:
                il1_id[i] = all_id[i][j]
    return il1_id

def SubhaloMatch(tng_subID):
    '''
    return matched SubhaloID in il1.
    No match return -1
    '''
    tng_GID = Groups_tng[tng_subID]
    il1_GID = tng_to_il1[tng_GID]
    if il1_GID == -1:
        return il1_GID
    try:
        if len(il1_HaloToSubdisk[il1_GID]) == 1:
            return (il1_HaloToSubdisk[il1_GID])[0]
        il1_subIDs = il1_HaloToSubdisk[il1_GID]
    except KeyError:
        il1_subIDs = il1_GroupSubhalo[il1_GID]

    return SubPositionMatch(tng_subID, il1_subIDs)
    
def SubPositionMatch(tng_subID, il1_subIDs):
    '''
    If one dark matter halo has several center galaxies, match them by subhalo position.
    return a list of their distance
    '''
    mass = abs(tng_DMMass[tng_subID] - il1_DMMass[il1_subIDs])
    ind = np.where(mass == mass.min())[0][0]
    return il1_subIDs[ind]


# In[2]:


#match disk halo between TNG & illustris-1
tngdisk = np.load('f:/Linux/localRUN/diskID_TNG.npy')
il1disk = np.load('f:/Linux/localRUN/diskID_il1.npy')

tng_barred = np.load('f:/Linux/localRUN/barredID_TNG.npy')
il1_barred = np.load('f:/Linux/localRUN/barredID_il1.npy')

Groups_tng = il.func.loadSubhalos('TNG', 99, "SubhaloGrNr")
Groups_il1 = il.func.loadSubhalos('il1', 135, "SubhaloGrNr")

tng_diskGIDs = Groups_tng[tngdisk]
il1_diskGIDs = Groups_il1[il1disk]

tng_barredGID = Groups_tng[tng_barred]
il1_barredGID = Groups_il1[il1_barred]

#load TNG groups informations
cm_TNG = il.func.loadhalos('TNG', 99, 'GroupPos')
tng_subCM = il.func.loadSubhalos('TNG', 99, 'SubhaloCM')
mass_TNG = il.func.loadhalos('TNG', 99, 'GroupMassType')[:, 1]

#build illustris halo as a tree
cm_il1 = il.func.loadhalos('il1', 135, 'GroupPos')
il1_subCM = il.func.loadSubhalos('il1', 135, 'SubhaloCM')
mass_il1 = il.func.loadhalos('il1', 135, 'GroupMassType')[:, 1]

#Is centre galaxy
tng_parent = il.func.loadSubhalos('TNG', 99, "SubhaloParent")
il1_parent = il.func.loadSubhalos('il1', 135, "SubhaloParent")
tng_centrediskID = tngdisk[tng_parent[tngdisk] == 0]
il1_centrediskID = il1disk[il1_parent[il1disk] == 0]

#SFR
tng_SFR = il.func.loadSubhalos('TNG', 99, "SubhaloSFRinHalfRad")
il1_SFR = il.func.loadSubhalos('il1', 135, "SubhaloSFRinHalfRad")

#DM mass
tng_DMMass = il.func.loadSubhalos('TNG', 99, "SubhaloMassInRadType")[:,1]
il1_DMMass = il.func.loadSubhalos('il1', 135, "SubhaloMassInRadType")[:,1]

#Stellar mass
tng_StellarMass = il.func.loadSubhalos('TNG', 99, "SubhaloMassInRadType")[:,4]
il1_StellarMass = il.func.loadSubhalos('il1', 135, "SubhaloMassInRadType")[:,4]

#Gas Fraction
tng_GasFraction = il.func.loadSubhalos('TNG', 99, "SubhaloMassInRadType")[:,0] / il.func.loadSubhalos('TNG', 99, "SubhaloMassInHalfRad")
il1_GasFraction = il.func.loadSubhalos('il1', 135, "SubhaloMassInRadType")[:,0] / il.func.loadSubhalos('il1', 135, "SubhaloMassInHalfRad")
tng_GasFraction[np.isnan(tng_GasFraction)] = 0
il1_GasFraction[np.isnan(il1_GasFraction)] = 0


# In[3]:


tree = KDTree(cm_il1, leaf_size = 2)


# In[4]:


#tngID, il1ID, dis, absMass
sameHalo = [[], [], [], []]
for groupID in np.unique(tng_diskGIDs):
    sameHalo[0].append([groupID])
    dis, index = tree.query(cm_TNG[groupID].reshape(1, -1), k=10)

    il1_id = []
    distance = []
    absMass = []
    for i in range(10):
        if abs(mass_TNG[groupID] - mass_il1[index[0, i]]) <= mass_TNG[groupID]*0.5:
            il1_id.append(index[0, i])
            distance.append(dis[0, i])
            absMass.append(abs(mass_TNG[groupID] - mass_il1[index[0, i]]))
        else:
            il1_id.append(-1)
            distance.append(dis[0, i])
            absMass.append(abs(mass_TNG[groupID] - mass_il1[index[0, i]]))

    sameHalo[1].append(il1_id)
    sameHalo[2].append(distance)
    sameHalo[3].append(absMass)

sameHalo = np.array(sameHalo)


# In[5]:


tng_to_il1 = np.load('f:/Linux/localRUN/tng_to_il1/tng_to_il1.npy').item()
tng_bar_hostGID_matched = np.load('f:/Linux/localRUN/tng_to_il1/tng_bar_hostGID_matched.npy')
tng_barred_matchedGID = np.load('f:/Linux/localRUN/tng_to_il1/tng_barred_matchedGID.npy')
tng_disk_matchedGID = np.load('f:/Linux/localRUN/tng_to_il1/tng_disk_matchedGID.npy')
tng_to_il1_matched_barred_centergalaxies = np.load('f:/Linux/localRUN/tng_to_il1/tng_to_il1_matched_barred_centergalaxies.npy').item()
tng_barred_il1_nobarGID = np.load('f:/Linux/localRUN/tng_to_il1/tng_barred_il1_nobarGID.npy')

#Build a dictionary il1_DMhaloID to il1_disk_SubhaloID, return a list
il1_HaloToSubdisk = {}
last = -1
for i in il1_centrediskID:
    tmp = Groups_il1[i]
    if tmp != last:
        il1_HaloToSubdisk[tmp] = [i]
    else:
        il1_HaloToSubdisk[tmp].append(i)
    last = tmp
    
#Build a dictionary for all il1_subhalo corresponding to their DMhalo, return a list
il1_GroupSubhalo = {}
last = -1
for i in range(len(Groups_il1)):
    if il1_parent[i] != 0:
        continue
    tmp = Groups_il1[i]
    if tmp != last:
        il1_GroupSubhalo[tmp] = [i]
    else:
        il1_GroupSubhalo[tmp].append(i)
    last = tmp

#tng DMhalo matched il1 DMhalo with no disk galaxy
tng_matched_nodiskGID=[]
for i in tng_to_il1.values():
    if (i != -1) & (i not in il1_diskGIDs):
        tng_matched_nodiskGID.append(i)


# In[6]:


#Dictionary of tng subhalo to il1 subhalo
tng_to_il1_subID = {}
for i in tng_barred:
    il1_sid = SubhaloMatch(i)
 
    tng_to_il1_subID[i] = SubhaloMatch(i)


# In[7]:


#TNG bar/unbar -- Gasfraction / StellatMass / SFR / AM
disk_SFR = tng_SFR[tngdisk]
bar_SFR = tng_SFR[tng_barred]

disk_Gas = tng_GasFraction[tngdisk]
bar_Gas = tng_GasFraction[tng_barred]

disk_SMass = tng_StellarMass[tngdisk]
bar_SMass = tng_StellarMass[tng_barred]


# In[13]:


#illustris-1 bar/unbar -- Gasfraction / StellatMass / SFR / AM
il1_disk_SFR = il1_SFR[il1disk]
il1_bar_SFR = il1_SFR[il1_barred]

il1_disk_Gas = il1_GasFraction[il1disk]
il1_bar_Gas = il1_GasFraction[il1_barred]

il1_disk_SMass = tng_StellarMass[il1disk]
il1_bar_SMass = il1_StellarMass[il1_barred]


# In[ ]:




