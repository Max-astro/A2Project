import numpy as np
import h5py
import json
import sys
sys.path.append('/home/zhouzb/')
import illustris_python as il

def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if simu == 'TNG' or simu == 'TNG100-1':
        ldir = '/Raid0/zhouzb/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = '/Raid0/zhouzb/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)

    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])


il1_snap = [120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [84, 67, 59, 50, 40, 33, 29, 25]

tng_diskID = np.load('/Raid0/zhouzb/rundata/diskID_4WP.npy')
il1_diskID = np.load('/Raid0/zhouzb/rundata/diskID_il1.npy')

gidHist = {}
for snap in il1_snap:
    gidHist[snap] = il.func.loadSubhalos('il1', snap, 'SubhaloGrNr')
print("Load success.")

#output
count = 0
data = {}
for haloID in il1_diskID:
    data[haloID] = {}
    prog = LoadMergHist('il1', haloID)[0]
    for snap in il1_snap:
        try:
            subID = prog[snap]
        except:
            subID = -1

        gid = -1
        if subID != -1:
            gid = gidHist[snap][subID]
        
        data[haloID][snap] = gid
    count += 1
    print("Done: %d / %d"%(count, len(il1_diskID)))

np.save('/Raid0/zhouzb/il1_GrNr.npy', data)
        

#---------TNG---------
gidHist = {}
for snap in tng_snap:
    gidHist[snap] = il.func.loadSubhalos('TNG', snap, 'SubhaloGrNr')
print("Load success.")

#output
count = 0
data = {}
for haloID in tng_diskID:
    data[haloID] = {}
    prog = LoadMergHist('TNG', haloID)[0]
    for snap in tng_snap:
        try:
            subID = prog[snap]
        except:
            subID = -1
            
        gid = -1
        if subID != -1:
            gid = gidHist[snap][subID]
        
        data[haloID][snap] = gid
    count += 1
    print("Done: %d / %d"%(count, len(tng_diskID)))

np.save('/Raid0/zhouzb/tng_GrNr.npy', data)
        

