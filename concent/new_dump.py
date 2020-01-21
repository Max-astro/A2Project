import numpy as np
import h5py
import json
import sys
sys.path.append('f:/Linux')
import illustris_python as il
import matplotlib.pyplot as plt
from struct import *

sys.path.append("C:/Users/qq651/OneDrive/Codes/A2project/")
from plotTools.plot import *

  

def dumpCutOff(path, simu, snapnum, halos, redshift):
    if simu == 'TNG' or simu == 'tng':
        fdir = 'g:/TNG/dark/'
    else:
        fdir = 'g:/Illustris-1/dark/'

    PtclPos = path + '.Pos'
    halosOut = path + '.halo'

    ptclFILE = open(PtclPos, 'wb')
    haloFILE = open(halosOut, 'wb')


    nPtcl = 0
    ptclFILE.write(pack('if', nPtcl, redshift))
    haloFILE.write(pack('i', len(halos)))

    done = 0
    for haloID in halos:
        try:
            f = h5py.File(fdir + 'snap_%d/cutout_%d.hdf5' % (snapnum, haloID), 'r')
        except IOError:
            print("Halo_%d in snap_%d. File not exist."%(haloID, snapnum))
            # f.close()
            return
        
        data = np.array(f['PartType1']['Coordinates'])
        haloFILE.write(pack('i', len(data)))  #nhalo
        nPtcl += len(data)
        
        for i in data:
            for j in i:
                ptclFILE.write(pack('f', j))
        f.close()

#------------------------------Ptcl pos----------------------
        done += 1
        print("Done: %d / %d"%(done, len(halos)))

#------------------------------Groupcat-----------------------------------        
    gcfields = ['Group_R_Crit200', 'GroupPos']
    gcat = il.func.loadhalos(simu, snapnum, fields=gcfields)


    for fds in gcfields:
        if fds == 'GroupPos':
            for ii in gcat[fds][halos]:
                haloFILE.write(pack('3f', ii[0], ii[1], ii[2]))
        elif fds == 'Group_R_Crit200':
            for ii in gcat[fds][halos]:
                haloFILE.write(pack('f', ii))
#------------------------------Groupcat----------------------------------------

    if (nPtcl > 2147483647):
        print("Particle number over INT_MAX!")
    else:
        print('All done, ptcl: ',nPtcl)

    ptclFILE.seek(0, 0)
    ptclFILE.write(pack('i', nPtcl,))

    haloFILE.close()
    ptclFILE.close()
    print('Total particle: ', nPtcl)


def findGid(simu, snap, ids):
    res = []
    GID = il.func.loadSubhalos(simu, snap, 'SubhaloGrNr')
    for haloID in ids:
        try:
            prog = LoadMergHist(simu, haloID)[0]
            subID = prog[snap]
        except:
            continue

        res.append(GID[subID])
    return np.unique(res)


il1_snap = [120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [84, 67, 59, 50, 40, 33, 29, 25]
Redshift = [0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0]


        
    
#------------------------TNG--------------------------
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy', allow_pickle=True)
disk_GID = findGid('TNG', 84, tng_diskID)
dumpCutOff('g:/ptcl/tng_84', 'TNG', 84, disk_GID, 0.2)

#------------------------TNIllustris-1 snap135--------------------------

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy', allow_pickle=True)
GID = il.func.loadSubhalos('il1', 135, 'SubhaloGrNr')
diskGID = GID[disk]
diskGID = np.unique(diskGID)

dumpPtcl('/Raid0/zhouzb/haloProfileIO/il1_135', 'il1', 135, diskGID, 0)


