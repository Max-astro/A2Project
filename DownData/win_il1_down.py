import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
import requests
sys.path.append('F:\Linux')
import illustris_python as il


def get(path, params=None, savedir=None):
    # make HTTP GET request to path
    headers = {"api-key":"44140abd45569b745a4e3e7cf1284ce0"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        if savedir != None:
            filename = savedir + filename
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r


def HaloProgenitors(haloID):
    '''
    haloID is the subhalo's ID in snap_099
    return a dict = {'SnapNum' : SubfindID}
    '''
    url = 'http://www.tng-project.org/api/Illustris-1/snapshots/135/subhalos/%d'%haloID
    try:
        info = get(url)
        sublink = get(info['trees']['sublink_mpb'], savedir='F:/Linux/data/il1_cutoff/sublink')
    except:
        print(sys.exc_info()[0])
        return -1

    f = h5py.File(sublink, 'r')
 
    #Find halo's Subfind ID with redshift(ie:SnapNum), and save the dict in '/Raid0/zhouzb/diskHalo_Sublink/'
    snap_num = np.array(f['SnapNum'])
    subfind_ID = np.array(f['SubfindID'])
    Progenitors_dict = {}
    for i in range(len(snap_num)):
        Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]

    return Progenitors_dict

#Calculate subhalo's A2 parameter, return a list of A2 inside the half mass radius

'''
snap_135 z=0
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.31
snap_108 z=0.4
snap_103 z=0.5
snap_85 z=1.0
snap_74 z=1.53
snap_68 z=2.0
'''
def FilepathList(path, suffix='.hdf5'):
    L = []
    for files in os.listdir(path):
        if os.path.splitext(files)[1] == '%s'%suffix:
            L.append(int((os.path.splitext(files)[0])))
    return L

disk = FilepathList('F:/Linux/data/135_4WP_ALL/','.png')
snap = [127, 120, 113, 103, 108, 95, 85, 74, 68]
errorHalo = []

for haloID in disk:
    Prog_dict = HaloProgenitors(haloID)
    if Prog_dict == -1:
        print('halo: %d Network ERROR, Try next'%haloID)
        errorHalo.append(haloID)
        continue
    else:  
        np.save('F:/Linux/data/il1_cutoff/npy/halo_%d.npy'%haloID, Prog_dict)
        errtime = 0
        #Download stellar particles' information in all selected snapshot z
        for z in snap:
            print('Now downloading halo %d in snap_%d'%(haloID, z))
            try:    
                subID = Prog_dict['%d'%z]
            except:
                continue
            cutoff_url = 'http://www.tng-project.org/api/Illustris-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime'%(z, subID)
            try:
                data = get(cutoff_url, savedir='F:/Linux/data/il1_cutoff/disk_%d/'%z)
            except:
                errtime += 1
                print("halo %d in snap_%d Fail:"%(haloID, z), sys.exc_info()[0])
                print("You need to reload this halo.")
                errorHalo.append(haloID)
                break
            else:
                print('halo %d in snap_%d downloaded'%(haloID, z))
                         
        print('halo %d in all snapshot download Completed'%haloID)

if len(errorHalo) == 0:
    print('All done.')
else:
    print('%d halo download faild'%len(errorHalo))
    print("Error halo's ID were saved in '/Raid0/zhouzb/downError.log.npy'.")
    np.save('F:/Linux/data/il1_cutoff/downError.log.npy', errorHalo)
        




    







