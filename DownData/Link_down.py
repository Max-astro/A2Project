import requests
import sys
import h5py
import numpy as np
import os

def get(path, params=None, savedir=None):
    # make HTTP GET request to path
    headers = {"api-key":"27d44ba55cd115b10f2dd9153589aff0"}
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
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/%haloID/sublink/simple.json"%haloID
    try:
        sublink = get(url, savedir='/home/sublink/')

    except:
        print(sys.exc_info()[0])
        return -1

    f = sublink
 
    #Find halo's Subfind ID with redshift(ie:SnapNum), and save the dict in '/Raid0/zhouzb/diskHalo_Sublink/'
    snap_num = np.array(f['SnapNum'])
    subfind_ID = np.array(f['SubfindID'])
    Progenitors_dict = {}
    for i in range(len(snap_num)):
        Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]

    f.close()
    return Progenitors_dict



'''
snap_91 z=0.1
snap_84 z=0.2
snap_78 z=0.3
snap_72 z=0.4
snap_67 z=0.5
snap_59 z=0.7
snap_50 z=1.0
snap_40 z=1.5
snap_33 z=2.0
'''


barred = np.load('F:/Linux/data/099fig/barredID.npy')
snap = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33]
errorHalo = []

for haloID in barred:
    Prog_dict = HaloProgenitors(haloID)
    if Prog_dict == -1:
        print('halo: %d Network ERROR, Try next'%haloID)
        errorHalo.append(haloID)
        continue
    else:    
        #Download stellar particles' information in all selected snapshot z
        for z in snap:
            print('Now download halo %d in snap_%d'%(haloID, z))
            try:
                subID = Prog_dict['%d'%z]
                cutoff_url = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime'%(z, subID)
                if os.path.isfile('F:/Linux/data/TNG/cutoff/disk_%d/cutout_%d.hdf5'%(z, subID)) == False:
                    get(cutoff_url, savedir='F:/Linux/data/TNG/cutoff/disk_%d/'%z)
            except:
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
    np.save('F:/Linux/data/TNG/errorID.npy', errorHalo)
        









