import requests
import sys
import h5py
import numpy as np

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"44140abd45569b745a4e3e7cf1284ce0"}

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


def SelectDisk(snap_num):
    '''
    Input Snapshot number, like snap_num = 99 (z=0)
    Select disk galaxies, return haloID of them.
    '''
    #select halo stellar particles > 40000
    with h5py.File('/Raid0/zhouzb/TNGdata/offsets_0%d.hdf5'%snap_num,'r') as offset:
        haloSBT = (np.array(offset['Subhalo']['SnapByType']))[:,4]
    #Total halo number, it also be an index of haloID

    halolen = haloSBT[1:]
    halolen = np.append(halolen, halolen.max()) - haloSBT
    

    with h5py.File('/Raid0/zhouzb/TNGdata/stellar_circs.hdf5','r') as cir:
        haloID = np.array(cir['Snapshot_%d'%snap_num]['SubfindID'])
        cir07frac = np.array(cir['Snapshot_%d'%snap_num]['CircAbove07Frac'])
        MassTensor = np.array(cir['Snapshot_%d'%snap_num]['MassTensorEigenVals'])
    #circularity parameter ϵ > 0.2
    cir_mask = cir07frac > 0.2
    
    #flatness of the galaxy is defined as the ratio M1=(M2M3)**0.5 , disk galaxy's flatness < 0.7
    MassTensor = MassTensor[cir_mask]
    haloID = haloID[cir_mask]

    flat = MassTensor[:,0]/(MassTensor[:,1]*MassTensor[:,2])**0.5
    flat_mask = flat < 0.7

    haloID = haloID[flat_mask]

    mas_mask = halolen[haloID] >= 40000
    haloID = haloID[mas_mask]
    return haloID

snapshot = [72, 67, 59, 50]

for snap_num in snapshot:
    print('Downloading snapshot %d'%snap_num)
    disk = SelectDisk(snap_num)
    halourl = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/'%snap_num
    #Subhalo's Masses and Coordinates url
    # MCurl = 'http://www.tng-project.org/api/Illustris-1/snapshots/103/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates'%haloID

    star_disk = []
    errorid = []
    for haloID in disk:
        MCurl = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime'%(snap_num, haloID)
        try:
            get(MCurl, savedir='/Raid0/zhouzb/cut_TNG/snap/snap_%d/'%snap_num)
        except:
            print("%d Fail:"%haloID, sys.exc_info()[0])
            errorid.append(haloID)
            continue
        else:
            star_disk.append(haloID)
            print('halo %d downloaded'%haloID)



    if len(errorid) == 0:
        print('All done. Total cutout number:')
        print(len(star_disk))

    else:
        print('Downloaded number: %d'%len(star_disk))
        print('Faild number: %d'%len(errorid))
        print('Error halo ID saved in: /Raid0/zhouzb/cut_TNG/snap/')
        np.save('/Raid0/zhouzb/cut_TNG/snap/snap_%d.error.npy',errorid)

