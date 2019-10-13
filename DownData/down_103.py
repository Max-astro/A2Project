import requests
import sys
import h5py
import numpy as np

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"0e689ed4aade028427f1f55aa3cd3488"}

def get(path, params=None):
    # make HTTP GET request to path
    headers = {"api-key":"0e689ed4aade028427f1f55aa3cd3488"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def finddisk(snapnum):
    '''Return a list of all Disk galaxies in this snapshot'''
    with h5py.File('/Raid0/zhouzb/stellar_circs.hdf5','r') as cir:
        return np.array(cir['Snapshot_%d'%snapnum]['SubfindID'])

disk = finddisk(103)


halourl = 'http://www.tng-project.org/api/Illustris-1/snapshots/103/subhalos/'
#Subhalo's Masses and Coordinates url
# MCurl = 'http://www.tng-project.org/api/Illustris-1/snapshots/103/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates'%haloID

star_disk = []
errorid = []
for haloID in disk:

    try:
        temp = get(halourl+str(haloID))
        if temp['len_stars'] >= 40000:
            MCurl = 'http://www.tng-project.org/api/Illustris-1/snapshots/103/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates'%haloID
            try:
                get(MCurl)
            except:
                print("%d Fail:"%haloID, sys.exc_info()[0])
                errorid.append(haloID)
                continue
            else:
                star_disk.append(haloID)
                print('halo %d downloaded'%haloID)
    except:
        continue


if len(errorid) == 0:
    print('All done. Total cutout number:')
    print(len(star_disk))

else:
    print('Downloaded number: %d'%len(star_disk))
    print('Faild number: %d'%len(errorid))
    print('Error halo ID:')
    for i in errorid:
        print(i)

