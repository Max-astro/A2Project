import requests
import sys
import h5py
import numpy as np
import os
import json

def get(path, params=None, savedir=None, filename=''):
    # make HTTP GET request to path
    headers = {"api-key":"6be79adf6f05af759ee1ada46249f9b4"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        if savedir != None:
            filename = savedir + filename + '.json'
        with open(filename, 'w') as f:
            json.dump(r.json(), f)
        return filename # return the filename string

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        if savedir != None:
            filename = savedir + filename
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def HaloProgenitors(haloID, simu):
    '''
    haloID is the subhalo's ID in snap_099
    return a dict = {'SnapNum' : SubfindID}
    '''
    if simu == 'TNG' or simu == 'tng':
        url = 'http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/%d/sublink/simple.json' % haloID
        savedir = '/Users/max/OneDrive/Documents/tng_DiskMerTree/'
    else:
        url = 'http://www.tng-project.org/api/Illustris-1/snapshots/135/subhalos/%d/sublink/simple.json' % haloID
        savedir = '/Users/max/OneDrive/Documents/il1_DiskMerTree'
    try:
        sublink = get(url,savedir=savedir,filename=str(haloID))
    except requests.exceptions.ConnectionError:
        print(sys.exc_info()[0])
        return - 1
    print("halo ", haloID, "finished.")
    return sublink

diskID = np.load('/Users/max/OneDrive/Documents/diskID_il1.npy')
l=0
d=0
for haloID in diskID:
    if os.path.exists('/Users/max/OneDrive/Documents/il1_DiskMerTree/%d.json' % haloID):
        d += 1
    else:
        HaloProgenitors(haloID, 'il1')
        l += 1
    print('done: ', d, ', err finish: ', l, ', all: ', l+d)