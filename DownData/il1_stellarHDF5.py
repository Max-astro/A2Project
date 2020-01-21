import requests
import sys
import h5py
import numpy as np
import os
import time

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

'''
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.3
snap_108 z=0.4
snap_103 z=0.5
snap_95 z=0.7
snap_85 z=1.0
snap_75 z=1.5
snap_68 z=2.0
'''


base_url = "http://www.tng-project.org/api/Illustris-1/"
params = {'stars':'Masses,Coordinates,Velocities,GFM_StellarFormationTime'}
sim_metadata = get(base_url)
snap = [127, 120, 113, 108, 103, 95, 85, 75, 68]
for snapnum in snap:
    print('Start download snap_%d'%snapnum)
    for i in range(sim_metadata['num_files_snapshot']):
        time_start=time.time()
        file_url = base_url + "files/snapshot-%d."%snapnum + str(i) + ".hdf5"
        if os.path.isfile('/Raid1/Illustris/Illustris-1/Snapshot/snap_%d/snap_0%d.%d.hdf5'%(snapnum,snapnum,i))==False:
            saved_filename = get(file_url, params, '/Raid1/Illustris/Illustris-1/Snapshot/snap_%d/' % snapnum)

            time_end=time.time()
            print('Part %d saved.' % i, 'time: %.3fs'%(time_end - time_start))

    print('Snap %d finished.\n' % snapnum)
