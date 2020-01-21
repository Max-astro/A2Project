import requests
import sys
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

il1_snap = [99, 92, 89, 82, 80, 78, 76, 73, 71, 70, 69]
tng_snap = [63, 56, 53, 47, 45, 43, 41, 38, 36, 35, 34]
Redshift = [0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9]






base_url = "http://www.tng-project.org/api/TNG100-1/"
params = {'stars':'Masses,Coordinates,Velocities,GFM_StellarFormationTime'}
sim_metadata = get(base_url)
snap = [72]

for snapnum in snap:
    for i in range(sim_metadata['num_files_snapshot']):
        file_url = base_url + "files/snapshot-%d."%snapnum + str(i) + ".hdf5"
        if os.path.isfile('j:/TNG/test/snap_%d/snap_0%d.%d.hdf5'%(snapnum, snapnum, i)) == False:
            saved_filename = get(file_url, params, 'j:/TNG/test/snap_%d/'%snapnum)
            print(saved_filename)
