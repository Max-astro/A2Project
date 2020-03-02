import requests
import sys
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

base_url = "http://www.tng-project.org/api/TNG100-1/"
params = {'stars':'Masses,Coordinates,Velocities,GFM_StellarFormationTime'}
sim_metadata = get(base_url)
snap = [67, 59, 50]
for snapnum in snap:
    for i in range(sim_metadata['num_files_snapshot']):
        file_url = base_url + "files/snapshot-%d."%snapnum + str(i) + ".hdf5"
        if os.path.isfile('/Raid1/Illustris/TNG-100/Snapshot/snap_%d/snap_0%d.%d.hdf5'%(snapnum, snapnum, i)) == False:
            saved_filename = get(file_url, params, '/Raid1/Illustris/TNG-100/Snapshot/snap_%d/'%snapnum)
            print(saved_filename)
