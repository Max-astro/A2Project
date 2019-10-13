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

tng_snap = [53, 47, 45, 43, 41, 38, 36, 35, 34]
Redshift = [0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9]

for snap in tng_snap:
    done = 0
    for i in range(448):
        if os.path.exists('/Raid1/Illustris/TNG-100/Groupcat/fof_subhalo_tab_0%d.%d.hdf5'%(snap, i)):
            continue
        url = "http://www.tng-project.org/api/TNG100-1/files/groupcat-%d.%d.hdf5" % (snap, i)
        try:
            get(url, savedir='/Raid1/Illustris/TNG-100/Groupcat/')
            done += 1
            print('Snap_%d done: %d / 448' % (snap, done))
        except:
            print('Snap_%d part %d failed.' % (snap,done))
    print('Snap %d Finished.'%snap)

