import requests
import sys
import os
import json
import numpy as np

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

def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if (simu == 'TNG') or (simu == 'tng'):
        ldir = 'f:/Linux/localRUN/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = 'f:/Linux/localRUN/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])


rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0])
il1_snap = [135, 120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [99, 84, 67, 59, 50, 40, 33, 29, 25]

# for i in tng_snap:
#     os.system("mkdir snap_%d"%i)

tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

done = 0
for haloID in tng_diskID:
    prog = LoadMergHist('TNG', haloID)[0]
    for snap in tng_snap:
        try:
            subID = prog[snap]
        except:
            continue

        MCurl = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/%d/cutout.hdf5?bhs=all' % (snap, subID)
        try:
            get(MCurl, savedir='f:/Linux/TNG_cutoff/bhs/snap_%d/' % snap)
            print("sub_%d done" % subID)
        except requests.exceptions.ConnectionError:
            print('halo_%d in snap_%d ERROR' % (haloID, snap))
    done += 1
    print('Halo %d Finished.' % haloID, ' Done: %d / %d' % (done, len(tng_diskID)))