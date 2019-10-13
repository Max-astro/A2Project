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
    if simu == 'TNG':
        ldir = '/home/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = '/home/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

# il1_snap = [99, 92, 89, 82, 80, 78, 76, 73, 71, 70, 69]
# tng_snap = [63, 56, 53, 47, 45, 43, 41, 38, 36, 35, 34]
tng_snap = [25, 29]
il1_snap = [60, 64]
# Redshift = [0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9]


tng_barID = np.load('/home/barredID_4WP_TNG.npy')
il1_barID = np.load('/home/barredID_il1.npy')

done = 0
for haloID in il1_barID:
    prog = LoadMergHist('il1', haloID)[0]
    for snap in il1_snap:
        try:
            subID = prog[snap]
        except:
            continue

        MCurl = 'http://www.tng-project.org/api/Illustris-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime' % (snap, subID)
        try:
            get(MCurl, savedir='/home/snap_%d/' % snap)
            print("sub_%d done" % subID)
        except requests.exceptions.ConnectionError:
            print('halo_%d in snap_%d ERROR' % (haloID, snap))
    done += 1
    print('Halo %d Finished.' % haloID, ' Done: %d / %d' % (done, len(il1_barID)))


