import requests
import sys
import os
import json
import numpy as np

def us_get(path, params=None):
    # make HTTP GET request to path
    headers = {"api-key":"6be79adf6f05af759ee1ada46249f9b4"}
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    return r

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

il1_snap = [120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [84, 67, 59, 50, 40, 33, 29, 25]

# rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0])


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
        subUrl = 'http://www.tng-project.org/api/Illustris-1/snapshots/%d/subhalos/%d/' % (snap, subID)
        urls = get(subUrl)
        dataUrl = urls['cutouts']['parent_halo'] + '?dm=Coordinates,Velocities'
        try:
            get(dataUrl, savedir='/home/snap_%d/' % snap)
            print("sub_%d done" % subID)
        except requests.exceptions.ConnectionError:
            print('halo_%d in snap_%d ERROR' % (haloID, snap))
    done += 1
    print('Halo %d Finished.' % haloID, ' Done: %d / %d' % (done, len(il1_barID)))


def getParticles(savedir, simu, snap, haloID, fields='?dm=Coordinates,Velocities', isSubhalo=True):
    baseUrl = 'http://www.tng-project.org/api/%s/snapshots/%d/subhalos/%d/' % (simu, snap, subID)
    if not isSubhalo:
        baseUrl = us_get(baseUrl)['cutouts']['parent_halo']
    dataUrl = baseUrl + fields

    try:
        get(dataUrl, savedir=savedir + str(snap))
            print("sub_%d done" % subID)
    except requests.exceptions.ConnectionError:
        print('ERROR')
        return -1
    
def getHaloPtcl(savedir, simu, snapList, haloID, fields='?dm=Coordinates,Velocities', isSubhalo=True):
    prog = LoadMergHist(simu, haloID)[0]
    for snap in snapList:
        try:
            subID = prog[snap]
            if subID != -1:
                getParticles(savedir, simu, snap, subID, fields, isSubhalo)
        except:
            continue
