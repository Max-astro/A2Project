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
    if simu == 'TNG' or simu == 'TNG100-1':
        ldir = '/Raid0/zhouzb/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = '/Raid0/zhouzb/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)

    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

def getParticles(savedir, simu, snap, haloID, fields='?dm=Coordinates,Velocities', isSubhalo=True):
    baseUrl = 'http://www.tng-project.org/api/%s/snapshots/%d/subhalos/%d/' % (simu, snap, haloID)
    if not isSubhalo:
        baseUrl = us_get(baseUrl)['cutouts']['parent_halo']
    subID = int(baseUrl.split('/')[-2])
    dataUrl = baseUrl + fields
    
    if not os.path.isdir(savedir):
        print('Savedir Invalid')
        return - 1
        
    savedir += '/snap_%d/'%snap
    if os.path.exists(savedir + 'cutout_%d.hdf5'%subID):
        print("Already download")
        return 0         
    print(dataUrl, savedir + 'cutout_%d.hdf5'%subID)
    try:
        get(path=dataUrl, savedir=savedir)
    except:
        print(sys.ext_info()[0])
        return -1

    return 0

def getHaloPtcl(savedir, simu, snapList, haloID, fields='?dm=Coordinates,Velocities', isSubhalo=True):
    prog = LoadMergHist(simu, haloID)[0]

    res = 0
    for snap in snapList:
        try:
            subID = prog[snap]
            if subID != -1:
                res = getParticles(savedir, simu, snap, subID, fields, isSubhalo)
        except:
            print(sys.exc_info()[0])
            res = -1
            continue
    if res == 0:
        print("%s halo: %d done"%(simu, haloID))
    else:
        print("%s halo: %d Download faild"%(simu, haloID))




il1_snap = [120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [84, 67, 59, 50, 40, 33, 29, 25]

# rs = np.array([0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0])


tng_diskID = np.load('/Raid0/zhouzb/rundata/diskID_4WP.npy')

done = 0
for haloID in tng_diskID:
    getHaloPtcl('/Raid0/zhouzb/cut_TNG/dark/', 'TNG100-1', tng_snap, haloID, isSubhalo=False)
    done += 1
    print('Halo done: ',done)
    print(' ')
  
 
"http://www.tng-project.org/api/TNG100-1/snapshots/84/halos/1/cutout.hdf5?dm=Coordinates,Velocities"