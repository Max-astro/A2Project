import requests
import sys
import os
import json
import numpy as np
import datetime

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
        ldir = 'g:/TNG/down/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = 'g:/TNG/down/il1_DiskMerTree/%d.json' % subhaloID
    with open(ldir) as f:
        data = json.load(f)
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

def getParticles(savedir, simu, snap, haloID, fields='?dm=Coordinates', isSubhalo=True):
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

def getHalo(savedir, simu, snap, haloID, fields='?dm=Coordinates'):
    baseUrl = 'http://www.tng-project.org/api/%s/snapshots/%d/halos/%d/cutout.hdf5' % (simu, snap, haloID)
    dataUrl = baseUrl + fields

    if not os.path.isdir(savedir):
        print('Savedir Invalid')
        return - 1

    savedir += '/snap_%d/'%snap
    if os.path.exists(savedir + 'cutout_%d.hdf5'%haloID):
        print("Already download")
        return 0
    print(dataUrl,'\n', savedir + 'cutout_%d.hdf5'%haloID)
    try:
        get(path=dataUrl, savedir=savedir)
    except:
        print(sys.exc_info()[0])
        return -1

    return 0


il1_snap = [120, 108, 95, 85, 75, 68, 64, 60]
tng_snap = [84, 67, 59, 50, 40, 33, 29, 25]


tng_diskID = np.load('g:/TNG/down/diskID_4WP.npy', allow_pickle=True)
il1_diskID = np.load('g:/Illustris-1/down/diskID_il1.npy', allow_pickle=True)

il1_GrNr = np.load('g:/Illustris-1/down/il1_GrNr.npy', allow_pickle=True).item()

done = 0
for haloID in il1_diskID:
    GrNr = il1_GrNr[haloID]
    for snap in il1_snap:
        gid = GrNr[snap]
        if gid == -1:
            continue
        if os.path.isfile('g:/Illustris-1/dark/snap_%d/cutout_%d.hdf5' % (snap, gid)):
            print("snap_%d halo_%d Already done."%(snap, gid))
            continue
        getHalo('g:/Illustris-1/dark/', 'Illustris-1', snap, gid, fields='?dm=Coordinates')

    done += 1    
    output = open('g:/Illustris-1/dark/log.txt','a+')
    
    print("Halo_%d finished. Done: %d\n" % (haloID, done))
    print(datetime.datetime.now().strftime('Time: %h.%d %T\n'))
    
    output.write("Halo_%d finished. Done: %d\n" % (haloID, done))
    output.write(datetime.datetime.now().strftime('Time: %h.%d %T\n'))
    output.close()

print('-------------------FINISHED!-------------------------')

