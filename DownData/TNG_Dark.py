import requests
import sys
import numpy as np

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
TNG Snapshot-Redshift:
snap_91 z=0.1
snap_84 z=0.2
snap_78 z=0.3
snap_72 z=0.4
snap_67 z=0.5
snap_59 z=0.7
snap_50 z=1.0
snap_40 z=1.5
snap_33 z=2.0
'''

SnapList = [91, 84, 78, 72, 67, 59, 50, 40, 33]
RedShift = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

errorHalo = []
diskID = np.load('f:/Linux/localRUN/diskID_TNG.npy')
for haloID in diskID[611:]:
    url = "http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/%d/sublink/simple.json" % haloID
    Prog = np.array(get(url)['Main'])
    for snap in SnapList:
        try:
            subID = (Prog[Prog[:,0] == snap])[0][1]
        except:
            continue

        print('Start downloading halo_%d in snap_%d'%(haloID, snap))
        cutoff_url = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/%d/cutout.hdf5?dm=Coordinates' % (snap, subID)
        try:
            get(cutoff_url, savedir='f:/Linux/TNG_cutoff/snap_%d/' % snap)
        except:
            errorHalo.append(haloID)
            print('halo_%d in snap_%d ERROR'%(haloID, snap))
            continue
    print('halo_%d finished' % haloID)

print('done, error number: %d' % len(errorHalo))
np.save('f:/Linux/localRUN/errorID.npy',errorHalo)
