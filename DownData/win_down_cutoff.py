import requests
import sys
import h5py
import numpy as np

def get(path, params=None, savedir=None):
    # make HTTP GET request to path
    headers = {"api-key":"44140abd45569b745a4e3e7cf1284ce0"}
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

def HaloProgenitors(haloID):
    '''
    haloID is the subhalo's ID in snap_099
    return a dict = {'SnapNum' : SubfindID}
    '''
    url = 'http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/%d'%haloID
    try:
        info = get(url)
        sublink = get(info['trees']['sublink_mpb'])
    except:
        print(sys.exc_info()[0])
        return -1

    with h5py.File(sublink, 'r') as f:
        #Find halo's Subfind ID with redshift(ie:SnapNum), and save the dict in '/Raid0/zhouzb/diskHalo_Sublink/'
        snap_num = np.array(f['SnapNum'])
        subfind_ID = np.array(f['SubfindID'])

    Progenitors_dict = {}
    for i in range(len(snap_num)):
        Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]

    return Progenitors_dict


def SelectDisk(snap_num):
    '''
    Input Snapshot number, like snap_num = 99 (z=0)
    Select disk galaxies, return haloID of them.
    '''
    #select halo stellar particles > 40000
    with h5py.File('/Raid0/zhouzb/TNGdata/offsets_0%d.hdf5'%snap_num,'r') as offset:
        haloSBT = (np.array(offset['Subhalo']['SnapByType']))[:,4]
    #Total halo number, it also be an index of haloID

    halolen = haloSBT[1:]
    halolen = np.append(halolen, halolen.max()) - haloSBT
    

    with h5py.File('/Raid0/zhouzb/TNGdata/stellar_circs.hdf5','r') as cir:
        haloID = np.array(cir['Snapshot_%d'%snap_num]['SubfindID'])
        cir07frac = np.array(cir['Snapshot_%d'%snap_num]['CircAbove07Frac'])
        MassTensor = np.array(cir['Snapshot_%d'%snap_num]['MassTensorEigenVals'])
    #circularity parameter ϵ > 0.2
    cir_mask = cir07frac > 0.2
    
    #flatness of the galaxy is defined as the ratio M1=(M2M3)**0.5 , disk galaxy's flatness < 0.7
    MassTensor = MassTensor[cir_mask]
    haloID = haloID[cir_mask]

    flat = MassTensor[:,0]/(MassTensor[:,1]*MassTensor[:,2])**0.5
    flat_mask = flat < 0.7

    haloID = haloID[flat_mask]

    mas_mask = halolen[haloID] >= 40000
    haloID = haloID[mas_mask]
    return haloID

'''
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


barred = np.load('F:/Linux/data/099fig/barredID.npy')
snap = [67, 59, 50, 40, 33]
errorHalo = []


for z in snap:
    print('Downloading snap_%d'%z)
    disk = SelectDisk(z)

    for haloID in disk:
        cutoff_url = 'http://www.tng-project.org/api/TNG100-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime'%(z, haloID)
        try:
            get(cutoff_url, savedir='F:\\Linux\\data\\TNG\\snap_%d\\'%z)
        except:
            print("halo %d in snap_%d Fail:"%(haloID, z), sys.exc_info()[0])

            errorHalo.append(haloID)
            break
        else:
            print('halo %d in snap_%d downloaded'%(haloID, z))

    print('halo %d in all snapshot download Completed'%haloID)

if len(errorHalo) == 0:
    print('All done.')
else:
    print('%d halo download faild'%len(errorHalo))
    print("Error halo's ID were saved in 'F:\\Linux\\data\\TNG\\downError.log.npy'.")
    np.save('F:\\Linux\\data\\TNG\\downError.log.npy', errorHalo)
        
