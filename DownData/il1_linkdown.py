import requests
import sys
import h5py
import numpy as np
import illustris_python as il


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


def SelectDisk(snap_num):
    '''
    Input Snapshot number, like snap_num = 99 (z=0)
    Select disk galaxies, return haloID of them.
    '''
    #select halo stellar particles > 40000
    halolen = il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',135,'SubhaloLenType')[:,4]

    with h5py.File('/Raid0/zhouzb/data/stellar_circs.hdf5','r') as cir:
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

def HaloProgenitors(haloID):
    '''
    haloID is the subhalo's ID in snap_099
    return a dict = {'SnapNum' : SubfindID}
    '''
    url = 'http://www.tng-project.org/api/Illustris-1/snapshots/135/subhalos/%d'%haloID
    try:
        info = get(url)
        sublink = get(info['trees']['sublink_mpb'], savedir='/Raid0/zhouzb/il1_cutoff/sub_hdf5/')
    except:
        print(sys.exc_info()[0])
        return -1

    f = h5py.File(sublink, 'r')
 
    #Find halo's Subfind ID with redshift(ie:SnapNum), and save the dict in '/Raid0/zhouzb/diskHalo_Sublink/'
    snap_num = np.array(f['SnapNum'])
    subfind_ID = np.array(f['SubfindID'])
    Progenitors_dict = {}
    for i in range(len(snap_num)):
        Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]

    return Progenitors_dict

#Calculate subhalo's A2 parameter, return a list of A2 inside the half mass radius
def subhaloA2(haloID, snapnum):
    with h5py.File('/Raid0/zhouzb/il1_cutoff/disk_%d/cutout_%d.hdf5'%(snapnum, haloID),'r') as f:
        coor = np.array(f['PartType4']['Coordinates'])
        mas = np.array(f['PartType4']['Masses'])
        vel = np.array(f['PartType4']['Velocities'])
        sf_time = np.array(cut['PartType4']['GFM_StellarFormationTime'])

    half_r = (il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',snapnum,'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',snapnum,'SubhaloPos')[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles


    r = coor - halo_position
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    #Calculate angular momentum
    for i in range(len(vel)):
        vel[i] = vel[i] * mas[i]
    J = np.cross(coor, vel)
    Jz = J.sum(axis = 0)

    #Set halo's angular momentum Jz as z axis
    coor = rothalo(coor, Jz)

    #r_list is a list about the distance between stellar particles and halo_position
    #r_sort a index about paritcles' distance
    r_list = ((coor[:,:2]**2).sum(1))**0.5
    r_sort = r_list.argsort()

    #a0 = Sigma(M[i])
    a0 = np.zeros(len(r_sort))
    ptr = 0
    for i in r_sort:
        if ptr == 0:
            a0[ptr] = mas[i]
        else:
            a0[ptr] = mas[i] + a0[ptr-1]
        ptr += 1

    #Position angle: Theta = arctan(y/x)
    #Creat a list of a_m(R)
    a2_R = np.zeros(len(r_sort))
    b2_R = np.zeros(len(r_sort))
    list_i = 0
    for star in r_sort:
        Theta = np.arctan(coor[star,0] / coor[star,1])
        #a_i = M[numb] * cos( m * Theta[numb] )  , m=2
        a_i = mas[star] * np.cos(2*Theta)
        #b_i = M[numb] * sin( m * Theta[numb] ) , m=2
        b_i = mas[star] * np.sin(2*Theta) 

        #a2_R[i] = a_i.sum(:i)
        if list_i > 0:
            a2_R[list_i] = a_i + a2_R[list_i - 1]
            b2_R[list_i] = b_i + b2_R[list_i - 1]
        else:
            a2_R[list_i] = a_i
            b2_R[list_i] = b_i

        list_i += 1
        
    #It's a list about A2 parameter, sorted by distance between halo and center
    A2_R = (a2_R**2 + b2_R**2)**0.5 / a0
    return A2_R


def plot_A2(Prog_dict, haloID):
    A2list = []
    halolen = []
    snaplist = []

    plt.figure()
    plt.xlabel('Red Shift')
    plt.ylabel('A2')

    #Calculate A2 parameter 
    tmp = subhaloA2(haloID, 0)
    A2 = np.array(tmp[0])
    A2 = A2[int(len(A2) / 100):]
    if A2.max() < 0.15:
        continue

    A2list.append(A2.max())
    snaplist.append(0)


    #Load halo's haloID with redshift
    SnapDict = Prog_dict
    #Load data with redshift
    for snap_num in snapshot:
        try:
            haloID = SnapDict['%d'%snap_num]
            tmp = subhaloA2(haloID, snap_num))
            snaplist.append(Redshift['%d'%snap_num])
        except:
            print('Halo %d no found in snapshot %d'%(diskID, snap_num))
            continue

        A2 = np.array(tmp)
        A2 = A2[int(len(A2) / 100):]
        A2list.append(A2.max())

        #if particle number less than 40000, use carmine and yellow
        halolen.append(il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1',snap_num,'SubhaloMassType')[haloID,4])



    #Normal Color: blue = StellarMass, red = GasFraction
    #Small halo color: cyan = StellarMass, carmine = GasFraction
    #if particle number less than 40000, use carmine and yellow
    mask = np.array(halolen) < 40000

    #Plot StellarMass in ax1 with black point, and GasFraction in ax2 with blue point
    # red = ax1.plot(x_point, mass_point, 'ob', color = 'r', label = 'StellarMass')
    plt.scatter(snaplist, A2list, marker='o', color='b')
    plt.plot(snaplist, A2list)
    if mask.any() != False:
        A2list = np.array(A2list[1:])
        snaplist = np.array(snaplist[1:])


        # carmine = ax1.plot(x_point[mask], mass_point[mask], 'ob', color = 'c', label = 'StellarMass, P<40000')
        yellow = plt.plot(snaplist[mask], A2list[mask], 'ob', color = 'y', label = 'Stellar particles <40000')
        plt.legend()
    plt.xlim(-0.1, 1.1)
    plt.ylim(0, 0.6)
    plt.savefig('/Raid0/zhouzb/fig_TNG/z-a2/%d.png'%diskID)
    print('fig_%d saved'%diskID)    


'''
snap_135 z=0
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.31
snap_108 z=0.4
snap_103 z=0.5
snap_85 z=1.0
snap_75 z=1.53
snap_68 z=2.0
'''


disk = SelectDisk(135)
snap = [135, 127, 120, 113, 103, 108, 95, 85, 74, 68]
errorHalo = []

for haloID in disk[:2]:
    Prog_dict = HaloProgenitors(haloID)
    if Prog_dict == -1:
        print('halo: %d Network ERROR, Try next'%haloID)
        errorHalo.append(haloID)
        continue
    else:  
        np.save('/Raid0/zhouzb/il1_cutoff/sublink/halo_%d.npy'%haloID, Prog_dict)
        errtime = 0
        #Download stellar particles' information in all selected snapshot z
        for z in snap:
            print('Now downloading halo %d in snap_%d'%(haloID, z))
            subID = Prog_dict['%d'%z]
            cutoff_url = 'http://www.tng-project.org/api/Illustris-1/snapshots/%d/subhalos/%d/cutout.hdf5?stars=Masses,Coordinates,Velocities,GFM_StellarFormationTime'%(z, subID)
            try:
                data = get(cutoff_url, savedir='/Raid0/zhouzb/il1_cutoff/disk_%d/'%z)
            except:
                errtime += 1
                print("halo %d in snap_%d Fail:"%(haloID, z), sys.exc_info()[0])
                print("You need to reload this halo.")
                errorHalo.append(haloID)
                break
            else:
                print('halo %d in snap_%d downloaded'%(haloID, z))

            #If no error, calculate A2 and plot picture
        # if errtime == 0:
        #     plot_A2(Prog_dict, haloID)



                         
        print('halo %d in all snapshot download Completed'%haloID)

if len(errorHalo) == 0:
    print('All done.')
else:
    print('%d halo download faild'%len(errorHalo))
    print("Error halo's ID were saved in '/Raid0/zhouzb/downError.log.npy'.")
    np.save('/Raid0/zhouzb/downError.log.npy', errorHalo)
        




    







