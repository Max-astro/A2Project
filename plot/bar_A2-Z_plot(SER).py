import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys
import requests
import illustris_python as il

def get(path, params=None, savedir=None):
    # make HTTP GET request to path
    headers = {"api-key":" 27d44ba55cd115b10f2dd9153589aff0"}
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

def rothalo(coor, Jz):
    '''
    set halo's angular momentum as z axis
    '''
    rot = RotMatrix(Norm(Jz))

    for i in range(len(coor)):
        coor[i] = np.dot(rot, coor[i])    

    return coor

#Vector Normalized
def Norm(vect):
    '''Input a vector, normalize it in-place'''
    mol=np.linalg.norm(vect)
    vect /= mol
    return vect

#Find the rotation matrix when input vector rotated to z_axis
def RotMatrix(vect):
    '''Input a vector, return a (3,3) rotation matrix '''
    x=vect[0]
    y=vect[1]
    z=vect[2]

    sinA= y/(z**2+y**2)**0.5
    cosA= z/(z**2+y**2)**0.5
    sinB= x/(x**2+(y*sinA+z*cosA)**2)**0.5
    cosB= (y*sinA+z*cosA)/(x**2+(y*sinA+z*cosA)**2)**0.5

    return np.array([[cosB, -sinA*sinB, -cosA*sinB], [0, cosA, -sinA], [sinB, sinA*cosB, cosA*cosB]])

#Calculate subhalo's A2 parameter, return a list of A2 inside the half mass radius
def subhaloA2(simu, haloID, snapnum):
    if snapnum >40:
        try:
            coor = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields='Coordinates')
            mas = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields='Masses')
            vel = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields='Velocities')
            sf_time = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields='GFM_StellarFormationTime')
        except:
            print(sys.exc_info()[0])
            print('halo %d in snap_%d load faild.'%(haloID,snapnum))
            print(' ')
    else:
        try:
            with h5py.File('/Raid0/zhouzb/cut_TNG/disk/cutout_%d.hdf5','r') as f:
                coor = np.array(f['PartType4']['Coordinates'])
                mas = np.array(f['PartType4']['Masses'])
                vel = np.array(f['PartType4']['Velocities'])
                sf_time = np.array(f['PartType4']['GFM_StellarFormationTime'])
        except:
            print(sys.exc_info()[0])
            print('halo %d in snap_%d load faild, data lost.'%(haloID,snapnum))
            print(' ')

    half_r = (il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snapnum,'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snapnum,'SubhaloPos')[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles


    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    coor[coor > 37500] -= 75000
    coor[coor < -37500] += 75000
    #Calculate angular momentum
    V = np.sum(vel * mas[:, np.newaxis], 0) / mas.sum()
    vel = vel - V
    L = np.sum(np.cross(coor, vel * mas[:, np.newaxis]), axis=0)

    #Set halo's angular momentum Jz as z axis
    coor = rothalo(coor, L)

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
    A2_R = (a2_R ** 2 + b2_R ** 2)** 0.5 / a0
    A2_R = A2_R[int(len(A2_R) / 100):]
    return A2_R.max()


def Progenitors_dictionary(haloID):
    with h5py.File('/Raid0/zhouzb/sublink/sublink_mpb_%d.hdf5'%haloID,'r') as f:
        snap_num = np.array(f['SnapNum'])
        subfind_ID = np.array(f['SubfindID'])

        Progenitors_dict = {}
        for i in range(len(snap_num)):
            Progenitors_dict['%d'%snap_num[i]] = subfind_ID[i]
    return Progenitors_dict


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

snapshot = [91, 84, 78, 72, 67, 59, 50, 40, 33]

ToList = {'99': 0, '91' : 1, '84' : 2, '78' : 3, '72' : 4, '67' : 5, '59' : 6, '50' : 7, '40' : 8, '33' : 9}
Redshift = {'99': 0, '91' : 0.1, '84' : 0.2, '78' : 0.3, '72' : 0.4, '67' : 0.5, '59' : 0.7, '50' : 1.0, '40' : 1.5, '33' : 2.0}
#Redshift = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
'''
GasFraction / StellarMass[0] = snap_99
GasFraction / StellarMass[1] = 91,
GasFraction / StellarMass[2] = 84,
GasFraction / StellarMass[3] = 78,
GasFraction / StellarMass[4] = 72,
GasFraction / StellarMass[5] = 67,
GasFraction / StellarMass[6] = 59,
GasFraction / StellarMass[7] = 50,
GasFraction / StellarMass[8] = 40,
GasFraction / StellarMass[9] = 33
'''


ids = np.load('/Raid0/zhouzb/TNG/barredID.npy')

#Painting
for diskID in ids:
    A2list = []
    halolen = []
    snaplist = []

    plt.figure()
    plt.xlabel('Red Shift')
    plt.ylabel('A2')

    #Load A2 parameter inside the stellar half mass radius*2, and pick A2.max() as x point
    tmp = subhaloA2(diskID, 99)
    A2 = np.array(tmp)
    A2 = A2[int(len(A2) / 100):]
    if A2.max() < 0.15:
        continue

    A2list.append(A2.max())
    snaplist.append(0)


    #Load halo's haloID with redshift
    SnapDict = Progenitors_dictionary(diskID)
    #Load data with redshift
    for snap_num in snapshot:
        try:
            haloID = SnapDict['%d'%snap_num]
        except:
            print('Halo %d progeniter in snapshot %d no found.'%(diskID, snap_num))
            continue
        
        try:
            tmp = subhaloA2(haloID, snap_num)
        except:
            print('halo %d in snap_%d error'%(haloID, snap_num),sys.exc_info()[0])
            print('')
            continue
        
        snaplist.append(Redshift['%d'%snap_num])

        A2 = np.array(tmp)
        A2 = A2[int(len(A2) / 100):]
        A2list.append(A2.max())

        #if particle number less than 40000, use carmine and yellow
        halolen.append((il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', snap_num, 'SubhaloLenType')[:,4])[haloID])



    #Normal Color: blue = StellarMass, red = GasFraction
    #Small halo color: cyan = StellarMass, carmine = GasFraction
    #if particle number less than 40000, use carmine and yellow
    mask = np.array(halolen) < 20000

    #Plot StellarMass in ax1 with black point, and GasFraction in ax2 with blue point
    # red = ax1.plot(x_point, mass_point, 'ob', color = 'r', label = 'StellarMass')
    plt.scatter(snaplist, A2list, marker='o', color='b')
    plt.plot(snaplist, A2list)
    if mask.any() != False:
        A2list = np.array(A2list[1:])
        snaplist = np.array(snaplist[1:])


        # carmine = ax1.plot(x_point[mask], mass_point[mask], 'ob', color = 'c', label = 'StellarMass, P<40000')
        yellow = plt.plot(snaplist[mask], A2list[mask], 'ob', color = 'y', label = 'Stellar particles <20000')
        plt.legend()
    
    #point the first time a2<0.15
    droppoint = 0
    for t in range(len(A2list)):
        if A2list[t] < 0.15:
            droppoint = t
            break
    if droppoint != 0:
        plt.scatter(snaplist[droppoint], A2list[droppoint], marker='o', color='r')
        plt.annotate('Drop point', xy=(snaplist[droppoint], A2list[droppoint]), xycoords='data', xytext=(+8, +18), textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))


    
    plt.xlim(-0.1, 2.1)
    plt.ylim(0, 0.8)
    plt.title('halo %d'%diskID)
    plt.savefig('/Raid0/zhouzb/fig_TNG/new/%d.png'%diskID)
    print('fig_%d saved'%diskID)

