import numpy as np
import h5py
import json
import sys
import csv
import os
sys.path.append('F:\Linux')
import illustris_python as il

def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if simu == 'TNG':
        ldir = 'f:/Linux/localRUN/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = 'f:/Linux/localRUN/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])

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


def loadCutoff(simu, haloID, snapnum, partType, fields):
    res = {}
    if simu == 'TNG' or simu == 'tng':
        f = h5py.File('g:/Linux/il1_cutoff_dm/tng/snap_%d/cutout_%d.hdf5' % (snapnum, haloID))
    else:
        f = h5py.File('g:/Linux/il1_cutoff_dm/snap_%d/cutout_%d.hdf5' % (snapnum, haloID))
    for i, field in enumerate(fields):
        res[field] = np.array(f['PartType%d' % partType][field])
    return res
            

def subhaloA2_Cutoff(simu, haloID, snapnum):
    if haloID == -1:
        return - 1
    try:
        data = loadCutoff(simu, snapnum=snapnum, haloID=haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
        coor= data['Coordinates']
        mas = data['Masses']
        vel = data['Velocities']
        sf_time = data['GFM_StellarFormationTime']
    except:
        print(sys.exc_info()[0])
        print('halo %d in snap_%d load faild.'%(haloID,snapnum))
        print(' ')
        return -1


    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = il.func.loadSubhalos(simu, snapnum, 'SubhaloPos')[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - halo_position
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
        if coor[star, 1] == 0:
            continue
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



il1_snap = [99, 92, 89, 82, 80, 78, 76, 73, 71, 70, 69]
tng_snap = [63, 56, 53, 47, 45, 43, 41, 38, 36, 35, 34]
Redshift = [0.6, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.6, 1.7, 1.8, 1.9]


# barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
count = 0
for subhaloID in il1_barID:
    f = open('f:/Linux/localRUN/subA2_il1_99to69.csv', 'a', newline='')
    f_csv = csv.writer(f)
    A2 = []
    prog = LoadMergHist('il1', subhaloID)[0]
    for snapnum in il1_snap:
        print('Processing halo_%d in snap_%d' % (subhaloID, snapnum))
        try:
            haloID = prog[snapnum]
            a2tmp = subhaloA2_Cutoff('il1', haloID, snapnum)
        except:
            print(sys.exc_info()[0])
            A2.append(-1)
            continue
        A2.append(a2tmp)

    A2.insert(0, subhaloID)
    f_csv.writerow(A2)
    count += 1
    print('Halo %d done. %d / %d'%(subhaloID, count, len(il1_barID)))
    f.close()

