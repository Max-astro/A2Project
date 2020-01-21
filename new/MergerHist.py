import numpy as np
import sys
import os
import json
import time
import csv
import illustris_python as il

# import h5py
# import illustris_python as il

'''
This code can determine whether the galaxy bar was form by merger event or secular evolution

'''

#Calculate subhalo's A2 parameter, return a list of A2 inside the half mass radius
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


def rothalo(coor, Jz):
    '''
    set halo's angular momentum as z axis
    '''
    rot = RotMatrix(Jz / np.linalg.norm(Jz))

    for i in range(len(coor)):
        coor[i] = np.dot(rot, coor[i])    

    return coor

def subhaloA2(haloID, snapnum, simu='TNG'):
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates','Masses','Velocities','GFM_StellarFormationTime'])
    coor = np.array(data['Coordinates'])
    mas = np.array(data['Masses'])
    vel = np.array(data['Velocities'])
    sf_time = np.array(data['GFM_StellarFormationTime'])

    half_r = (il.func.loadSubhalos(simu,snapnum,'SubhaloHalfmassRadType'))[haloID,4]
    halo_position = coor[0] + 0.00000001

    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - halo_position
    r[r > 37500] -= 75000
    r[r < -37500] += 75000
    dis = np.linalg.norm(r, axis=1)
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    #Calculate angular momentum
    V = np.sum(vel * mas[:, np.newaxis], 0) / mas.sum()
    vel = vel - V
    Jz = np.sum(np.cross(coor, vel * mas[:, np.newaxis]), axis=0)

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
    A2_R = (a2_R ** 2 + b2_R ** 2)** 0.5 / a0
    return A2_R
    

def LoadMergHist(simu, subhaloID):
    '''
    return subhalo's main progenitor and merger history with snapshot
    '''
    if simu == 'TNG':
        ldir = '/Raid0/zhouzb/merg_data/tng_DiskMerTree/%d.json' % subhaloID
    else:
        ldir = '/Raid0/zhouzb/merg_data/il1_DiskMerTree/%d.json' % subhaloID
    
    with open(ldir) as f:
        data = json.load(f)
    
    Main = np.array(data['Main'])
    return dict(zip(Main[:, 0], Main[:, 1])), np.array(data['Mergers'])


def isSecular(progID, SnapList, A2list, Mergers):
    for i in range(7, -1, -1):
        if (A2list[i] > 0.15) & (A2list[i + 1] <0.15) & (A2list[i] != 1) & (progID[i] != -1):
            #If no merge event in A2 > 0.15 snapshot, then this halo bar is evolved by secular evolution
            subhalos = (Mergers[Mergers[:, 0] == SnapList[i]])[:,1]
            if subhalos.size == 0:
                break
            else:
                MergMass = il.func.loadSubhalos('TNG', SnapList[i], 'SubhaloMass')[subhalos]
                Mass = il.func.loadSubhalos('TNG', SnapList[i], 'SubhaloMass')[progID[i]]

                for m in MergMass:
                    if Mass > 500 * m:
                        continue
                    elif (Mass > 4 * m) & (Mass < 500 * m):
                        return [SnapList[i], 'MinorMerge']
                    else:
                        return [SnapList[i], 'Merge']
            break
    return [-1, 'Secular']


'''
Illustris-1 Snapshot-Redshift:
snap_127 z=0.1
snap_120 z=0.2
snap_113 z=0.3
snap_108 z=0.4
snap_103 z=0.5
snap_95 z=0.7
snap_85 z=1.0
snap_75 z=1.5
snap_68 z=2.0

SnapList = [135, 127, 120, 113, 108, 103, 95, 85, 75, 68]
RedShift = [0, 0.1 , 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
'''


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
SnapList = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33]
RedShift = [0, 0.1 , 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]

#Load barred subhaloID, determine the type of origin of bar.
tng_barredID = np.load('/Raid0/zhouzb/npydata_TNG/bigID.npy')

#Write csv file
cfile = open('/Raid0/zhouzb/TNG_MergSnap.csv','a+', newline='')
csv_write = csv.writer(cfile, dialect = 'excel')

barFromedSnap = []
secular = 0
for haloID in tng_barredID:
    time_start=time.time()
    Main, Mergers = LoadMergHist('TNG', haloID)
    A2list = []
    for snap in SnapList[1:]:
        try:
            if Main[snap] == -1:
                A2list.append(1)
            else:
                A2 = subhaloA2(Main[snap], snap)
        except KeyError:
            A2list.append(0)
        else:
            A2 = A2[int(len(A2) / 100):]
            A2list.append(A2.max())

    haloMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMass')[haloID]
    progID = []
    for i in SnapList:
        try:
            progID.append(Main[i])
        except KeyError:
            progID.append(-1)
    oriSnap = isSecular(progID, SnapList, A2list, Mergers)
    csv_write.writerow([haloID, oriSnap[0], oriSnap[1]])

    time_end = time.time()

        

