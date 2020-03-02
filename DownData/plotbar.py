import numpy as np
import h5py
import os
import illustris_python as il
import matplotlib.pyplot as plt

# def FilepathList(path, suffix='.hdf5'):
#     L = []
#     for files in os.listdir(path):
#         if os.path.splitext(files)[1] == '%s'%suffix:
#             L.append(int(os.path.splitext(files)[0]))
#     return L    

# ids = FilepathList('/Raid0/zhouzb/TNG_a2/small/', '.npy')

#Vector Normalized
def Norm(vect):
    '''Input a vector, normalize it in-place'''
    mol=(vect[0]**2+vect[1]**2+vect[2]**2)**0.5
    for i in range(3):
        vect[i]=vect[i]/mol

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

def rothalo(coor, Jz):
    '''
    set halo's angular momentum as z axis
    '''
    rot = RotMatrix(Norm(Jz))

    for i in range(len(coor)):
        coor[i] = np.dot(rot, coor[i])    

    return coor

def SelectDisk(snap_num):
    '''
    Input Snapshot number, like snap_num = 99 (z=0)
    Select disk galaxies, return haloID of them.
    '''
    #select halo stellar particles > 20000
    stellarlen = il.groupcat.loadSubhalos('/Raid1/Illustris/Illustris-1/', 99, 'SubhaloLenType')[:,4]

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

    mas_mask = stellarlen[haloID] > 20000
    haloID = haloID[mas_mask]
    return haloID


StellarHalfmassRads = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloHalfmassRadType')[:,4]
SubhaloPositions = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloPos')
Stellarlen = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100', 99, 'SubhaloLenType')[:,4]
#Count the number of barred galaxies
barred = 0

#Calculate each disk halo's A2 parameter
ids = SelectDisk(99)
for haloID in ids:
    coor = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', 99, haloID, partType=4, fields='Coordinates')
    mas = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', 99, haloID, partType=4, fields='Masses')
    vel = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', 99, haloID, partType=4, fields='Velocities')
    half_r = StellarHalfmassRads[haloID]
    halo_position = SubhaloPositions[haloID]

    sf_time = il.snapshot.loadSubhalo('/Raid1/Illustris/TNG-100', 99, haloID, partType=4, fields='GFM_StellarFormationTime')
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

    np.save('/Raid0/zhouzb/TNG_a2/snap_99/%d.npy'%haloID,np.vstack((A2_R,r_list[r_sort])))
    print('halo_%d datalen: %d'%(haloID,len(A2_R)))
    #Plot barred galaxies
    if A2_R.max() > 0.15:
        barred += 1

        plt.rcParams['figure.figsize'] = (21.6, 14.4)
        plt.scatter(coor[:,0], coor[:,1], marker='.')
        if Stellarlen[haloID] > 20000:
            plt.savefig('/Raid0/zhouzb/TNG_a2/barfig/big/%d.png'%haloID)
        else:
            plt.savefig('/Raid0/zhouzb/TNG_a2/barfig/small/%d.png'%haloID)
        plt.close()
        print('Halo %d is barred, figer saved.')

print('Disk galaxies number: %d'%len(ids))
print('Barred: %d'%len(barred))
print('Bar fraction: ')
print(len(barred) / len(ids))

