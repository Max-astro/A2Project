import numpy as np
import h5py
import illustris_python as il

basePath = '/Raid1/Illustris/TNG-100'

# def groups(snap, haloID, fields, data):
    # '''
    # parameter:   
    # fields : ['Group', 'Header', 'Offsets', 'Subhalo']
    # data: halo's information that you want to extract
    # '''
    # with h5py.File('/Raid1/Illustris/TNG-100/Groupcat/fof_subhalo_tab_%03d.0.hdf5'%snap,'r') as f:
    #     info = np.array(f['%s'%fields]['%s'%data])
    # for n in range(1,8):
    #     with h5py.File('/Raid1/Illustris/TNG-100/Groupcat/fof_subhalo_tab_%03d.%d.hdf5'%(snap, n),'r') as f:
    #         info = np.concatenate((info, np.array(f['%s'%fields]['%s'%data])))
    # if haloID == -1:
    #     return info
    # else:
    #     return info[haloID]


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
    #select halo stellar particles > 40000
    with h5py.File('/Raid0/zhouzb/TNGdata/offsets_0%d.hdf5'%snap_num,'r') as offset:
        haloSBT = (np.array(offset['Subhalo']['SnapByType']))[:,4]
    #Total halo number, it also be an index of haloID
    ids = np.arange(len(haloSBT))

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
Now calculate all of this snapshot's disk halos A2parameter
Each disk halo will be return two list : 1. A2 ; 2. The corresponding radius. They are sorted by distance between stellar particle and halo center
'''
snap_num = 50
ids = SelectDisk(99)
StellarHalfmassRads = (il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snap_num,'SubhaloHalfmassRadType'))[:,4]
SubhaloPositions = il.groupcat.loadSubhalos('/Raid1/Illustris/TNG-100',snap_num,'SubhaloPos')

#Calculate each disk halo's A2 parameter
for haloID in ids:
    with h5py.File('/Raid0/zhouzb/cut_TNG/disk_%d/cutout_%d.hdf5'%(snap_num,haloID),'r') as cut:
        coor = np.array(cut['PartType4']['Coordinates'])
        mas = np.array(cut['PartType4']['Masses'])
        vel = np.array(cut['PartType4']['Velocities'])
        sf_time = np.array(cut['PartType4']['GFM_StellarFormationTime'])

    half_r = StellarHalfmassRads[haloID]
    halo_position = SubhaloPositions[haloID]

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

    np.save('/Raid0/zhouzb/TNG_a2/disk_%d/%d.npy'%(snap_num, haloID),np.vstack((A2_R,r_list[r_sort])))
    print('halo_%d datalen: %d'%(haloID,len(A2_R)))
    
print('All done')



