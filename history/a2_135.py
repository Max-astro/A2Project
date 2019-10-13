import numpy as np
import h5py
import os
import illustris_python as il

#Vector Normalized
def Norm(vect):
    '''Input a vector, normalize it in-place'''
    mol=(vect[0]**2+vect[1]**2+vect[2]**2)**0.5
    for i in range(3):
        vect[i]=vect[i]/mol

    return vect

#Normalized particles' coordinates in halo by minus halo center of mass
def HaloNorm(cm,coor):
    '''Input: (coordinates of halo center , particle's coordinates), alter the particle's coordinates in-place'''
    for i in range(3):
        coor[i] = coor[i] - cm[i]

    mol=(coor[0]**2+coor[1]**2+coor[2]**2)**0.5
    for i in range(3):
        coor[i]=coor[i]/mol

    return coor

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

#Get subhalo's star particle number from groups
#Subhalolen: [halo0, halo0+1, halo0+1+2,....,Sigma(halon)]
def groups_135(data):
    '''
    parameter: data= 'len' or 'cm'. 'len' return SubhaloID corresponded stellar numbers, shape:(halonumber+1,). 'cm' return Subhalo's center of mass's 3D-coordinates, shape:(halonumber, 3).
    '''
    if data == 'len' :
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.0.hdf5','r') as f:
            Subhalolen = np.array(f['Offsets']['Subhalo_SnapByType'][:,4])
        for n in range(1,8):
            with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.%d.hdf5'%n,'r') as f:
                Subhalolen = np.concatenate((Subhalolen, np.array(f['Offsets']['Subhalo_SnapByType'][:,4])))
        return Subhalolen

    elif data == 'cm':
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.0.hdf5','r') as f:
            SubhaloCM_list = np.array(f['Subhalo']['SubhaloCM'])
        for n in range(1,8):
            with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_135.%d.hdf5'%n,'r') as f:
                SubhaloCM_list = np.concatenate((SubhaloCM_list, np.array(f['Subhalo']['SubhaloCM'])))
        return SubhaloCM_list

#Get subhalo's star particles
#Parameter: (halo's ID, This halo's particle number in group catalog)
#return array.shape: (This halo's stellar particle number, 4) .  which [0] is the stellar particle's mass, [1,2,3] are particles's 3D-coordinates
def gethalo(haloID, Subhalolen):
    '''Input:(haloID, a list about all halo's particles number), return the halo's stellar particles's mass and coordinates, shape(4, stellar number)'''
    start = Subhalolen[haloID]
    end = Subhalolen[haloID + 1]
    #haloinfo.shape = (This halo's stellar particle number, 4) .  which [0] is the stellar particle's mass, [1,2,3] are particles's 3D-coordinates

    part = 0
    while(end > 0):
        with h5py.File('/Raid1/Illustris/Illustris-1/Snapshot/snap_135/snap_135.%d.hdf5'%part,'r') as H5File:
            #Partlen is a number which counts this part of file's stellar particles
            Partlen = len(np.array(H5File['PartType4']['Masses']))
            if start - Partlen < 0:

                #start/end point in the same file
                stars_mass = np.array(H5File['PartType4']['Masses'][start : end]) if end <= Partlen else np.array(H5File['PartType4']['Masses'][start :])
                stars_coor = np.array(H5File['PartType4']['Coordinates'][start : end,:]) if end <= Partlen else np.array(H5File['PartType4']['Coordinates'][start :,:])
                end -= Partlen

                #Need to load more than 1 part of snapshot file
                while(end > 0):
                    part += 1
                    with h5py.File('/Raid1/Illustris/Illustris-1/Snapshot/snap_135/snap_135.%d.hdf5'%part,'r') as H5File:
                        Partlen = len(np.array(H5File['PartType4']['Masses']))
                        stars_mass = np.concatenate((stars_mass, np.array(H5File['PartType4']['Masses'][: end]))) if end <= Partlen else np.concatenate((stars_mass, np.array(H5File['PartType4']['Masses'][:])))
                        stars_coor = np.concatenate((stars_coor, np.array(H5File['PartType4']['Coordinates'][: end,:]))) if end <= Partlen else np.concatenate((stars_coor, np.array(H5File['PartType4']['Coordinates'][:,:])))
                        end -= Partlen

                return np.vstack((stars_mass, stars_coor.T))

            else:
                start -= Partlen
                end -= Partlen
                part += 1



#Select disk galasy SubhaloID
#return a list of disk galaxies's ID
#def finddisk_135():
    # '''Return a list of all Disk galaxies in this snapshot'''
    # with h5py.File('/Raid0/zhouzb/stellar_circs.hdf5','r') as cir:
    #     return np.array(cir['Snapshot_135']['SubfindID'])

#load subhalo_0's data:
# mass center coordinate : halo_cm
# stellar numbers : halo_starlen
# stellars' coordinates : halo_starcoor
#

#DiskHalo_ID = finddisk_135()
#load 50 Subhalo's stellar particles
#DiskHalo_ID = DiskHalo_ID[0]
#load group catalog information in groups_135.hdf5 files
halolen_list = groups_135(data = 'len')
haloCM_list = groups_135(data = 'cm')

new = np.load('/Raid0/zhouzb/data/snap_135_diskID.npy','r')
for haloID in new:
    # Load halo stellar mass and stellar particles' 3D-coordinates
    halostellar = gethalo(haloID, halolen_list)
    #select stellar particles > 40000

    # Load stellar numbers
    halo_particlen = halolen_list[haloID + 1] - halolen_list[haloID]
    # Put all halo information in appropriate list
    halo_SPMass = halostellar[0,:]
    halo_starcoors = halostellar[1:,:].T
    halo_cm = haloCM_list[haloID,:]
    if len(halo_starcoors) < 40000:
        continue

    #Load halo_0's angular mommentum as Z axis
    angfile = h5py.File('/Raid0/zhouzb/data/angular_momentum.hdf5','r')
    halo_J=np.array(angfile['Subhalo_Jstars'])
    #halo_0's angular momentum
    vect = halo_J[haloID]


    #Normalize
    vect=Norm(vect)
    #halo_0's rotation matrix
    halo_rot = RotMatrix(vect)

    #set halocm as (0,0,0) and normlized star particles' coordinates
    halo_starcoors = halo_starcoors - halo_cm
    #Set all particles in halo_0 into new axis
    for i in range(halo_particlen):

        #Rotated star particles into new coordinate, which halo_J is the Z-axis
        halo_starcoors[i] = np.dot(halo_rot, halo_starcoors[i])


    #Calculate distance r
    #r_list is a list about distance of star particles to (0,0,0)
    r_list = np.zeros(halo_particlen)
    for i in range(halo_particlen):
        r_list[i] = (halo_starcoors[i][0]**2 + halo_starcoors[i][1]**2)**0.5
    #r_sort : Index about r_list from small to big
    r_sort = r_list.argsort()

    #a0 = Sigma(M[i])
    a0 = np.zeros(len(r_sort))
    ptr = 0
    for i in r_sort:
        if ptr == 0:
            a0[ptr] = halo_SPMass[i]
        else:
            a0[ptr] = halo_SPMass[i] + a0[ptr-1]
        ptr += 1

    #Position angle: Theta = arctan(y/x)
    #Creat a list of a_m(R)
    a2_R = np.zeros(len(r_sort))
    b2_R = np.zeros(len(r_sort))
    list_i = 0
    for numb in r_sort:

        Theta = np.arctan(halo_starcoors[numb][0]/halo_starcoors[numb][1])
        #a_i = M[numb] * cos( m * Theta[numb] )  , m=2
        a_i = halo_SPMass[numb] * np.cos(2*Theta)
        #b_i = M[numb] * sin( m * Theta[numb] ) , m=2
        b_i = halo_SPMass[numb] * np.sin(2*Theta)

    #a2_R[i] = a_i.sum(:i)
        if list_i > 0:
            a2_R[list_i] = a_i + a2_R[list_i - 1]
            b2_R[list_i] = b_i + b2_R[list_i - 1]
        else:
            a2_R[list_i] = a_i
            b2_R[list_i] = b_i

        list_i += 1

    #A2_R = (a2**2 + b2**2)**0.5 / a0
    A2_R = np.zeros(len(a2_R))
    for i in range(len(a2_R)):
        A2_R[i] = ((a2_R[i]**2 + b2_R[i]**2)**0.5 / a0[i])

    #Output haloID, A2_R, r_list(r = x**2+y**2)
#    with h5py.File('/Raid0/zhouzb/A2_test50/%d.hdf5'%haloID,'w') as op:
#        op.create_dataset('A2_R', data=A2_R)
#        op.create_dataset('r_list', data=r_list)

    np.save('/Raid0/zhouzb/a2_135/%d.npy'%haloID,np.vstack((A2_R,r_list[r_sort])))
    print('halo_%d datalen: %d'%(haloID,len(A2_R)))

print('All clear')
