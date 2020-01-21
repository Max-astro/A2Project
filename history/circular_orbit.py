import numpy as np
import h5py
import illustris_python as il

basePath = '/Raid1/Illustris/Illustris-1/'

def groups(snap, haloID, fields, data):
    '''
    parameter:   
    fields : ['Group', 'Header', 'Offsets', 'Subhalo']
    data: halo's information that you want to extract
    '''
    with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_%d.0.hdf5'%snap,'r') as f:
        info = np.array(f['%s'%fields]['%s'%data])
    for n in range(1,8):
        with h5py.File('/Raid1/Illustris/Illustris-1/Groupcat/groups_%d.%d.hdf5'%(snap, n),'r') as f:
            info = np.concatenate((info, np.array(f['%s'%fields]['%s'%data])))
    if haloID == -1:
        return info
    else:
        return info[haloID]

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


#return subhalo's angular momentum
def Ang(haloID):
    coor = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
    mas = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Masses')
    vel = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Velocities')
    half_r = (groups(135, haloID, 'Subhalo','SubhaloHalfmassRadType'))[4]
    halo_position = groups(135, haloID, 'Subhalo', 'SubhaloPos')

    sf_time = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='GFM_StellarFormationTime')
    is_star = (sf_time>=0.0) # don't use wind particles


    r = coor - halo_position
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    vel = vel[(inside & is_star)]
    mas = mas[(inside & is_star)]
    coor = r[(inside & is_star)]

    for i in range(len(vel)):
        vel[i] = vel[i] * mas[i]

    J = np.cross(coor, vel)
    ang = J.sum(axis = 0)
    return ang




def sethalo(haloID):
    ''' 
    Set halo's angular momentum as it z axis, and set most bound particle as halo's center
    Return halo particles' new coordinate
    '''
    coor = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
    # vel = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Velocities')
    center = groups(135, haloID, 'Subhalo','SubhaloPos')
    with h5py.File('/Raid0/zhouzb/data/angular_momentum.hdf5','r') as f:
        tmp = np.array(f['Subhalo_JstarsInRad'])
    Jz = tmp[haloID]

    coor -= center

    #Align the Z axis with angular momentum vector
    #halo's rotation matrix
    rot = RotMatrix(Norm(Jz))

    for i in range(len(coor)):
        coor[i] = np.dot(rot, coor[i])    

    return coor



def Circ_Frac(haloID):
    cm = groups(135, haloID, 'Subhalo','SubhaloCM')
    # 6.67×10-11N·m²/kg²
    G = 6.67*10**(-11)

    coor = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
    mas = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Masses')
    vel = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Velocities')


    #Make a list to save ϵ = J / J(E) for each stellar particles
    Epsion = np.zeros(len(mas))

    #r10 is 10 times the stellar half-mass radius ['SubhaloHalfmassRadType']
    r10 = 10 * (groups(135, haloID, 'Subhalo','SubhaloHalfmassRadType'))[4]
    J_cir = (G * mas * r10)**0.5

    #Calculate each stellar particles' specific angular momentum Jz
    coor -= cm
    Jz = np.cross(coor,vel)  


    Epsion = Jz / J_cir

    CircAbove07Frac = len(Epsion[Epsion > 0.7]) / len(Epsion)
    return CircAbove07Frac
    
def MassTensorEigenVals(haloID):
    #Load subhalo particles' data
    #r = sethalo(haloID)
    r = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
    mas = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Masses')
    half_r = (groups(135, haloID, 'Subhalo','SubhaloHalfmassRadType'))[4]
    sf_time = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='GFM_StellarFormationTime')
    #cm = groups(135, haloID, 'Subhalo','SubhaloCM')
    # don't use wind particles
    r -= r[0]
    r = r[sf_time >= 0.0]
    mas = mas[sf_time >= 0.0]
    #calculate particle's distance from halo center
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*2)

    r = r[inside]

    Mi = ((mas[inside] * r[:,0]**2).sum())**0.5 / (mas[inside].sum())**0.5
    Mj = ((mas[inside] * r[:,1]**2).sum())**0.5 / (mas[inside].sum())**0.5
    Mk = ((mas[inside] * r[:,2]**2).sum())**0.5 / (mas[inside].sum())**0.5

    M = np.array([Mi, Mj, Mk])
    M.sort()
    return M

#halo's stellar particles mass center inside the 2 times of half mass radius
def stellar_cm(haloID):
    coor = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
    half_r = (groups(135, haloID, 'Subhalo','SubhaloHalfmassRadType'))[4]



# def half_r(haloID):
#     coor = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Coordinates')
#     mas = il.snapshot.loadSubhalo('/Raid1/Illustris/Illustris-1/', 135, haloID, partType=4, fields='Masses')  

#     #Most bound particle is the halo's center
#     coor -= coor[0]  
#     M = mas.sum()
#     r = ((coor**2).sum(1))**0.5
#     index = r.argsort()

#     M_sum = 0
#     for i in index:
#         M_sum += mas[i]
#         if M_sum > M*0.5:
#             return r[i]


