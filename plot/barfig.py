import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import illustris_python as il

# tng_barID = np.load('/Raid0/zhouzb/TNG/barredID_4WP_TNG.npy')
# il1_barID = np.load('/Raid0/zhouzb/il1_data/barredID_il1.npy')


tng_big = 172650
il1_big = 206714
tng_sml = 31347
il1_sml = 382757

# tng_sml = 428289
# tng_smID = 407021
# il1_smID = 449253



bbox = dict(boxstyle="round", fc="1")
# plt.annotate('text', (0.5, 0.5), xytext = (0.5, 0.5), textcoords = 'offset points', bbox = bbox)#, arrowprops = arrowprops)



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

def coorInfo(simu, snapnum, haloID):
    # data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
    data = np.load('f:/Linux/localRUN/fig1/%d.npy'%haloID, allow_pickle=True).item()
    coor= data['Coordinates']
    mas = data['Masses']
    vel = data['Velocities']
    sf_time = data['GFM_StellarFormationTime']
    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*6)

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
    return coor


tng_sMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMassType')[:, 4] / 0.6774
il1_sMass = il.func.loadSubhalos('il1', 135, 'SubhaloMassType')[:, 4] / 0.704

tng_half_r = (il.func.loadSubhalos('TNG', 99, 'SubhaloHalfmassRadType'))[:,4]
il1_half_r = (il.func.loadSubhalos('il1', 135, 'SubhaloHalfmassRadType'))[:, 4]



def fig4():
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    coor = coorInfo('TNG', 99, tng_big)
    ax[0][0].hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    ax[0][0].set_title('TNG-100 subhalo %d' % tng_big, fontsize=15)
    ax[0][0].annotate(r'$M_*\sim10^{11.5}$', (10, 25), xytext = (10, 26), textcoords = 'offset points', bbox = bbox, fontsize=14)
    ax[0][0].set_xlim(-40, 40)
    ax[0][0].set_ylim(-40, 40)

    coor = coorInfo('il1', 135, il1_big)
    ax[0][1].hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    ax[0][1].set_title('Illustris-1 subhalo %d' % il1_big, fontsize=15)
    ax[0][1].annotate(r'$M_* \sim 10^{11.5}$', (10, 25), xytext = (10, 26), textcoords = 'offset points', bbox = bbox, fontsize=14)
    ax[0][1].set_xlim(-40, 40)
    ax[0][1].set_ylim(-40, 40)

    coor = coorInfo('TNG', 99, tng_sml)
    ax[1][0].hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    ax[1][0].set_title('TNG-100 subhalo %d' % tng_sml, fontsize=15)
    ax[1][0].annotate(r'$M_* \sim 10^{10.8}$', (6.7, 17.6), xytext = (7, 15), textcoords = 'offset points', bbox = bbox, fontsize=14)
    ax[1][0].set_xlim(-25, 25)
    ax[1][0].set_ylim(-25, 25)

    coor = coorInfo('il1', 135, il1_sml)
    ax[1][1].hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    ax[1][1].set_title('Illustris-1 subhalo %d' % il1_sml, fontsize=15)
    ax[1][1].annotate(r'$M_* \sim 10^{10.9}$', (6.7, 17.6), xytext = (7, 15), textcoords = 'offset points', bbox = bbox, fontsize=14)
    ax[1][1].set_xlim(-25, 25)
    ax[1][1].set_ylim(-25, 25)

    fig.savefig('f:/Linux/local_result/fig1.pdf')


def plotBarFig(simu, snapnum, haloID):
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
    coor= data['Coordinates']
    mas = data['Masses']
    vel = data['Velocities']
    sf_time = data['GFM_StellarFormationTime']
    if simu == 'tng' or simu == 'TNG':
        half_r = tng_half_r[haloID]
        sMass = tng_sMass[haloID]
    else:
        half_r = il1_half_r[haloID]
        sMass = il1_sMass[haloID]

    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*6)

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
    #plot
    plt.rcParams['figure.figsize'] = (36, 36)
    plt.tick_params(labelsize=70)
    plt.hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')

    bbox = dict(boxstyle="round", fc="1")
    plt.annotate(r'$M_*$ = %.2f' % sMass, (22, 22), xytext=(22, 22), textcoords='offset points', bbox=bbox)
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)
    if simu == 'TNG' or simu == 'tng':
        savedir = '/Raid0/zhouzb/barfig/tng/%d.png' % haloID
    else:
        savedir = '/Raid0/zhouzb/barfig/il1/%d.png' % haloID
    plt.savefig(savedir)
    plt.close()

simu = 'TNG'
snapnum = 99
for haloID in barID:
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
    coor= data['Coordinates']
    mas = data['Masses']
    vel = data['Velocities']
    sf_time = data['GFM_StellarFormationTime']
    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*6)

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
    #plot
    plt.rcParams['figure.figsize'] = (36, 36)
    plt.tick_params(labelsize=20)
    plt.hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    plt.savefig('/Raid0/zhouzb/barfig/tng/%d.png'%haloID)
    plt.close()

simu = 'il1'
snapnum = 99
for haloID in il1_barID:
    data = il.func.loadgalaxy(simu, snapnum, haloID, partType=4, fields=['Coordinates', 'Masses', 'Velocities', 'GFM_StellarFormationTime'])
    coor= data['Coordinates']
    mas = data['Masses']
    vel = data['Velocities']
    sf_time = data['GFM_StellarFormationTime']
    half_r = (il.func.loadSubhalos(simu, snapnum, 'SubhaloHalfmassRadType'))[haloID,4]
    is_star = (sf_time>=0.0) # don't use wind particles

    r = coor - coor[0]
    dis = ((r**2).sum(1))**0.5
    inside = dis < (half_r*6)

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
    #plot
    plt.xticks()
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')
    plt.rcParams['figure.figsize'] = (36, 36)
    plt.hist2d(coor[::2, 0], coor[::2, 1], bins=500, norm=LogNorm(0.035), cmap='jet')
    plt.savefig('/Raid0/zhouzb/barfig/il1/%d.png'%haloID)
    plt.close()

