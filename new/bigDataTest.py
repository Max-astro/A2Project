import numpy as np
import sys
import h5py
import time
import illustris_python as il

f = h5py.File('/Raid1/Illustris/Illustris-1/snap_ics.hdf5', 'r')
particles = il.func.loadgalaxy('il1', 135, 344358, partType=1, fields='ParticleIDs')
#Pick 1000 particles ramdomly
pickInd = np.random.randint(0, len(particles), 1000)
seleID =particles[pickInd]

ids = np.array(f['PartType1']['ParticleIDs'])

#Find one
load_start = time.time()
f1=np.where(ids[:] == seleID[0])[0]
load_end = time.time()
print('Find one, Time use: ', load_end - load_start)

#find 10
load_start = time.time()
f2=[]
for i in seleID[:10]:
    f2.append(np.where(ids[:] == i)[0][0])
load_end = time.time()
print('Find 10p, Time use: ', load_end - load_start)
print(' ')
print(f1)
print(f2)
# #find 100
# load_start = time.time()
# f3=np.where(ids[:] == seleID[:100])
# load_end = time.time()
# print('Find 100p, Time use: ', load_end - load_start)

# #find all 1000 particle
# load_start = time.time()
# f4=np.where(ids == seleID[:100])
# load_end = time.time()
# print('Find a total subhalo(1000p), Time use: ', load_end - load_start)
# print('')

