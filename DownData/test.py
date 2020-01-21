import h5py
import numpu as np
for snap in [50, 59, 67]:
    for i in range(448):
        try:
            with h5py.File('/home/snap_%d/snap_0%d.%d.hdf5'%(snap, snap, i), 'r') as f:
                a = np.array(f['PartType4']['Masses'])[0]
        except:
            print('snap_%d part_%d broken'%(snap, i))


import numpu as np
import os
nums = np.load('/home/nums.npy')
for i in nums:
    if os.path.isfile('/home/snap_67/snap_067.%d.hdf5'%i):
        os.system("rm -f /home/snap_67/snap_067.%d.hdf5"%i) 