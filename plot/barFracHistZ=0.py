import numpy as np
import h5py
import json
import sys
import seaborn as sns
sys.path.append('F:\Linux')
import illustris_python as il
from datetime import datetime
sys.path.append(r"C:/Users/qq651/OneDrive/Codes/A2project/plot/")
from Tools import *

#Galaxies morphlogy info
tng2il1 = np.load('F:/Linux/localRUN/Match/tng2il1_allsub.npy',allow_pickle=1).item()
bar2bar = np.load('F:/npy/bar2bar.npy',allow_pickle=1).item()
bar2disk = np.load('f:/npy/bar2no.npy',allow_pickle=1).item()

il1_barID = np.load('f:/Linux/localRUN/barredID_il1.npy')
il1_diskID = np.load('f:/Linux/localRUN/diskID_il1.npy')

tng_barID = np.load('f:/Linux/localRUN/barredID_4WP_TNG.npy')
tng_diskID = np.load('f:/Linux/localRUN/diskID_4WP.npy')

tng_barID_total = np.load('f:/Linux/localRUN/barredID_TNG.npy')
tng_diskID_total = np.load('f:/Linux/localRUN/diskID_TNG.npy')

#halo info
il1_sMass = il.func.loadSubhalos('il1', 135, 'SubhaloMassType')[:, 4] / 0.704
il1_sMass = np.log10(il1_sMass * 10 ** 10)
il1_sMass[np.isinf(il1_sMass)] = 0
tng_sMass = il.func.loadSubhalos('TNG', 99, 'SubhaloMassType')[:, 4] / 0.6774
tng_sMass = np.log10(tng_sMass * 10 ** 10)
tng_sMass[np.isinf(tng_sMass)] = 0


#
bins = np.linspace(10, 12, 25)
barhigh = HistValAndBin(il1_sMass[il1_diskID], bins, more=1)
sns.barplot(bins, barhigh)



