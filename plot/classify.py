import numpy as np
import os


path = 'f:/Linux/scp data/selected/'

flyby = []
for files in os.listdir(path + 'flyby'):
    flyby.append(int(files.split('.')[0]))

secular = []
for files in os.listdir(path + 'secular'):
    secular.append(int(files.split('.')[0]))

merge = []
for files in os.listdir(path + 'merge'):
    merge.append(int(files.split('.')[0]))

np.savez('f:/Linux/localRUN/haloID_Classify', flyby = flyby, secular = secular, merge = merge)
