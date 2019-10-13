import os

snapList = [99, 92, 89, 82, 76, 73, 71, 69]

fname_pre = "http://www.illustris-project.org/api/Illustris-1/files/groupcat-"

for snap in snapList:
    for i in range(8):
        url = "http://www.tng-project.org/api/Illustris-1/files/groupcat-%d.%d.hdf5"%(snap,i)
        os.system("wget -nc --content-disposition --header=\"API-Key: 6be79adf6f05af759ee1ada46249f9b4\" %s"%url)
        print("snap_%d part %d done"%(snap,i))
