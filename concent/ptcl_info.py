from struct import *

def dumpToBin(path, N, Mass, rs, coor, vel):
    f = open(path, 'wb')
    f.write(pack('i', N))
    f.write(pack('d', Mass))
    f.write(pack('f', rs))
    for i in coor:
        for j in i:
            f.write(pack('f', j))
    for i in vel:
        for j in i:
            f.write(pack('f', j))
    f.close()

def dumpPPT(path, numbers, Idata, Fdata, coor):
    f = open(path, 'wb')
    #Ngr, GMostBound
    f.write(pack('i', numbers))
    for i in range(numbers):
        f.write(pack('i', Idata[i]))
    #
    for line in range(5):
        for i in range(numbers):
            f.write(pack('f', Fdata[line][i]))
    #GCenter Ngr*3 dimension
    for line in range(2): 
        for i in range(numbers):
            for p in range(3):
                f.write(pack('f', coor[line][i][p]))
    f.close()