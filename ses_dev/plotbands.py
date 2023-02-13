#!/usr/bin/python3.11
from matplotlib import pyplot as plt


B = [None for i in range(0, 3)]
nkpts, nbands = None, None
x, ys = [], []

with open("bands", 'r') as infile:
    #read in B
    line = infile.readline().split()
    B[0] = [float(l) for l in line]
    line = infile.readline().split()
    B[1] = [float(l) for l in line]
    line = infile.readline().split()
    B[2] = [float(l) for l in line]

    #get counting information
    nkpts, nbands = [int(l) for l in infile.readline().split()]
    
    #now read in the k-points and make nbands arrays of band energies
    #also need to keep track of distance between k-points
    cx = 0.0
    lx, ly, lz = 0.0, 0.0, 0.0
    x = [None for i in range(0, nkpts)]
    y = [[None for i in range(0, nkpts)] for j in range(0, nbands)]
    for i in range(0, nkpts):
        _, ka, kb, kc = [float(l) for l in infile.readline().split()] ##[0]=wgt
        ##transform k from direct to recip space
        kx = ka*B[0][0] + kb*B[1][0] + kc*B[2][0]
        ky = ka*B[0][1] + kb*B[1][1] + kc*B[2][1]
        kz = ka*B[0][2] + kb*B[1][2] + kc*B[2][2]

        ##set this x point as the distance between current and last k point
        cx += (  (kx-lx)**2 + (ky-ly)**2 + (kz-lz)**2  )**(1/2)  
        x[i] = cx  

        ##set the current coords to the last coords
        lx, ly, lz = kx, ky, kz

        ##now, move on to the actual bands
        for j in range(0, nbands):
            y[j][i] = float(infile.readline())


    infile.close()

#Plot the band structure: one plot for each band
for i in range(0, nbands):
    plt.plot(x, [y_ - 6.708 for y_ in y[i]], color="#0058fd")
#    plt.plot(x, y[i], color="#eb7c00")

plt.show()




