#!/usr/bin/python3.11

B = [None for i in range(0, 3)]
nkpts, nbands = None, None

vbm, cbm = -1000, +1000
ef = +7.5

for fileind in range(0, 2):
    with open("bands-" + str(fileind), 'r') as infile:
        #read in B
        line = infile.readline().split()
        B[0] = [float(l) for l in line]
        line = infile.readline().split()
        B[1] = [float(l) for l in line]
        line = infile.readline().split()
        B[2] = [float(l) for l in line]

        #get counting information
        nkpts, nbands = [int(l) for l in infile.readline().split()]
   
        #now need to check the band extrema of each k-point 
        for i in range(0, nkpts):
            _, ka, kb, kc = [float(l) for l in infile.readline().split()] ##[0]=wgt

            for j in range(0, nbands):
                be = float(infile.readline())
                if(be < ef):
                    vbm = vbm if be < vbm else be
                if(be > ef):
                    cbm = cbm if be > cbm else be

        infile.close()

#Print the gap
print(cbm - vbm, "eV")
