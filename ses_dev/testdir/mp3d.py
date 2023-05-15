#!/usr/bin/python3.11

from numpy import empty
from matplotlib import pyplot as plt


def showmat(A):
    for i in range(0, A.shape[0]):
        for j in range(0, A.shape[1]):
            for k in range(0, A.shape[2]):
                if(A[i][j][k] == -10):
                    continue
                    ax.text(i, j, k, '.', ha="center", va="center")
                else:
                    ax.text(i, j, k, A[i][j][k], ha="center", va="center", fontsize=6)
    return

ndivs = [4, 5, 3]  ##[rows, cols]
inds = 0 ##slow index
indm = 1 ##medium index
indf = 2 ##fast index
lens = ndivs[inds] ##length of slow dim
lenm = ndivs[indm] ##length of medium dim
lenf = ndivs[indf] ##length of fast dim

A = empty(shape=(max(ndivs)+1, max(ndivs)+1, max(ndivs)+1), dtype=int)
A.fill(-10)

##slow, med, fast indices, followed by len(slow, med, fast) dims
##i have no reference for these - I managed to get them through "skill" (i.e. brute force, for the most part)
##note that these can be optimized a fair bit
def acc3_ee(s, m, f, ns, nm, nf):
    return s*nm*nf + m*nf + f
    return (s*nm + m)*nf + f
def acc3_eo(s, m, f, ns, nm, nf):
    return s*nm*nf + (m+1)*nf - (f+1)
    return (((s*nm + m + 1)*nf) - 1) - f
def acc3_oe(s, m, f, ns, nm, nf):
    return (s+1)*nm*nf -(m*nf + f + 1)
def acc3_oo(s, m, f, ns, nm, nf):
    return (s+1)*nm*nf -((m+1)*nf - (f+1) + 1)



#a path through the irreducible k point grid that minimizes the maximum distance between any two points
#i'm pretty sure that this is the optimal solution to the 3D trav. sals. prob. for a reciprocal space M.P.
#grid ... i'll take my comp sci degree now, thanks

#assume that lens is even
for s in range(0, lens//2):
    if(s%2 == 0): ##s is even
        for m in range(0, lenm):
            if(m%2 == 0): ##m is even
                for f in range(0, lenf):
                    A[s][m][f] = acc3_ee(s, m, f, lens, lenm, lenf)
            else:         ##m is odd
                for f in range(0, lenf):
                    A[s][m][f] = acc3_eo(s, m, f, lens, lenm, lenf)
    else:         ##s is odd
        for m in range(0, lenm):
            if(m%2 == 0): ##m is even
                for f in range(0, lenf):
                    A[s][m][f] = acc3_oe(s, m, f, lens, lenm, lenf)
            else:         ##m is odd
                for f in range(0, lenf):
                    A[s][m][f] = acc3_oo(s, m, f, lens, lenm, lenf)

#correction if lens is odd, assuming lenm is even
if(lens%2 != 0):
    s = lens//2
    if((lens-1)%4 != 0): ##true for 3, 7, 11, 15, ...
        for m in range(lenm-1, lenm//2, -1):
            if(m%2 == 0): ##m is even
                for f in range(0, lenf):
                    A[s][m][f] = acc3_oe(s, m, f, lens, lenm, lenf)
            else:         ##m is odd
                for f in range(0, lenf):
                    A[s][m][f] = acc3_oo(s, m, f, lens, lenm, lenf)
    else:                 ##true for 1, 5,  9, 13, ...
        for m in range(0, lenm//2):
            if(m%2 == 0): ##m is even
                for f in range(0, lenf):
                    A[s][m][f] = acc3_ee(s, m, f, lens, lenm, lenf)
            else:         ##m is odd
                for f in range(0, lenf):
                    A[s][m][f] = acc3_eo(s, m, f, lens, lenm, lenf)

#correction if lenm is odd, assuming lenf is even
if(lenm%2 != 0):
    m = lenm//2
    if((lens-1)%4 != 0):  ##true for 3, 7, 11, 15, ...
        if((lenm-1)%4 != 0): ##true for 3, 7, 11, 15, ...
            print("un")
            for f in range(0, lenf//2):
                A[s][m][f] = acc3_ee(s, m, f, lens, lenm, lenf) ##oo, ee
        else:                ##true for 1, 5,  9, 13, ...
            print("do")
            for f in range(lenf-1, lenf//2, -1):
                A[s][m][f] = acc3_oe(s, m, f, lens, lenm, lenf) ##oe, eo
    else:                 ##true for 1, 5,  9, 11, 15, ...
        if((lenm-1)%4 != 0): ##true for 3, 7, 11, 15, ...
            print("tr")
            for f in range(lenf-1, lenf//2, -1):
                A[s][m][f] = acc3_oe(s, m, f, lens, lenm, lenf) ##eo, oe
        else:                ##true for 1, 5,  9, 15, ...
            print("qa")
            for f in range(0, lenf//2):
                A[s][m][f] = acc3_ee(s, m, f, lens, lenm, lenf) ##ee, oo

#correction if lenf is odd (fills the gamma point)
if(lenf%2 != 0):
    f = lenf//2
    A[s][m][f] = lens*lenm*lenf//2










#plotting junk---------------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xticks([i for i in range(0, max(ndivs)+1)])
ax.set_yticks([i for i in range(0, max(ndivs)+1)])
ax.set_zticks([i for i in range(0, max(ndivs)+1)])
ax.set_xlim(0, max(ndivs)+1)
ax.set_ylim(0, max(ndivs)+1)
ax.set_zlim(0, max(ndivs)+1)
showmat(A)
plt.show()

