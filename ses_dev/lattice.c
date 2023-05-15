//Implementation of functions for working on the real/recip space lattice
#include <string.h>
#include "def.h"

//Change of basis x = Ay, stores result in x
//A should be stored with its mathamatical columns as rows
static forceinline void _cob(const fp A[restrict 3][3],
                             fp x[restrict 3], const fp y[restrict 3]){
    x[0] = A[0][0]*y[0] + A[1][0]*y[1] + A[2][0]*y[2];
    x[1] = A[0][1]*y[0] + A[1][1]*y[1] + A[2][1]*y[2];
    x[2] = A[0][2]*y[0] + A[1][2]*y[1] + A[2][2]*y[2];
}

//Sets cartesian coordinates of entire lattice
void tocrdsc(const lat lattice){
    uchar i;
    ushort j, acc;

    for(i = (uchar)0, acc = (ushort)0; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < (ushort)(lattice.speccounts[i]); ++j, ++acc){
            _cob(lattice.A, lattice.sites[acc].crdsc, 
                 lattice.sites[acc].crdsf);
        }
    }
}

void tocrdsf(const fp Ai[restrict 3][3], const lat lattice){
    uchar i;
    ushort j, acc;

    for(i = (uchar)0, acc = (ushort)0; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < (ushort)(lattice.speccounts[i]); ++j, ++acc){
            _cob(Ai, lattice.sites[acc].crdsf, lattice.sites[acc].crdsc);
        }
    }
}

//moves FRACTIONAL COORDS ONLY back into the range [0, 1).
//optionally shifts the images as well
void backtocell(const lat lattice, const uchar imgs){
    uchar i; ushort j, k, acc;
#define atom lattice.sites[acc] 
#define img lattice.sites[acc].imgs[k]

    //don't mess with imgs, just do main atoms
    if(!imgs){
        for(i = (uchar)0, acc = 0u; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
            while(atom.crdsf[0] < ZERO) atom.crdsf[0] += ONE;
            while(atom.crdsf[0] >= ONE) atom.crdsf[0] -= ONE;
            while(atom.crdsf[1] < ZERO) atom.crdsf[1] += ONE;
            while(atom.crdsf[1] >= ONE) atom.crdsf[1] -= ONE;
            while(atom.crdsf[2] < ZERO) atom.crdsf[2] += ONE;
            while(atom.crdsf[2] >= ONE) atom.crdsf[2] -= ONE;
        }
        }   
    }
    //also adjust the periodic images
    else{
        for(i = (uchar)0, acc = 0u; i < lattice.nspecs; ++i){
        for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
            ///main
            while(atom.crdsf[0] < ZERO) atom.crdsf[0] += ONE;
            while(atom.crdsf[0] >= ONE) atom.crdsf[0] -= ONE;
            while(atom.crdsf[1] < ZERO) atom.crdsf[1] += ONE;
            while(atom.crdsf[1] >= ONE) atom.crdsf[1] -= ONE;
            while(atom.crdsf[2] < ZERO) atom.crdsf[2] += ONE;
            while(atom.crdsf[2] >= ONE) atom.crdsf[2] -= ONE;
            ///imgs
            for(k = (ushort)0; k < atom.nimgs; ++k){
                while(img->crdsf[0] < ZERO) img->crdsf[0] += ONE;
                while(img->crdsf[0] >= ONE) img->crdsf[0] -= ONE;
                while(img->crdsf[1] < ZERO) img->crdsf[1] += ONE;
                while(img->crdsf[1] >= ONE) img->crdsf[1] -= ONE;
                while(img->crdsf[2] < ZERO) img->crdsf[2] += ONE;
                while(img->crdsf[2] >= ONE) img->crdsf[2] -= ONE;
            }
        }
        }   
    }
#undef atom
#undef img
}

//dot product
static forceinline fp _dot(const fp a[restrict 3], const fp b[restrict 3]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//cross product
static forceinline void _crs(const fp a[restrict 3], const fp b[restrict 3],
                             fp res[restrict 3]){
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2]; ///not a mistake, just implicit mul by -1
    res[2] = a[0]*b[1] - a[1]*b[0];
}

//Changes A -> B w/ A the row vectors {a1, a2, a3} -> {b1, b2, b3}
void torecipbasis(fp A[restrict 3][3]){
    fp scv;
    fp tmp[3][3] = {{A[0][0], A[0][1], A[0][2]},
                    {A[1][0], A[1][1], A[1][2]},
                    {A[2][0], A[2][1], A[2][2]}};
    
    //this algo is really just setting bi = 2pi/V * (aj cross ak) w/ V = vol(A)
    ///calc scaling factor
    _crs(tmp[1], tmp[2], A[0]);
    scv = TWOPI / _dot(tmp[0], A[0]);
    ///cross products
    _crs(tmp[1], tmp[2], A[0]);
    _crs(tmp[2], tmp[0], A[1]);
    _crs(tmp[0], tmp[1], A[2]);
    ///apply scaling factor
    A[0][0] *= scv; A[0][1] *= scv; A[0][2] *= scv;
    A[1][0] *= scv; A[1][1] *= scv; A[1][2] *= scv;
    A[2][0] *= scv; A[2][1] *= scv; A[2][2] *= scv;
}


//changes A -> A^-1 assuming A is a 3x3 matrix
void toinverse(fp A[restrict 3][3]){
    fp det33, B[3][3];
#define det22(tl, tr, bl, br) (tl*br - tr*bl)

    //get inv(a) by dividing the co-factor matrix by determinant
    det33 = A[0][0]*det22(A[1][1], A[1][2], A[2][1], A[2][2]) -
            A[0][1]*det22(A[1][0], A[1][2], A[2][0], A[2][2]) +
            A[0][2]*det22(A[1][0], A[1][1], A[2][0], A[2][1]);
    det33 = ONE / det33; ///finite for physically possible lattices

    //note that the transpose is implicit
    B[0][0] = +det22(A[1][1], A[2][1], A[1][2], A[2][2])*det33;
    B[0][1] = -det22(A[0][1], A[2][1], A[0][2], A[2][2])*det33;
    B[0][2] = +det22(A[0][1], A[1][1], A[0][2], A[1][2])*det33;
    B[1][0] = -det22(A[1][0], A[2][0], A[1][2], A[2][2])*det33;
    B[1][1] = +det22(A[0][0], A[2][0], A[0][2], A[2][2])*det33;
    B[1][2] = -det22(A[0][0], A[1][0], A[0][2], A[1][2])*det33;
    B[2][0] = +det22(A[1][0], A[2][0], A[1][1], A[2][1])*det33;    
    B[2][1] = -det22(A[0][0], A[2][0], A[0][1], A[2][1])*det33;
    B[2][2] = +det22(A[0][0], A[1][0], A[0][1], A[1][1])*det33;

    memcpy(A, B, 9u*sizeof(fp));

#undef det22
}



