//Implementation of Monkhurst-Pack k-point grid generation

#include "def.h"

static forceinline ushort _ceil(const ushort a, const ushort b){
    return a/b + (a%b != (ushort)0);
}

/*
*  various ways of accessing the 1D representation of the MP grid
*  s, m, f are the slow, medium, and fast indices and n* are the lengths of
*  their dimensions
*/
static forceinline ushort _indee(const ushort s, const ushort m, const ushort f,
                                 const ushort nm, const ushort nf){
    return (s*nm + m)*nf + f;
}
static forceinline ushort _indeo(const ushort s, const ushort m, const ushort f,
                                 const ushort nm, const ushort nf){
    return (s*nm + m + (ushort)1)*nf - f - (ushort)1;
}
static forceinline ushort _indoe(const ushort s, const ushort m, const ushort f,
                                 const ushort nm, const ushort nf){
    return (s*nm + nm - m)*nf - f - (ushort)1;
}
static forceinline ushort _indoo(const ushort s, const ushort m, const ushort f,
                                 const ushort nm, const ushort nf){
    return (s*nm + nm - m - (ushort)1)*nf + f;
}

//Gives the fractional coordinates of monkhorst-pack k-points along b_i
//for an index i and number of subdivisions ni
//(altered from the original in that this works for 0-indexing)
static forceinline fp _mpfrac(const ushort i, const ushort ni){
    return (fp)((ushort)2*i + (ushort)1 - ni) / (fp)((ushort)2*ni);
}

/*
*  sets a m.p. grid path that minimizes the reciprocal space distance between
*  all k-points in the ibz.  i'm pretty sure this is the optimal solution to
*  the trav. sales. prob. for a 3d staggered grid w/ inversion symmetry ...
*  ... i'll take my comp sci degree now, thanks
*  as mentioned, this APPLIES INVERSION SYMMETRY!  
*/
//B is the recip lattice, sds are the number of subdivisions along b1, b2, b3,
//and kpts is the pre-allocated array of kpoints, size 
//sds[0]*sds[1]*sds[2]/2 + 1
void mpgridfi(const fp B[restrict 3][3], const ushort sds[restrict 3],
              kpt*restrict kpts){
    fp db2[3];
    ushort s, m, f, idxs, idxm, idxf, lens, lenm, lenf; 
    ushort map[3];
    fp w;
    ushort ind;

    //first, need to figure out which dimensions are slow, medium, and fast
    //(fastest index is the one with the smallest distance between kp values)
    db2[0] = (B[0][0]*B[0][0] + B[0][1]*B[0][1] + B[0][2]*B[0][2]) / (fp)sds[0];
    db2[1] = (B[1][0]*B[1][0] + B[1][1]*B[1][1] + B[1][2]*B[1][2]) / (fp)sds[1];
    db2[2] = (B[2][0]*B[2][0] + B[2][1]*B[2][1] + B[2][2]*B[2][2]) / (fp)sds[2];
    map[0] = (ushort)0; map[1] = (ushort)1; map[2] = (ushort)2;
    if(db2[0] > db2[1]){SWAP(fp, db2[0], db2[1]); SWAP(ushort, map[0], map[1]);}
    if(db2[1] > db2[2]){SWAP(fp, db2[1], db2[2]); SWAP(ushort, map[1], map[2]);}
    if(db2[0] > db2[1]){SWAP(fp, db2[0], db2[1]); SWAP(ushort, map[0], map[1]);}
    idxs = map[2]; idxm = map[1]; idxf = map[0];
    lens = sds[idxs]; lenm = sds[idxm]; lenf = sds[idxf];

    //most of the k-points have a weight = 2/nkpts b/c of inversion symmetry
    w = (fp)2.0 / (fp)(sds[0]*sds[1]*sds[2]);

    //now the fun part - actually making the k-points
#define kptcs kpts[ind].crds[idxs]
#define kptcm kpts[ind].crds[idxm]
#define kptcf kpts[ind].crds[idxf]
    ///slow index (assume length of slow index is even)
    for(s = (ushort)0; s < lens/(ushort)2; ++s){
        if(s%(ushort)2 == (ushort)0){ ///s even
            for(m = (ushort)0; m < lenm; ++m){
                if(m%(ushort)2 == (ushort)0){ ///m even
                    for(f = (ushort)0; f < lenf; ++f){
                        ind = _indee(s, m, f, lenm, lenf);
                        kptcs = _mpfrac(s, lens);
                        kptcm = _mpfrac(m, lenm);
                        kptcf = _mpfrac(f, lenf);
                        kpts[ind].wgt = w;
                    }
                }
                else{                         ///m odd
                    for(f = (ushort)0; f < lenf; ++f){
                        ind = _indeo(s, m, f, lenm, lenf);
                        kptcs = _mpfrac(s, lens);
                        kptcm = _mpfrac(m, lenm);
                        kptcf = _mpfrac(f, lenf);
                        kpts[ind].wgt = w;
                    }
                }
            }
        }
        else{                         ///s odd
            for(m = (ushort)0; m < lenm; ++m){
                if(m%(ushort)2 == (ushort)0){ ///m even
                    for(f = (ushort)0; f < lenf; ++f){
                        ind = _indoe(s, m, f, lenm, lenf);
                        kptcs = _mpfrac(s, lens);
                        kptcm = _mpfrac(m, lenm);
                        kptcf = _mpfrac(f, lenf);
                        kpts[ind].wgt = w;
                    }
                }
                else{                         ///m odd
                    for(f = (ushort)0; f < lenf; ++f){
                        ind = _indoo(s, m, f, lenm, lenf);
                        kptcs = _mpfrac(s, lens);
                        kptcm = _mpfrac(m, lenm);
                        kptcf = _mpfrac(f, lenf);
                        kpts[ind].wgt = w;
                    }
                }
            }
        }
    }

    if(lens%(ushort)2 == (ushort)0) goto end;

    ///medium index (correction if slow index length is odd, assume length of
    ///medium index is even)
    s = lens/(ushort)2; ///s is positioned at b_slow = 0 (aligned w/ gamma)
    if((lens-(ushort)1)%(ushort)4 == 0){ ////case lens = 1, 5, 9, ...
        for(m = (ushort)0; m < lenm/(ushort)2; ++m){
            if(m%(ushort)2 == 0){ ////m even
                for(f = (ushort)0; f < lenf; ++f){
                    ind = _indee(s, m, f, lenm, lenf);
                    kptcs = ZERO; ////_mpfrac(s, lens);
                    kptcm = _mpfrac(m, lenm);
                    kptcf = _mpfrac(f, lenf);
                    kpts[ind].wgt = w;
                }
            }
            else{                 ////m odd
                for(f = (ushort)0; f < lenf; ++f){
                    ind = _indeo(s, m, f, lenm, lenf);
                    kptcs = ZERO; ////_mpfrac(s, lens);
                    kptcm = _mpfrac(m, lenm);
                    kptcf = _mpfrac(f, lenf);
                    kpts[ind].wgt = w;
                }
            }
        }
    }
    else{                                ////case lens = 3, 7, 11, ...
        for(m = _ceil(lenm, (ushort)2); m < lenm; ++m){
            if(m%(ushort)2 == 0){ ////m even
                for(f = (ushort)0; f < lenf; ++f){
                    ind = _indoe(s, m, f, lenm, lenf);
                    kptcs = ZERO; ////_mpfrac(s, lens);
                    kptcm = _mpfrac(m, lenm);
                    kptcf = _mpfrac(f, lenf);
                    kpts[ind].wgt = w;
                }
            }
            else{                 ////m odd
                for(f = (ushort)0; f < lenf; ++f){
                    ind = _indoo(s, m, f, lenm, lenf);
                    kptcs = ZERO; ////_mpfrac(s, lens);
                    kptcm = _mpfrac(m, lenm);
                    kptcf = _mpfrac(f, lenf);
                    kpts[ind].wgt = w;
                }
            }
        }
    }

    if(lenm%(ushort)2 == (ushort)0) goto end;
    
    ///fast index (correction if medium index length is odd, assume length of
    ///fast index is even)
    m = lenm/(ushort)2;
    if((lens-(ushort)1)%(ushort)4 == (ushort)0){ ////case lens = 1, 5, 9, ...
        if((lenm-(ushort)1)%(ushort)4 == (ushort)0){ ////case lenm = 1, 5, 9, 
            for(f = (ushort)0; f < lenf/(ushort)2; ++f){
                ind = _indee(s, m, f, lenm, lenf);
                kptcs = ZERO; ////_mpfrac(s, lens);
                kptcm = ZERO; ////_mpfrac(m, lenm);
                kptcf = _mpfrac(f, lenf);
                kpts[ind].wgt = w;
            }
        }
        else{                                        ////case lenm = 3, 7, 11,
            for(f = _ceil(lenf, (ushort)2); f < lenf; ++f){
                ind = _indoe(s, m, f, lenm, lenf);
                kptcs = ZERO; ////_mpfrac(s, lens);
                kptcm = ZERO; ////_mpfrac(m, lenm);
                kptcf = _mpfrac(f, lenf);
                kpts[ind].wgt = w;
            }
        }
    }  
    else{                                        ////case lens = 3, 7, 11, ...
        if((lenm-(ushort)1)%(ushort)4 == (ushort)0){ ////case lenm = 1, 5, 9, 
            for(f = _ceil(lenf, (ushort)2); f < lenf; ++f){
                ind = _indoe(s, m, f, lenm, lenf);
                kptcs = ZERO; ////_mpfrac(s, lens);
                kptcm = ZERO; ////_mpfrac(m, lenm);
                kptcf = _mpfrac(f, lenf);
                kpts[ind].wgt = w;
            }
        }
        else{                                        ////case lenm = 3, 7, 11,
            for(f = (ushort)0; f < lenf/(ushort)2; ++f){
                ind = _indee(s, m, f, lenm, lenf);
                kptcs = ZERO; ////_mpfrac(s, lens);
                kptcm = ZERO; ////_mpfrac(m, lenm);
                kptcf = _mpfrac(f, lenf);
                kpts[ind].wgt = w;
            }
        }
    }

    if(lenf%(ushort)2 == 0) goto end;

    ///finally, the gamma point, which we miss if lens, lenm, lenf are all odd
    ind = lens*lenm*lenf / (ushort)2;
    kptcs = ZERO; ////_mpfrac(s, lens);
    kptcm = ZERO; ////_mpfrac(m, lenm);
    kptcf = ZERO; ////_mpfrac(f, lenf); ////if uncommented (why??), set f first
    kpts[ind].wgt = HALF*w; ///gamma has no symmetry, so it's weight is 1/nkpts
#undef kptcs 
#undef kptcm
#undef kptcf

    end: //wow vic, this seems pointless.  why not just return instead of using
         //a label that only returns?  its because if some poor soul wants to
         //do further symmetry operations, they'll probably have to do some 
         //clean-up here
    return;
}





